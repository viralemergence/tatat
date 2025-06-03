from argparse import ArgumentParser
from collections import defaultdict
from contextlib import contextmanager
from os import chdir, environ, getcwd
from pathlib import Path
from shutil import rmtree
import sqlite3
import subprocess
from tempfile import mkdtemp
from typing import Any, Iterator, Union

@contextmanager
def temporarily_change_working_directory(new_directory: Path):
    starting_directory = getcwd()
    try:
        chdir(new_directory)
        yield
    finally:
        chdir(starting_directory)

class InputFastaManager:
    def __init__(self, assembly_fasta: Path, outdir: Path, sqlite_db: Path, transcriptome: str, prefix_column: str) -> None:
        self.assembly_fasta = assembly_fasta
        self.outdir = outdir
        self.sqlite_db = sqlite_db
        self.transcriptome = transcriptome
        self.prefix_column = prefix_column

    def run(self) -> Path:
        transcriptome_filtered_transcript_ids = self.extract_filtered_transcript_ids(self.sqlite_db, self.transcriptome)
        transcript_prefix_mapping = self.extract_transcript_prefix_mapping(self.prefix_column, self.sqlite_db)

        tmp_prefixed_fasta = self.write_temporary_prefixed_fasta(self.assembly_fasta, self.outdir,
                                                                 transcriptome_filtered_transcript_ids, transcript_prefix_mapping)
        return tmp_prefixed_fasta

    @staticmethod
    def extract_filtered_transcript_ids(sqlite_db: Path, transcriptome: str) -> set[int]:
        print(f"\nExtracting transcript ids that belong to transcriptome: {transcriptome}")
        with sqlite3.connect(sqlite_db, timeout=600) as connection:
            cursor = connection.cursor()
            sql_query = ("SELECT t.uid "
                         "FROM transcripts t "
                         "LEFT OUTER JOIN samples s ON t.sample_uid = s.uid "
                         f"WHERE s.transcriptome = '{transcriptome}'")
            cursor.execute(sql_query)
            return {row[0] for row in cursor.fetchall()}

    @staticmethod
    def extract_transcript_prefix_mapping(prefix_column: str, sqlite_db: Path) -> dict[str]:
        print("Extracting transcript id to prefix mapping")
        with sqlite3.connect(sqlite_db, timeout=600) as connection:
            cursor = connection.cursor()
            sql_query = f"SELECT uid,{prefix_column} FROM transcripts"
            cursor.execute(sql_query)
            return {row[0]: row[1] for row in cursor.fetchall()}
    
    @classmethod
    def write_temporary_prefixed_fasta(cls, assembly_fasta: Path, outdir: Path,
                                       filtered_transcript_ids: set[int], transcript_prefix_mapping: dict[str]) -> Path:
        print("Writing temporary prefixed fasta")
        temp_dir = mkdtemp(dir=outdir)
        outfile = Path(temp_dir) / assembly_fasta.name
        with outfile.open("w") as outhandle:
            for fasta_seq in cls.fasta_chunker(assembly_fasta):
                transcript_id = int(fasta_seq[0][1:])
                if transcript_id not in filtered_transcript_ids:
                    continue

                prefix = transcript_prefix_mapping[transcript_id]
                new_header = ">" + "_".join([prefix, "prefix", str(transcript_id)])
                outhandle.write(f"{new_header}\n")
                for line in fasta_seq[1:]:
                    outhandle.write(f"{line}\n")
        return outfile

    @staticmethod
    def fasta_chunker(fasta_path: Path) -> Iterator[list[str]]:
        fasta_seq = []
        first_chunk = True
        with fasta_path.open() as inhandle:
            for line in inhandle:
                line = line.strip()
                if not line.startswith(">"):
                    fasta_seq.append(line)
                else:
                    if first_chunk:
                        fasta_seq.append(line)
                        first_chunk = False
                        continue
                    yield fasta_seq
                    fasta_seq = [line]
            if fasta_seq:
                yield fasta_seq

class EvigeneManager:
    def __init__(self, assembly_fasta: Path, outdir: Path, cpus: int, memory: int) -> None:
        self.assembly_fasta = assembly_fasta
        self.outdir = outdir
        self.cpus = cpus
        self.memory = memory

    # Evigene assembly classifier
    def run_assembly_classifier(self, phetero: Union[None, int], minaa: Union[None, int]) -> None:
        soft_link_path = self.set_soft_link_path(self.outdir, self.assembly_fasta)
        self.make_softlink(self.assembly_fasta, soft_link_path)

        with temporarily_change_working_directory(self.outdir):
            self.run_evigene(soft_link_path, self.cpus, self.memory, phetero, minaa)

    @staticmethod
    def set_soft_link_path(outdir: Path, assembly_fasta_path: Path) -> Path:
        return outdir / assembly_fasta_path.name

    @staticmethod
    def make_softlink(original_path: Path, soft_link_path: Path) -> None:
        ln_command = ["ln", "-s", f"{original_path}", f"{soft_link_path}"]
        p = subprocess.Popen(ln_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        while p.poll() is None and (line := p.stdout.readline()) != "":
            print(line.strip())
        p.wait()

        if p.poll() != 0:
            raise Exception(p.stderr.readlines())

    @staticmethod
    def run_evigene(soft_link_path: Path, cpus: int, memory: int, phetero: Union[None, int], minaa: Union[None, int]) -> None:
        environment_variables = environ.copy()
        evigene_path = environment_variables["EVIGENE"]

        evigene_command = [f"{evigene_path}/scripts/prot/tr2aacds.pl", "-NCPU", f"{cpus}", "-MAXMEM", f"{memory}",
                           "-log", "-cdna", f"{soft_link_path}"]
        if phetero:
            evigene_command.extend([f"-pHetero={phetero}"])
        if minaa:
            evigene_command.extend([f"-MINAA={minaa}"])

        p = subprocess.Popen(evigene_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        while p.poll() is None and (line := p.stdout.readline()) != "":
            print(line.strip())
        p.wait()

        if p.poll() != 0:
            print(p.stderr.readlines())
            raise Exception("Evigene (tr2aacds.pl) did not complete successfully")
        
    # Append evigene assembly classifier output to metadata
    def run_metadata_appender(self, sqlite_db: Path) -> None:
        transcript_paths = self.set_transcript_paths(self.outdir, self.assembly_fasta)
        transcript_classes = self.extract_transcript_classes(transcript_paths)

        with sqlite3.connect(sqlite_db) as connection:
            self.update_transcripts_table_with_evigene_info(connection, transcript_classes)

    @classmethod
    def set_transcript_paths(cls, outdir: Path, assembly_fasta_path: Path) -> list[Path]:
        transcript_paths = [cls.set_okay_transcript_path(outdir, assembly_fasta_path),
                            cls.set_okalt_transcript_path(outdir, assembly_fasta_path),
                            cls.set_drop_transcript_path(outdir, assembly_fasta_path)]
        return transcript_paths

    @staticmethod
    def set_okay_transcript_path(outdir: Path, assembly_fasta_path: Path) -> Path:
        return outdir / "okayset" / f"{assembly_fasta_path.stem}.okay.tr"

    @staticmethod
    def set_okalt_transcript_path(outdir: Path, assembly_fasta_path: Path) -> Path:
        return outdir / "okayset" / f"{assembly_fasta_path.stem}.okalt.tr"

    @staticmethod
    def set_drop_transcript_path(outdir: Path, assembly_fasta_path: Path) -> Path:
        return outdir / "dropset" / f"{assembly_fasta_path.stem}.drop.tr"

    @staticmethod
    def extract_transcript_classes(transcript_paths: list[Path]) -> defaultdict[dict[Any]]:
        transcript_classes = dict()

        for transcript_class_path in transcript_paths:
            with transcript_class_path.open() as inhandle:
                for line in inhandle:
                    line = line.strip()
                    if line[0] != ">":
                        continue
                    sequence_id = line.split(" ")[0][1:].split("_prefix_")[-1]
                    try:
                        class_drop_info = line.split(" ")[1][:-1]
                        transcript_class = class_drop_info.split(",")[0].replace("evgclass=", "")
                        okay_drop_flag = class_drop_info.split(",")[1]
                    except IndexError:
                        # This IndexError appears to be entirely driven by transcript headers in the
                        # drop file to which class info is not added. However, not all the headers
                        # are wrong, so the file is still processed
                        transcript_class = None
                        okay_drop_flag = None

                    if okay_drop_flag is None:
                        evigene_pass = None
                    elif okay_drop_flag == "okay":
                        evigene_pass = 1
                    else:
                        evigene_pass = 0

                    transcript_classes[sequence_id] = {"transcript_class": transcript_class,
                                                    "evigene_pass": evigene_pass
                                                    }
        return transcript_classes

    @staticmethod
    def update_transcripts_table_with_evigene_info(connection: sqlite3.Connection,
                                                   transcript_classes: defaultdict[dict[Any]]) -> None:
        values = [(data["transcript_class"], data["evigene_pass"], uid) for uid, data in transcript_classes.items()]
        cursor = connection.cursor()
        sql_statement = "UPDATE transcripts SET transcript_class = ?, evigene_pass = ? WHERE uid = ?"
        cursor.executemany(sql_statement, values)
        connection.commit()

class CdsMetadataManager:
    def __init__(self, assembly_fasta: Path, outdir: Path, sqlite_db: Path) -> None:
        self.assembly_fasta = assembly_fasta
        self.outdir = outdir
        self.sqlite_db = sqlite_db

    def run(self) -> None:
        cds_paths = self.set_cds_paths(self.outdir, self.assembly_fasta)
        cds_metadata = self.extract_cds_metadata(cds_paths)

        with sqlite3.connect(self.sqlite_db, timeout=600) as connection:
            self.insert_cds_info_to_cds_table(connection, cds_metadata)

    @classmethod
    def set_cds_paths(cls, outdir: Path, assembly_fasta_path: Path) -> list[Path]:
        transcript_paths = [cls.set_okay_cds_path(outdir, assembly_fasta_path),
                            cls.set_okalt_cds_path(outdir, assembly_fasta_path)]
        return transcript_paths

    @staticmethod
    def set_okay_cds_path(outdir: Path, assembly_fasta_path: Path) -> Path:
        return outdir / "okayset" / f"{assembly_fasta_path.stem}.okay.cds"

    @staticmethod
    def set_okalt_cds_path(outdir: Path, assembly_fasta_path: Path) -> Path:
        return outdir / "okayset" / f"{assembly_fasta_path.stem}.okalt.cds"

    @classmethod
    def extract_cds_metadata(cls, cds_paths: list[Path]) -> list[dict[Any]]:
        cds_metadata = list()

        for cds_path in cds_paths:
            with cds_path.open() as inhandle:
                for line in inhandle:
                    line = line.strip()
                    try:
                        if line[0] != ">":
                            continue
                    except IndexError:
                        continue

                    transcript_id = cls.extract_transcript_id(line)
                    evigene_class = cls.extract_evigene_class(line)
                    strand = cls.extract_strand(line)
                    start, end = cls.extract_cds_coordinates(line)
                    cds_len = end - start + 1

                    cds_metadata.append(
                        {
                        "transcript_id": transcript_id,
                        "evigene_class": evigene_class,
                        "strand": strand,
                        "start": start,
                        "end": end,
                        "cds_len": cds_len
                        }
                    )
        return cds_metadata

    @staticmethod
    def extract_transcript_id(line: str) -> int:
        return int(line.split(" ")[0][1:].split("_prefix_")[-1].replace("utrorf", ""))

    @staticmethod
    def extract_evigene_class(line: str) -> str:
        return line.split(" ")[-1].replace("evgclass=", "").split(",")[0]
    
    @staticmethod
    def extract_strand(line: str) -> str:
        return line.split(" ")[4].split("=")[1][:-1]
    
    @staticmethod
    def extract_cds_coordinates(line: str) -> tuple[int, int]:
        coordinates = line.split(" ")[5].split("=")[1][:-1].split("-")
        first_coordinate = int(coordinates[0])
        second_coordinate = int(coordinates[1])

        if first_coordinate < second_coordinate:
            start = first_coordinate
            end = second_coordinate
        else:
            start = second_coordinate
            end = first_coordinate

        return start, end

    @staticmethod
    def insert_cds_info_to_cds_table(connection: sqlite3.Connection, cds_metadata: list[dict[Any]]) -> None:
        values = [tuple(data.values()) for data in cds_metadata]
        cursor = connection.cursor()
        sql_statement = ("INSERT INTO cds "
                         "(transcript_uid, evigene_class, strand, start, end, length) "
                         "VALUES (?,?,?,?,?,?)")
        cursor.executemany(sql_statement, values)
        connection.commit()

    def run_update_transcript_cds_ids(self) -> None:
        with sqlite3.connect(self.sqlite_db, timeout=600) as connection:
            transcript_cds_id_mapping = self.extract_transcript_cds_id_mapping(connection)
            self.update_transcripts_table_with_cds_ids(connection, transcript_cds_id_mapping)

    @staticmethod
    def extract_transcript_cds_id_mapping(connection: sqlite3.Connection) -> dict[str]:
        cursor = connection.cursor()
        sql_query = f"SELECT uid,transcript_uid FROM cds"
        cursor.execute(sql_query)
        cds_transcript_id_mapping = {row[0]: row[1] for row in cursor.fetchall()}

        transcript_cds_id_mapper = defaultdict(list)
        for cds_id, transcript_id in cds_transcript_id_mapping.items():
            cds_id = str(cds_id)
            transcript_cds_id_mapper[transcript_id].append(cds_id)
        return {transcript_id: ";".join(cds_ids) for transcript_id, cds_ids in transcript_cds_id_mapper.items()}

    @staticmethod
    def update_transcripts_table_with_cds_ids(connection: sqlite3.Connection,
                                                   transcript_cds_id_mapping: dict[str]) -> None:
        values = [(cds_ids, uid) for uid, cds_ids in transcript_cds_id_mapping.items()]
        cursor = connection.cursor()
        sql_statement = "UPDATE transcripts SET cds_ids = ? WHERE uid = ?"
        cursor.executemany(sql_statement, values)
        connection.commit()

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-assembly_fasta", type=str, required=True)
    parser.add_argument("-outdir", type=str, required=True)
    parser.add_argument("-sqlite_db", type=str, required=True)
    parser.add_argument("-transcriptome", type=str, required=True)
    parser.add_argument("-prefix_column", type=str, required=True)
    parser.add_argument("-cpus", type=int, required=False, default=1)
    parser.add_argument("-mem", type=int, required=False, default=1_000)
    parser.add_argument("-run_evigene", action="store_true", required=False)
    parser.add_argument("-phetero", type=int, required=False)
    parser.add_argument("-minaa", type=int, required=False)
    parser.add_argument("-run_transcript_metadata_appender", action="store_true", required=False)
    parser.add_argument("-run_cds_and_metadata", action="store_true", required=False)
    parser.add_argument("-update_transcript_cds_ids", action="store_true", required=False)
    args = parser.parse_args()

    outdir = Path(args.outdir) / args.transcriptome

    if args.run_evigene:
        if outdir.is_dir():
            rmtree(outdir)
        outdir.mkdir()

        ifm = InputFastaManager(Path(args.assembly_fasta), outdir,
                                Path(args.sqlite_db), args.transcriptome, args.prefix_column)
        tmp_fasta = ifm.run()
        assembly_fasta = tmp_fasta
    else:
        tmp_fasta = None
        assembly_fasta = Path(args.assembly_fasta)

    em = EvigeneManager(assembly_fasta, outdir, args.cpus, args.mem)

    if args.run_evigene:
        print("\nRunning evigene assembly classifier")
        em.run_assembly_classifier(args.phetero, args.minaa)

    if args.run_transcript_metadata_appender:
        print("\nRunning transcript metadata appender")
        em.run_metadata_appender(Path(args.sqlite_db))

    if tmp_fasta:
        tmp_fasta.unlink()
        tmp_fasta.parent.rmdir()

    if args.run_cds_and_metadata:
        print("\nRunning CDS and metadata extractor")
        cmm = CdsMetadataManager(assembly_fasta, outdir, Path(args.sqlite_db))
        cmm.run()

    if args.update_transcript_cds_ids:
        print("\nUpdating transcripts table with CDS ids")
        cmm = CdsMetadataManager(assembly_fasta, outdir, Path(args.sqlite_db))
        cmm.run_update_transcript_cds_ids()

    print("\nFinished\n")