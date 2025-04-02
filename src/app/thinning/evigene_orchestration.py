from argparse import ArgumentParser
from collections import defaultdict
from contextlib import contextmanager
from csv import DictReader, DictWriter
from os import chdir, environ, getcwd
from pathlib import Path
from shutil import copyfileobj
import subprocess
from tempfile import mkdtemp, NamedTemporaryFile
from typing import Any, Union

@contextmanager
def temporarily_change_working_directory(new_directory: Path):
    starting_directory = getcwd()
    try:
        chdir(new_directory)
        yield
    finally:
        chdir(starting_directory)

class PrefixedFastaManager:
    def __init__(self, assembly_fasta: Path, outdir: Path, prefix_column: str, metadata: Path) -> None:
        self.assembly_fasta = assembly_fasta
        self.outdir = outdir
        self.prefix_column = prefix_column
        self.metadata = metadata

    def run(self) -> Path:
        sequence_prefix_mapping = self.extract_sequence_prefix_mapping(self.prefix_column, self.metadata)
        tmp_prefixed_fasta = self.write_temporary_prefixed_fasta(self.assembly_fasta, self.outdir, sequence_prefix_mapping)
        return tmp_prefixed_fasta

    @staticmethod
    def extract_sequence_prefix_mapping(prefix_column: str, metadata: Path) -> dict[str]:
        print("\nExtracting sequence id to prefix mapping")
        sequence_prefix_mapping = {}
        with metadata.open() as inhandle:
            reader = DictReader(inhandle)
            for data in reader:
                sequence_id = data["sequence_id"]
                sequence_prefix = data[prefix_column]
                sequence_prefix_mapping[sequence_id] = sequence_prefix
        return sequence_prefix_mapping
    
    @staticmethod
    def write_temporary_prefixed_fasta(assembly_fasta: Path, outdir: Path, sequence_prefix_mapping: dict[str]) -> Path:
        print("Writing temporary prefixed fasta")
        temp_dir = mkdtemp(dir=outdir)
        outfile = Path(temp_dir) / assembly_fasta.name
        with assembly_fasta.open() as inhandle, outfile.open("w") as outhandle:
            for line in inhandle:
                line = line.strip()
                if line[0] != ">":
                    outhandle.write(line + "\n")
                    continue
                sequence_id = line[1:]
                prefix = sequence_prefix_mapping[sequence_id]
                new_header = ">" + "_".join([prefix, "prefix", sequence_id]) # TODO: REMOVE SEQ!
                outhandle.write(new_header + "\n")
        return outfile

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
    def run_metadata_appender(self, metadata_path: Path) -> None:
        transcript_paths = self.set_transcript_paths(self.outdir, self.assembly_fasta)
        transcript_classes = self.extract_transcript_classes(transcript_paths)

        tmpfile_path = self.write_appended_metadata_to_tempfile(transcript_classes, metadata_path, self.outdir)
        self.copy_file(tmpfile_path, metadata_path)
        tmpfile_path.unlink()

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
                        transcript_class = ""
                        okay_drop_flag = ""
                    evigene_pass = True if okay_drop_flag == "okay" else False

                    transcript_classes[sequence_id] = {"transcript_class": transcript_class,
                                                    "evigene_pass": evigene_pass
                                                    }
        return transcript_classes

    @staticmethod
    def write_appended_metadata_to_tempfile(transcript_classes: defaultdict[dict[Any]], metadata_path: Path, outdir: Path) -> Path:
        print("Starting to write to tmpfile\n(This may take a moment)")
        new_fields = list(transcript_classes[next(iter(transcript_classes))].keys())

        with metadata_path.open() as inhandle, NamedTemporaryFile(dir=outdir, mode="w", delete=False) as tmpfile:
            reader = DictReader(inhandle)
            write_field_names = reader.fieldnames + [field for field in new_fields if field not in reader.fieldnames]

            writer = DictWriter(tmpfile, fieldnames=write_field_names)
            writer.writeheader()
            for data in reader:
                sequence_id = data["sequence_id"]
                try:
                    transcript_class = transcript_classes[sequence_id]
                    data.update(transcript_class)
                except KeyError:
                    print(f"Key missing: {sequence_id}")
                writer.writerow(data)
        return Path(tmpfile.name)

    @staticmethod
    def copy_file(infile_path: Path, outfile_path: Path) -> None:
        print("Starting to copy tmpfile to metadata file")
        with infile_path.open("rb") as inhandle, outfile_path.open("wb") as outhandle:
            copyfileobj(inhandle, outhandle)

class CdsMetadataManager:
    def __init__(self, assembly_fasta: Path, outdir: Path, transcript_metadata: Path, cds_metadata: Path) -> None:
        self.assembly_fasta = assembly_fasta
        self.outdir = outdir
        self.transcript_metadata_path = transcript_metadata
        self.cds_metadata_path = cds_metadata

    def run(self) -> None:
        cds_paths = self.set_cds_paths(self.outdir, self.assembly_fasta)
        cds_metadata = self.extract_cds_metadata(cds_paths)
        self.write_cds_metadata(self.cds_metadata_path, cds_metadata)

        transcript_cds_id_mapping = self.extract_transcript_cds_id_mapping(cds_metadata)

        tmpfile_path = self.write_appended_metadata_to_tempfile(transcript_cds_id_mapping, self.transcript_metadata_path, self.outdir)
        self.copy_file(tmpfile_path, self.transcript_metadata_path)
        tmpfile_path.unlink()

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

        cds_id = 0
        for cds_path in cds_paths:
            with cds_path.open() as inhandle:
                for line in inhandle:
                    line = line.strip()
                    try:
                        if line[0] != ">":
                            continue
                    except IndexError:
                        continue

                    cds_id += 1
                    transcript_id = cls.extract_transcript_id(line)
                    evigene_class = cls.extract_evigene_class(line)
                    strand = cls.extract_strand(line)
                    start, end = cls.extract_cds_coordinates(line)
                    cds_len = end - start + 1

                    cds_metadata.append(
                        {
                        "cds_id": str(cds_id),
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
    def extract_transcript_id(line: str) -> str:
        return line.split(" ")[0][1:].split("_prefix_")[-1].replace("utrorf", "")

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
    def set_expected_evigene_classes() -> set[str]:
        expected_classes = {"althi", "althi1", "althia2", "althinc", "altmfrag", "altmfraga2", "altmid", "altmida2",
                            "main", "maina2", "mainnc",
                            "noclass", "noclassa2", "noclassa2nc", "noclassnc",
                            "parthi", "parthi1", "parthia2",
                            "perfdupl", "perffrag",
                            "smallorf"}
        additional_classes = {"altmidfrag", "altmidfraga2"}
        expected_classes = expected_classes | additional_classes

        return expected_classes

    @staticmethod
    def write_cds_metadata(cds_metadata_path: Path, cds_metadata: list[dict[Any]]) -> None:
        print("\nWriting CDS metadata")
        with cds_metadata_path.open("w") as outhandle:
            writer = DictWriter(outhandle, fieldnames=cds_metadata[0].keys())
            writer.writeheader()
            for data in cds_metadata:
                writer.writerow(data)

    @staticmethod
    def extract_transcript_cds_id_mapping(cds_metadata: list[dict[Any]]) -> dict[str]:
        transcript_cds_id_mapper = defaultdict(list)
        for data in cds_metadata:
            cds_id = data["cds_id"]
            transcript_id = data["transcript_id"]
            transcript_cds_id_mapper[transcript_id].append(cds_id)
        return {transcript_id: ";".join(cds_ids) for transcript_id, cds_ids in transcript_cds_id_mapper.items()}

    @staticmethod
    def write_appended_metadata_to_tempfile(transcript_cds_id_mapping: dict[str], metadata_path: Path, outdir: Path) -> Path:
        print("Starting to write appended metadata to tmpfile\n")
        new_fields = ["cds_ids"]

        with metadata_path.open() as inhandle, NamedTemporaryFile(dir=outdir, mode="w", delete=False) as tmpfile:
            reader = DictReader(inhandle)
            write_field_names = reader.fieldnames + [field for field in new_fields if field not in reader.fieldnames]

            writer = DictWriter(tmpfile, fieldnames=write_field_names)
            writer.writeheader()
            for data in reader:
                transcript_id = data["sequence_id"]
                try:
                    cds_ids = transcript_cds_id_mapping[transcript_id]
                    data.update({"cds_ids": cds_ids})
                except KeyError:
                    pass
                writer.writerow(data)
        return Path(tmpfile.name)

    @staticmethod
    def copy_file(infile_path: Path, outfile_path: Path) -> None:
        print("Starting to copy tmpfile to metadata file")
        with infile_path.open("rb") as inhandle, outfile_path.open("wb") as outhandle:
            copyfileobj(inhandle, outhandle)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-assembly_fasta", type=str, required=True)
    parser.add_argument("-outdir", type=str, required=True)
    parser.add_argument("-cpus", type=int, required=False, default=1)
    parser.add_argument("-mem", type=int, required=False, default=1_000)
    parser.add_argument("-run_evigene", action="store_true", required=False)
    parser.add_argument("-phetero", type=int, required=False)
    parser.add_argument("-minaa", type=int, required=False)
    parser.add_argument("-prefix_column", type=str, required=False)
    parser.add_argument("-metadata", type=str, required=False)
    parser.add_argument("-run_metadata_appender", action="store_true", required=False)
    parser.add_argument("-cds_metadata", type=str, required=False)
    args = parser.parse_args()

    if args.prefix_column and args.metadata:
        pfm = PrefixedFastaManager(Path(args.assembly_fasta), Path(args.outdir), args.prefix_column, Path(args.metadata))
        tmp_prefixed_fasta = pfm.run()
    elif args.prefix_column and not args.metadata:
        raise Exception("Prefix detected, but not metadata.\nCancelling run")
    else:
       tmp_prefixed_fasta = None

    if tmp_prefixed_fasta:
        assembly_fasta = tmp_prefixed_fasta
    else:
        assembly_fasta = Path(args.assembly_fasta)

    em = EvigeneManager(assembly_fasta, Path(args.outdir), args.cpus, args.mem)

    if args.run_evigene:
        print("\nRunning evigene assembly classifier")
        em.run_assembly_classifier(args.phetero, args.minaa)
    if args.run_metadata_appender:
        print("\nRunning metadata appender")
        em.run_metadata_appender(Path(args.metadata))

    if tmp_prefixed_fasta:
        tmp_prefixed_fasta.unlink()
        tmp_prefixed_fasta.parent.rmdir()

    if args.cds_metadata:
        print("\nRunning CDS metadata extractor")
        cmm = CdsMetadataManager(assembly_fasta, Path(args.outdir), Path(args.metadata), Path(args.cds_metadata))
        cmm.run()

    print("\nFinished\n")