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
        transcript_classes = defaultdict(dict)

        for transcript_class_path in transcript_paths:
            with transcript_class_path.open() as inhandle:
                while (line := inhandle.readline().strip()):
                    if line[0] != ">":
                        continue
                    sequence_id = line.split(" ")[0][1:].split("_prefix_")[-1]
                    try:
                        class_drop_info = line.split(" ")[1][:-1]
                        transcript_class = class_drop_info.split(",")[0].replace("evgclass=", "")
                        okay_drop_flag = class_drop_info.split(",")[1]
                    except IndexError:
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
                except KeyError:
                    print(f"Key missing: {sequence_id}")
                    continue
                data.update(transcript_class)
                writer.writerow(data)
        return Path(tmpfile.name)

    @classmethod
    def copy_file(cls, infile_path: Path, outfile_path: Path) -> None:
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

    print("\nFinished\n")