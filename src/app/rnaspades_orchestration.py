import argparse
from pathlib import Path
import subprocess

class FastqAssemblyPathManager:
    def __init__(self, fastq_dir: Path, unique_identifier: str, assembly_dir: Path, collated_dir: Path) -> None:
        self.ui_fastq_paths = self.extract_ui_fastq_paths(fastq_dir, unique_identifier)
        self.assembly_dir = self.generate_assembly_dir(assembly_dir, unique_identifier)
        self.output_collated_path = collated_dir / f"{unique_identifier}.fasta"

    @classmethod
    def extract_ui_fastq_paths(cls, fastq_dir: Path, unique_identifier: str) -> list[Path]:
        fastq_paths = cls.get_file_list(fastq_dir)
        return [path for path in fastq_paths if unique_identifier in path.stem]

    @staticmethod
    def get_file_list(directory: Path) -> list[Path]:
        return sorted([file for file in directory.iterdir()])

    @staticmethod
    def generate_assembly_dir(assembly_dir: Path, unique_identifier: str) -> Path:
        if not assembly_dir.exists():
            raise Exception(f"assembly_dir arg directory does not exist: {assembly_dir}")
        assembly_dir = assembly_dir / f"{unique_identifier}"
        if not assembly_dir.exists():
            print(f"Unique identifier assembly outdir does not exist.\nGenerating now: {assembly_dir}")
            assembly_dir.mkdir()
        return assembly_dir

class RnaspadesManager:
    def __init__(self, fastq_paths: list[Path], assembly_dir: Path, output_collated_path: Path, cpus: int, memory: int) -> None:
        self.fastq_paths = fastq_paths
        self.assembly_dir = assembly_dir
        self.output_collated_path = output_collated_path
        self.cpus = cpus
        self.memory = memory

    def run(self) -> None:
        rnaspades_assembly_path = self.run_rnaspades()
        self.move_rnaspades_assembly(rnaspades_assembly_path, self.output_collated_path)

    def run_rnaspades(self) -> Path:
        print("Starting rnaSPAdes")
        if len(self.fastq_paths) == 2:
            rnaspades_command = ["rnaspades.py",
                               "-1", self.fastq_paths[0],
                               "-2", self.fastq_paths[1],
                               "-o", self.assembly_dir,
                               "-t", f"{self.cpus}",
                               "-m", f"{self.memory}"]
        if len(self.fastq_paths) == 1:
            rnaspades_command = ["rnaspades.py",
                               "-s", self.fastq_paths[0],
                               "-o", self.assembly_dir,
                               "-t", f"{self.cpus}",
                               "-m", f"{self.memory}"]

        p = subprocess.Popen(rnaspades_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        while p.poll() is None and (line := p.stdout.readline()) != "":
            print(line.strip())
        p.wait()
        print(f"Exit code: {p.poll()}")
        for line in p.stderr:
            print(line)

        if p.poll() != 0:
            raise Exception("rnaSPAdes did not complete successfully")

        rnaspades_fasta_path = Path(f"{self.assembly_dir}/transcripts.fasta")
        return rnaspades_fasta_path

    @staticmethod
    def move_rnaspades_assembly(rnaspades_assembly_path: Path, output_collated_path: Path) -> None:
        print("Moving rnaSPAdes assembly")

        mv_command = ["mv", rnaspades_assembly_path, output_collated_path]
        result = subprocess.run(mv_command, capture_output=True, text=True)
        print(result.stdout)
        print(result.stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-fastq_dir", type=str, required=True)
    parser.add_argument("-assembly_dir", type=str, required=True)
    parser.add_argument("-collated_dir", type=str, required=True)
    parser.add_argument("-unique_identifier", type=str, required=True)
    parser.add_argument("-cpus", type=int, required=True)
    parser.add_argument("-memory", type=int, required=True)
    args = parser.parse_args()

    fapm = FastqAssemblyPathManager(Path(args.fastq_dir), args.unique_identifier, Path(args.assembly_dir), Path(args.collated_dir))

    rm = RnaspadesManager(fapm.ui_fastq_paths, fapm.assembly_dir, fapm.output_collated_path, args.cpus, args.memory)
    rm.run()