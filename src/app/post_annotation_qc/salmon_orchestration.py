import argparse
from pathlib import Path
import subprocess

class FastqPathManager:
    def __init__(self, fastq_dir: Path, sra: str) -> None:
        self.fastq_dir = fastq_dir
        self.sra = sra
        self.sra_fastq_paths = self.extract_sra_fastq_paths(self.fastq_dir, self.sra)

    @classmethod
    def extract_sra_fastq_paths(cls, fastq_dir: Path, sra: str) -> list[Path]:
        fastq_paths = cls.get_file_list(fastq_dir)
        return [path for path in fastq_paths if sra in path.stem]

    @staticmethod
    def get_file_list(directory: Path) -> list[Path]:
        return sorted([file for file in directory.iterdir()])

class SalmonManager:
    def __init__(self, fastq_paths: list[Path], outdir: Path, collated_dir: Path,
                 sra_number: str, salmon_index: Path, cpus: int) -> None:
        self.fastq_paths = fastq_paths
        self.outdir = self.generate_sra_outdir(outdir, sra_number)
        self.collated_dir = collated_dir
        self.sra_number = sra_number
        self.salmon_index = salmon_index
        self.cpus = cpus

    @staticmethod
    def generate_sra_outdir(outdir: Path, sra_number: str) -> Path:
        if not outdir.exists():
            raise Exception(f"outdir arg directory does not exist: {outdir}")
        sra_count_outdir = outdir / f"{sra_number}"
        if not sra_count_outdir.exists():
            print(f"SRA count outdir does not exist.\nGenerating now: {sra_count_outdir}")
            sra_count_outdir.mkdir()
        return sra_count_outdir

    def run(self) -> None:
        salmon_counts_path = self.run_salmon()
        self.move_salmon_counts(salmon_counts_path, self.collated_dir, self.sra_number)

    def run_salmon(self) -> Path:
        print("Starting Salmon")
        
        if len(self.fastq_paths) == 2:
            salmon_command = ["salmon", "quant",
                              "-i", self.salmon_index,
                              "-l", "A", "-p", f"{self.cpus}", "--validateMappings",
                              "-1", self.fastq_paths[0],
                              "-2", self.fastq_paths[1],
                              "-o", self.outdir
                              ]
        if len(self.fastq_paths) == 1:
            salmon_command = ["salmon", "quant",
                              "-i", self.salmon_index,
                              "-l", "A", "-p", f"{self.cpus}", "--validateMappings",
                              "-r", self.fastq_paths[0],
                              "-o", self.outdir
                              ]

        p = subprocess.Popen(salmon_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        while p.poll() is None and (line := p.stderr.readline()) != "":
            print(line.strip())
        p.wait()
        print(f"Exit code: {p.poll()}")
        for line in p.stderr:
            print(line.strip())

        if p.poll() != 0:
            raise Exception("Salmon did not complete successfully")

        salmon_counts_path = Path(f"{self.outdir}/quant.sf")
        return salmon_counts_path

    @staticmethod
    def move_salmon_counts(salmon_counts_path: Path, collated_dir: Path, sra_number: str) -> None:
        print("Moving salmon counts")
        new_salmon_counts_name = f"{sra_number}_salmon.txt"
        destination_path = f"{collated_dir}/{new_salmon_counts_name}"

        mv_command = ["mv", salmon_counts_path, destination_path]
        result = subprocess.run(mv_command, capture_output=True, text=True)
        print(result.stdout)
        print(result.stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-fastq_dir", type=str, required=True)
    parser.add_argument("-sra", type=str, required=True)
    parser.add_argument("-outdir", type=str, required=True)
    parser.add_argument("-collated_dir", type=str, required=True)
    parser.add_argument("-salmon_index", type=str, required=True)
    parser.add_argument("-cpus", type=int, required=True)
    args = parser.parse_args()

    fpm = FastqPathManager(Path(args.fastq_dir), args.sra)

    sm = SalmonManager(fpm.sra_fastq_paths, Path(args.outdir), Path(args.collated_dir),
                       args.sra, Path(args.salmon_index), args.cpus)
    sm.run()