from argparse import ArgumentParser
from pathlib import Path
import subprocess

class FastqPathManager:
    def __init__(self, fastq_dir: Path, unique_identifier: str, outdir: Path) -> None:
        self.ui_fastq_paths = self.extract_ui_fastq_paths(fastq_dir, unique_identifier)
        self.output_fastq_paths = self.set_output_paths(self.ui_fastq_paths, outdir)

    @classmethod
    def extract_ui_fastq_paths(cls, fastq_dir: Path, unique_identifier: str) -> list[Path]:
        fastq_paths = cls.get_file_list(fastq_dir)
        return [path for path in fastq_paths if unique_identifier in path.stem]

    @staticmethod
    def get_file_list(directory: Path) -> list[Path]:
        return sorted([file for file in directory.iterdir()])

    @staticmethod
    def set_output_paths(input_fastqs: list[Path], outdir: Path) -> list[Path]:
        output_fastqs = []
        for fastq in input_fastqs:
            outpath_stem = fastq.stem.replace(".fastq", "")
            outpath = outdir / f"{outpath_stem}_fastp.fastq.gz"
            output_fastqs.append(outpath)
        return output_fastqs

class FastpManager:
    def __init__(self, input_fastqs: list[Path], output_fastqs: list[Path]) -> None:
        self.input_fastqs = input_fastqs
        self.output_fastqs = output_fastqs

    def run_fastp(self, r1_adapter: str, r2_adapter: str) -> None:
        if len(self.input_fastqs) == 2:
            fastp_command = ["fastp", "-V",
                            "-i", self.input_fastqs[0], "-I", self.input_fastqs[1],
                            "-o", self.output_fastqs[0], "-O", self.output_fastqs[1],
                            "--adapter_sequence", r1_adapter,
                            "--adapter_sequence_r2", r2_adapter]
        if len(self.input_fastqs) == 1:
            fastp_command = ["fastp", "-V",
                            "-i", self.input_fastqs[0],
                            "-o", self.output_fastqs[0],
                            "--adapter_sequence", r1_adapter]

        p = subprocess.Popen(fastp_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        while p.poll() is None and (line := p.stderr.readline()) != "":
            print(line.strip())
        p.wait()
        print(f"Exit code: {p.poll()}")

        if p.poll() != 0:
            raise Exception("Fastp did not complete successfully")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--fastq_dir", type=str, required=True)
    parser.add_argument("-o", "--outdir", type=str, required=True)
    parser.add_argument("-u", "--unique_identifier", type=str, required=True)
    parser.add_argument("-r1a", "--r1_adapter", type=str, required=True)
    parser.add_argument("-r2a", "--r2_adapter", type=str, required=True)
    args = parser.parse_args()

    fpm = FastqPathManager(Path(args.fastq_dir), args.unique_identifier, Path(args.outdir))

    fm = FastpManager(fpm.ui_fastq_paths, fpm.output_fastq_paths)
    fm.run_fastp(args.r1_adapter, args.r2_adapter)