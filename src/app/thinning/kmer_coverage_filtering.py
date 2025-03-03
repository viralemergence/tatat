from argparse import ArgumentParser
from csv import reader
from pathlib import Path
from typing import Any, Iterator

class KmerCoverageManager:
    def __init__(self, assembly_path: Path, coverage_cutoff: float, outdir: Path) -> None:
        self.assembly_path = assembly_path
        self.coverage_cutoff = coverage_cutoff
        self.outpath = self.set_outpath(self.assembly_path, self.coverage_cutoff, outdir)

    @staticmethod
    def set_outpath(assembly_path: Path, coverage_cutoff: float, outdir: Path) -> Path:
        return outdir / f"{assembly_path.stem}_{coverage_cutoff}_filtered.fasta"
    
    def run_filter(self):
        removed = 0
        with self.outpath.open("w") as outhandle:
            for fasta_seq in self.fasta_chunker(self.assembly_path):
                fasta_header = fasta_seq[0].replace(">", "")
                assembly_info = self.extract_assembly_info(fasta_header)

                if assembly_info["cov"] < self.coverage_cutoff:
                    removed += 1
                    continue
                for line in fasta_seq:
                    outhandle.write(line + "\n")
        print(f"Sequences removed: {removed}")

    def run_qc(self):
        coverages = []
        lengths = []
        with self.outpath.open("w") as outhandle:
            for fasta_seq in self.fasta_chunker(self.assembly_path):
                fasta_header = fasta_seq[0].replace(">", "")
                assembly_info = self.extract_assembly_info(fasta_header)

                coverages.append(assembly_info["cov"])
                lengths.append(assembly_info["length"])

    @staticmethod
    def extract_assembly_info(fasta_header: str) -> dict[Any]:
        assembly_info = fasta_header.split("_")
        assembly_dict = dict(zip(assembly_info[::2], assembly_info[1::2]))
        assembly_dict["length"] = int(assembly_dict["length"])
        assembly_dict["cov"] = round(float(assembly_dict["cov"]), 2)
        return assembly_dict

    @staticmethod
    def fasta_chunker(fasta_path: Path) -> Iterator[list[str]]:
        fasta_seq = []
        first_chunk = True
        with fasta_path.open() as inhandle:
            reader_iterator = reader(inhandle)
            for line in reader_iterator:
                line = line[0]
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

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-assembly_path", type=str, required=True)
    parser.add_argument("-coverage_cutoff", type=float, required=True)
    parser.add_argument("-outdir", type=str, required=True)
    args = parser.parse_args()

    kcm = KmerCoverageManager(Path(args.assembly_path), args.coverage_cutoff, Path(args.outdir))
    kcm.run_filter()