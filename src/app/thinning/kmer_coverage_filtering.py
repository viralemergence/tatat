from argparse import ArgumentParser
from csv import reader
import matplotlib.pyplot as plt # type: ignore
import numpy as np # type: ignore
import pandas as pd # type: ignore
from pathlib import Path
import seaborn as sns # type: ignore
from typing import Any, Iterator, Union

class KmerCoverageManager:
    def __init__(self, assembly_path: Path, coverage_cutoff: float) -> None:
        self.assembly_path = assembly_path
        self.coverage_cutoff = coverage_cutoff
    
    def run_filter(self, outdir: Path) -> None:
        print("Starting kmer coverage filtering")
        print(f"kmer coverage cutoff: {self.coverage_cutoff}")
        outpath = self.set_outpath(self.assembly_path, self.coverage_cutoff, outdir)

        removed = 0
        with outpath.open("w") as outhandle:
            for i, fasta_seq in enumerate(self.fasta_chunker(self.assembly_path), 1):
                fasta_header = fasta_seq[0].replace(">", "")
                assembly_info = self.extract_assembly_info(fasta_header)

                if assembly_info["cov"] < self.coverage_cutoff:
                    removed += 1
                    continue
                for line in fasta_seq:
                    outhandle.write(line + "\n")

        percent_removed = round((removed/i)*100, 2)
        retained = i - removed
        print(f"Sequences removed: {removed}")
        print(f"Percent of sequences removed: {percent_removed}%")
        print(f"Sequences retained: {retained}\n")

    @staticmethod
    def set_outpath(assembly_path: Path, coverage_cutoff: float, outdir: Path) -> Path:
        return outdir / f"{assembly_path.stem}_{coverage_cutoff}_filtered.fasta"

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

    @staticmethod
    def extract_assembly_info(fasta_header: str) -> dict[Any]:
        assembly_info = fasta_header.split("_")
        assembly_dict = dict(zip(assembly_info[::2], assembly_info[1::2]))
        assembly_dict["length"] = int(assembly_dict["length"])
        assembly_dict["cov"] = round(float(assembly_dict["cov"]), 2)
        return assembly_dict

    def run_qc(self, qc_dir: Path, sampling_size: Union[int, None]) -> None:
        print("Starting kmer coverage QC")
        qc_graph_path = self.set_qc_graph_path(qc_dir, self.assembly_path, self.coverage_cutoff)

        coverage_and_length_data = self.extract_coverages_and_lengths(self.assembly_path)
        coverage_and_length_data = self.transform_coverages_and_lengths(coverage_and_length_data, sampling_size)

        print("Starting QC graph")
        self.generate_jointplot(coverage_and_length_data, self.coverage_cutoff, qc_graph_path)

    @staticmethod
    def set_qc_graph_path(qc_dir: Path, assembly_path: Path, coverage_cutoff: float) -> Path:
        return qc_dir / f"{assembly_path.stem}_{coverage_cutoff}_kmer_qc.jpeg"
    
    @classmethod
    def extract_coverages_and_lengths(cls, assembly_path: Path) -> pd.DataFrame:
        coverages = []
        lengths = []
        for fasta_seq in cls.fasta_chunker(assembly_path):
            fasta_header = fasta_seq[0].replace(">", "")
            assembly_info = cls.extract_assembly_info(fasta_header)
            coverages.append(assembly_info["cov"])
            lengths.append(assembly_info["length"])
        return pd.DataFrame({"kmer-Coverage": coverages, "Length": lengths})
    
    @staticmethod
    def transform_coverages_and_lengths(coverage_and_length_data: pd.DataFrame, sampling_size: int) -> pd.DataFrame:
        if not sampling_size is None:
            coverage_and_length_data = coverage_and_length_data.sample(n=sampling_size)
        coverage_and_length_data = np.log10(coverage_and_length_data)
        return coverage_and_length_data
    
    @staticmethod
    def generate_jointplot(coverage_and_length_data: pd.DataFrame, coverage_cutoff: float, qc_graph_path: Path) -> None:
        histogram = sns.jointplot(coverage_and_length_data, x="kmer-Coverage", y="Length", kind="hist",
                                  marginal_kws={"bins": 100}, bins=50)
        plt.axvline(x=np.log10(coverage_cutoff), color="black")
        fig = histogram.figure
        fig.savefig(qc_graph_path, bbox_inches="tight", dpi=300)
        plt.close(fig)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-assembly_path", type=str, required=True)
    parser.add_argument("-coverage_cutoff", type=float, required=True)
    parser.add_argument("-outdir", type=str, required=False)
    parser.add_argument("-qc_dir", type=str, required=False)
    parser.add_argument("-sampling_size", type=int, required=False)
    args = parser.parse_args()

    print()
    kcm = KmerCoverageManager(Path(args.assembly_path), args.coverage_cutoff)
    if not args.outdir is None:
        kcm.run_filter(Path(args.outdir))
    if not args.qc_dir is None:
        kcm.run_qc(Path(args.qc_dir), args.sampling_size)