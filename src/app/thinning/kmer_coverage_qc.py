from argparse import ArgumentParser
from csv import DictReader, DictWriter
import matplotlib.pyplot as plt # type: ignore
import numpy as np # type: ignore
import pandas as pd # type: ignore
from pathlib import Path
import seaborn as sns # type: ignore
from shutil import copyfileobj
from tempfile import NamedTemporaryFile
from typing import Union

class KmerCoverageManager:
    def __init__(self, metadata_path: Path, coverage_cutoff: float) -> None:
        self.metadata_path = metadata_path
        self.coverage_cutoff = coverage_cutoff
    
    def run_update_kmer_coverage_pass(self) -> None:
        print("\nStarting 'kmer_coverage_pass' flagging\n(This may take a moment)")
        print(f"kmer coverage cutoff: {self.coverage_cutoff}")

        tmpfile_path = self.write_appended_kmer_coverage_pass_metadata_to_tempfile(self.metadata_path, self.coverage_cutoff)
        self.copy_file(tmpfile_path, self.metadata_path)
        tmpfile_path.unlink()

    @staticmethod
    def write_appended_kmer_coverage_pass_metadata_to_tempfile(metadata_path: Path, coverage_cutoff: float) -> Path:
        new_fields = ["kmer_coverage_pass"]
        with metadata_path.open() as inhandle, NamedTemporaryFile(dir=metadata_path.parent, mode="w", delete=False) as tmpfile:
            reader = DictReader(inhandle)
            write_field_names = reader.fieldnames + [field for field in new_fields if field not in reader.fieldnames]

            writer = DictWriter(tmpfile, fieldnames=write_field_names)
            writer.writeheader()

            for data in reader:
                data["kmer_coverage_pass"] = True if float(data["kmer_coverage"]) > coverage_cutoff else False
                writer.writerow(data)

        return Path(tmpfile.name)

    @classmethod
    def copy_file(cls, infile_path: Path, outfile_path: Path) -> None:
        print("Starting to copy tmpfile to metadata file")
        with infile_path.open("rb") as inhandle, outfile_path.open("wb") as outhandle:
            copyfileobj(inhandle, outhandle)

    def run_qc(self, qc_dir: Path, sampling_size: Union[int, None]) -> None:
        qc_graph_path = self.set_qc_graph_path(qc_dir, self.metadata_path, self.coverage_cutoff)

        print("\nExtracting kmer coverage and assembly length for QC graph")
        coverage_and_length_data = self.extract_coverages_and_lengths(self.metadata_path)
        self.report_count_below_kmer_cutoff(coverage_and_length_data, self.coverage_cutoff)
        coverage_and_length_data = self.transform_coverages_and_lengths(coverage_and_length_data, sampling_size)

        print("\nGenerating QC graph")
        self.generate_jointplot(coverage_and_length_data, self.coverage_cutoff, qc_graph_path)

    @staticmethod
    def set_qc_graph_path(qc_dir: Path, metadata_path: Path, coverage_cutoff: float) -> Path:
        return qc_dir / f"{metadata_path.stem}_{coverage_cutoff}_kmer_qc.jpeg"
    
    @classmethod
    def extract_coverages_and_lengths(cls, metadata_path: Path) -> pd.DataFrame:
        coverages = []
        lengths = []
        with metadata_path.open() as inhandle:
            reader = DictReader(inhandle)
            for data in reader:
                coverages.append(float(data["kmer_coverage"]))
                lengths.append(int(data["length"]))
        return pd.DataFrame({"kmer-Coverage": coverages, "Length": lengths})
    
    @staticmethod
    def report_count_below_kmer_cutoff(coverage_and_length_data: pd.DataFrame, coverage_cutoff: float) -> None:
        count_below_cutoff = (coverage_and_length_data["kmer-Coverage"] < coverage_cutoff).sum()
        total = len(coverage_and_length_data)
        count_below_cutoff_percent = round((count_below_cutoff/total)*100, 2)
        print(f"Number of assemblies below cutoff: {count_below_cutoff}")
        print(f"Total number of assemblies: {total}")
        print(f"Percent of assemblies below cutoff: {count_below_cutoff_percent}%")

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
    parser.add_argument("-metadata_path", type=str, required=True)
    parser.add_argument("-coverage_cutoff", type=float, required=True)
    parser.add_argument("-update_kmer_coverage_pass", action="store_true", required=False)
    parser.add_argument("-qc_dir", type=str, required=False)
    parser.add_argument("-sampling_size", type=int, required=False)
    args = parser.parse_args()

    kcm = KmerCoverageManager(Path(args.metadata_path), args.coverage_cutoff)
    if args.update_kmer_coverage_pass:
        kcm.run_update_kmer_coverage_pass()
    if args.qc_dir:
        kcm.run_qc(Path(args.qc_dir), args.sampling_size)
    print("Finished")