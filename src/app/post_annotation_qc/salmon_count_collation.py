import argparse
from csv import reader
import pandas as pd # type: ignore
from pathlib import Path

class SalmonCountCollator:
    def __init__(self, counts_dir: Path, outdir: Path) -> None:
        self.counts_dir = counts_dir
        self.outdir = outdir

    def run(self) -> None:
        count_paths = self.get_file_list(self.counts_dir)
        gene_counts = self.extract_gene_counts(count_paths).round(2)
        outpath = self.outdir / "collated_salmon_counts.csv"
        gene_counts.to_csv(outpath)

    @staticmethod
    def get_file_list(directory: Path) -> list[Path]:
        return sorted([file for file in directory.iterdir()])
    
    @staticmethod
    def extract_gene_counts(count_paths: list[Path]) -> pd.DataFrame:
        gene_counts = {}
        for count_path in count_paths:
            sra = count_path.stem.replace("_salmon", "")
            with count_path.open() as inhandle:
                reader_iterator = reader(inhandle, delimiter="\t")
                header = next(reader_iterator)
                gene_counts[sra] = {line[0]: float(line[3]) for line in reader_iterator}
        return pd.DataFrame(gene_counts)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-counts_dir", type=str, required=True)
    parser.add_argument("-outdir", type=str, required=True)
    args = parser.parse_args()

    sc = SalmonCountCollator(Path(args.counts_dir), Path(args.outdir))
    sc.run()