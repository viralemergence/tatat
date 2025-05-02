from argparse import ArgumentParser
from pathlib import Path
import re
from typing import Iterator

class NcbiGeneExtractor:
    def __init__(self, ncbi_cds_fasta: Path, outpath: Path) -> None:
        self.ncbi_cds_fasta = ncbi_cds_fasta
        self.outpath = outpath

    def run(self) -> None:
        genes = self.extract_gene_symbols(self.ncbi_cds_fasta)
        with self.outpath.open("w") as outhandle:
            for gene in genes:
                outhandle.write(f"{gene}\n")

    @classmethod
    def extract_gene_symbols(cls, ncbi_cds_fasta: Path) -> set[str]:
        genes = set()

        gene_pattern = r"\[gene=(.*?)\]"

        for fasta_seq in cls.fasta_chunker(ncbi_cds_fasta):
            header = fasta_seq[0]
            gene = re.search(gene_pattern, header).group()[6:-1]
            genes.add(gene)
        return genes

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

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-ncbi_cds_fasta", type=str, required=True)
    parser.add_argument("-outpath", type=str, required=True)
    args = parser.parse_args()

    nge = NcbiGeneExtractor(Path(args.ncbi_cds_fasta), Path(args.outpath))
    nge.run()