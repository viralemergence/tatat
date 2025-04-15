from argparse import ArgumentParser
from Bio.Align import PairwiseAligner, substitution_matrices # type: ignore
from collections import defaultdict
from copy import deepcopy
import matplotlib.pyplot as plt # type: ignore
import pandas as pd # type: ignore
from pathlib import Path
import re
from seaborn import jointplot # type: ignore
from typing import Iterator

class PairwiseGeneAligner:
    def __init__(self, tatat_aa_fasta: Path, ncbi_cds_fasta: Path, ncbi_aa_fasta: Path) -> None:
        self.tatat_aa_fasta = tatat_aa_fasta
        self.ncbi_cds_fasta = ncbi_cds_fasta
        self.ncbi_aa_fasta = ncbi_aa_fasta

    def run(self) -> None:
        # Extract all aa sequences needed for analysis
        tatat_aa_sequences = self.extract_tatat_aa_sequences(self.tatat_aa_fasta)
        longest_ncbi_genes = self.extract_longest_ncbi_genes(self.ncbi_cds_fasta)
        ncbi_aa_sequences = self.extract_ncbi_aa_sequences(self.ncbi_aa_fasta, longest_ncbi_genes)

        # Perform pairwise alignments and return as pd.Dataframe
        print("\nStarting pairwise alignments\n(This may take awhile)")
        tatat_ncbi_alignment_data = self.calculate_tatat_ncbi_alignment_data(tatat_aa_sequences, ncbi_aa_sequences)

        # Count the number of alignments falling into specific sequence identity bins, then convert to proportions
        bin_counts = self.bin_sequence_identities(tatat_ncbi_alignment_data)
        print(bin_counts)

        print("\nStarting graphing")
        outdir = Path("/src/data/ncbi_gene_comparison")

        # Make plot for sequence identity bin counts
        self.plot_bin_counts(bin_counts, outdir)

        # Make bivariate histogram for gene lengths
        self.plot_sequence_lengths(tatat_ncbi_alignment_data, outdir)

    @classmethod
    def extract_tatat_aa_sequences(cls, tatat_aa_fasta: Path) -> dict[str]:
        return {fasta_seq[0].split(";")[1]: "".join(fasta_seq[1:]).replace("*", "")
                for fasta_seq in cls.fasta_chunker(tatat_aa_fasta)}

    @classmethod
    def extract_longest_ncbi_genes(cls, ncbi_cds_fasta: Path) -> dict[str]:
        longest_genes = defaultdict(lambda: defaultdict(int))

        gene_pattern = r"\[gene=(.*?)\]"
        protein_accession_pattern = r"\[protein_id=(.*?)\]"

        for fasta_seq in cls.fasta_chunker(ncbi_cds_fasta):
            header = fasta_seq[0]
            gene = re.search(gene_pattern, header).group()[6:-1]
            try:
                protein_accession = re.search(protein_accession_pattern, header).group()[12:-1]
            except AttributeError:
                continue
            sequence_len = len("".join(fasta_seq[1:]))

            longest_gene_len = longest_genes[gene]["length"]
            if sequence_len > longest_gene_len:
                new_data = {"length": sequence_len,
                            "protein_accession": protein_accession}
                longest_genes[gene].update(new_data)

        return {data["protein_accession"]: gene for gene, data in longest_genes.items()}
            
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

    @classmethod
    def extract_ncbi_aa_sequences(cls, ncbi_aa_fasta: Path, ncbi_genes: dict[str]) -> dict[str]:
        aa_sequences = dict()

        for fasta_seq in cls.fasta_chunker(ncbi_aa_fasta):
            header = fasta_seq[0]
            protein_accession = header.split(" ")[0][1:]
            try:
                gene = ncbi_genes[protein_accession]
            except KeyError:
                continue
            aa_seq = "".join(fasta_seq[1:])
            aa_sequences[gene] = aa_seq
        return aa_sequences

    @classmethod
    def calculate_tatat_ncbi_alignment_data (cls, tatat_aa_sequences: dict[str], ncbi_aa_sequences: dict[str]) -> pd.DataFrame:
        tatat_ncbi_alignment_data = list()
        baseline_aligner = cls.set_baseline_aligner()

        for tatat_gene, tatat_seq in tatat_aa_sequences.items():
            if tatat_gene.startswith("LOC"):
                continue
            try:
                ncbi_seq = ncbi_aa_sequences[tatat_gene]
            except KeyError:
                continue

            aligner = deepcopy(baseline_aligner)
            try:
                alignments = aligner.align(tatat_seq, ncbi_seq)
            except ValueError:
                continue

            alignment = alignments[0]

            match = 0
            for aa_aligned in zip(alignment[0], alignment[1]):
                if aa_aligned[0] == aa_aligned[1]:
                    match += 1

            len_similarity = round(len(tatat_seq)/ len(ncbi_seq), 2)
            if len_similarity > 0.8:
                ncbi_match_close_len = round((match/len(ncbi_seq))*100, 2)
            else:
                ncbi_match_close_len = None

            data = {"gene": tatat_gene, "match_count": match,
                    "tatat_len": len(tatat_seq), "ncbi_len": len(ncbi_seq),
                    "tatat_match": round((match/len(tatat_seq))*100, 2),
                    "ncbi_match": round((match/len(ncbi_seq))*100, 2),
                    "ncbi_match_close_len": ncbi_match_close_len
                    }
            tatat_ncbi_alignment_data.append(data)
        return pd.DataFrame(tatat_ncbi_alignment_data)

    @staticmethod
    def set_baseline_aligner() -> PairwiseAligner:
        aligner = PairwiseAligner()
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        return aligner

    @staticmethod
    def bin_sequence_identities(tatat_ncbi_alignment_data: pd.DataFrame) -> pd.DataFrame:
        bins = [0, 20, 40, 60, 80, 100]
        columns = ["tatat_match", "ncbi_match", "ncbi_match_close_len"]
        bin_counts = tatat_ncbi_alignment_data[columns].apply(lambda x: round(x.value_counts(bins=bins, sort=False, normalize=True), 2))
        bin_counts["ncbi_match_close_len"] = round(bin_counts["ncbi_match_close_len"] / bin_counts["ncbi_match_close_len"].sum(), 2)
        new_indices = ["0 - 20%", "20 - 40%", "40 - 60%", "60 - 80%", "80 - 100%"]
        bin_counts.index = new_indices
        return bin_counts

    @staticmethod
    def plot_bin_counts(bin_counts: pd.DataFrame, outdir: Path) -> None:
        bplot = bin_counts.T.plot.bar(stacked=True)

        for c in bplot.containers:
            labels = [round(v.get_height(),2) if v.get_height() > 0.6 else '' for v in c]
            bplot.bar_label(c, labels=labels, label_type='center')

        bplot.set_ylabel("Proportion of Sequence Alignments")
        bplot.set_xlabel("Sequence Identity Denominator")
        bplot.legend(title="Sequence Identity", bbox_to_anchor=(1.05, 1), loc="upper left")
        x_labels = ["TATAT", "NCBI", "NCBI\n(80% Length\nSimilarity)"]
        bplot.set_xticklabels(x_labels, rotation=45, ha="center")

        out_plot = outdir / "global_alignment_sequence_identities_bar.png"
        plt.savefig(out_plot, bbox_inches="tight")
        plt.close()

    @staticmethod
    def plot_sequence_lengths(tatat_ncbi_alignment_data: pd.DataFrame, outdir: Path) -> None:
        bin_number = 15
        joint = jointplot(data=tatat_ncbi_alignment_data,
                          x="tatat_len", y="ncbi_len", kind="hist",
                          bins=bin_number, marginal_kws={"bins": bin_number},
                          log_scale=True)
        joint.ax_joint.plot([100, 20_000], [100, 20_000], ls="--", color="black")
        joint.set_axis_labels("TATAT Sequence Lengths", "NCBI Sequence Lengths")
        out_plot = outdir / "gene_lengths_joint.png"
        plt.savefig(out_plot)
        plt.close()

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-tatat_aa_fasta", type=str, required=True)
    parser.add_argument("-ncbi_cds_fasta", type=str, required=True)
    parser.add_argument("-ncbi_aa_fasta", type=str, required=True)
    args = parser.parse_args()

    pga = PairwiseGeneAligner(Path(args.tatat_aa_fasta), Path(args.ncbi_cds_fasta), Path(args.ncbi_aa_fasta))
    pga.run()
    print("Finished\n")