import argparse
import pandas as pd # type: ignore
from pathlib import Path
import seaborn as sns # type: ignore
import matplotlib.pyplot as plt # type: ignore
from sklearn.manifold import MDS # type: ignore
from sklearn.metrics import pairwise_distances # type: ignore
import sqlite3

class SalmonCountMDS:
    def __init__(self, counts: Path, sqlite_db: Path, outdir: Path) -> None:
        self.counts = pd.read_csv(counts, index_col=0)
        self.metadata = self.extract_sample_metadata(sqlite_db)
        self.outdir = outdir

    @staticmethod
    def extract_sample_metadata(sqlite_db: Path) -> pd.DataFrame:
        with sqlite3.connect(sqlite_db) as connection:
            cursor = connection.cursor()
            cursor.execute("SELECT * FROM samples")
            results = cursor.fetchall()
            columns = [description[0] for description in cursor.description]
            return pd.DataFrame(results, columns=columns)

    def run(self) -> None:
        # Filter gene counts removing genes with essentially no count,
        # or housekeeping genes with overwhelming counts
        gene_counts = self.filter_gene_counts(self.counts)

        # Calculate MDS data for subsequent plotting and add metadata
        df = self.calculate_mds_data(gene_counts)
        df["Tissue"] = self.metadata["tissue"]
        df["Gender"] = self.metadata["gender"]

        # Generate MDS plot and save
        self.generate_mds_graph(df, self.outdir)

    @staticmethod
    def filter_gene_counts(gene_counts: pd.DataFrame) -> pd.DataFrame:
        gene_counts = gene_counts[gene_counts.var(axis=1) > 0]
        gene_counts = gene_counts[gene_counts.sum(axis=1) >= 100]
        gene_counts = gene_counts[gene_counts.mean(axis=1) <= 10_000]
        return gene_counts

    @staticmethod
    def calculate_mds_data(gene_counts: pd.DataFrame) -> pd.DataFrame:
        distances = pairwise_distances(gene_counts.T)
        mds = MDS(n_components=2, dissimilarity="precomputed", random_state=0)
        data_transformed = mds.fit_transform(distances)
        return pd.DataFrame(data_transformed, columns=["D1", "D2"])

    @classmethod
    def generate_mds_graph(cls, df: pd.DataFrame, outdir: Path) -> None:
        plt.rcParams["svg.fonttype"] = "none"
        fig, ax = plt.subplots(figsize=(3, 3))
        scatter = sns.scatterplot(x="D1", y="D2", data=df, style="Tissue", hue="Tissue", s=100)
        ax.ticklabel_format(axis="x", style="scientific", scilimits=(0,0), useMathText=True)
        ax.ticklabel_format(axis="y", style="scientific", scilimits=(0,0), useMathText=True)

        sns.move_legend(scatter, "upper left", bbox_to_anchor=(1,1))

        cls.add_tissue_connecting_lines(df)
        cls.add_gender_symbol(scatter, df)

        outpath = outdir / "MDS.svg"
        fig.savefig(outpath, bbox_inches="tight", dpi=300)
        plt.close(fig)

    @staticmethod
    def add_tissue_connecting_lines(df: pd.DataFrame) -> None:
        for tissue in df["Tissue"].unique():
            if tissue == "TT":
                continue
            if tissue == "OV":
                x = list(df.loc[(df["Tissue"] == tissue) | (df["Tissue"] == "TT"), "D1"])
                y = list(df.loc[(df["Tissue"] == tissue) | (df["Tissue"] == "TT"), "D2"])
            else:
                x = list(df.loc[df["Tissue"] == tissue, "D1"])
                y = list(df.loc[df["Tissue"] == tissue, "D2"])
            plt.plot([x[0], x[1]], [y[0], y[1]], color="grey", alpha=0.5, zorder=-1)

    @staticmethod
    def add_gender_symbol(scatter, df: pd.DataFrame) -> None:
        spacer = 1000
        for _, datum in df.iterrows():
            if datum["Gender"] == "male":
                text = "\u2642"
            if datum["Gender"] == "female":
                text = "\u2640"
            scatter.text(datum["D1"]+(spacer/2), datum["D2"]-spacer, text, fontsize=10, fontweight="bold")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-counts", type=str, required=True)
    parser.add_argument("-sqlite_db", type=str, required=True)
    parser.add_argument("-outdir", type=str, required=True)
    args = parser.parse_args()

    scp = SalmonCountMDS(Path(args.counts), Path(args.sqlite_db), Path(args.outdir))
    scp.run()