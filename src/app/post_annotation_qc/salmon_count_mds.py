import argparse
from csv import reader
import pandas as pd # type: ignore
from pathlib import Path
import seaborn as sns # type: ignore
import matplotlib.pyplot as plt # type: ignore
from sklearn.manifold import MDS # type: ignore
from sklearn.metrics import pairwise_distances # type: ignore

class SalmonCountMDS:
    def __init__(self, counts: Path, metadata: Path) -> None:
        self.counts = pd.read_csv(counts, index_col=0)
        self.metadata = self.extract_metadata(metadata)

    @staticmethod
    def extract_metadata(metadata: Path) -> dict[str]:
        with metadata.open() as inhandle:
            reader_iterator = reader(inhandle, delimiter=",")
            header = next(reader_iterator)
            return {line[0]: line[2] for line in reader_iterator}

    def run(self) -> None:
        gene_counts = self.counts

        gene_counts = gene_counts[gene_counts.var(axis=1) > 0]
        gene_counts = gene_counts[gene_counts.sum(axis=1) >= 100]

        distances = pairwise_distances(gene_counts.T)

        mds = MDS(n_components=2, dissimilarity="precomputed", random_state=0)
        data_transformed = mds.fit_transform(distances)

        df = pd.DataFrame(data_transformed, columns=["D1", "D2"])

        df["Tissue"] = self.metadata.values()
        
        fig, ax = plt.subplots(figsize=(3, 3))
        scatter = sns.scatterplot(x="D1", y="D2", data=df, style="Tissue", hue="Tissue")
        sns.move_legend(scatter, "upper left", bbox_to_anchor=(1,1))

        fig.savefig("/src/data/mds.png", bbox_inches="tight", dpi=300)
        plt.close(fig)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-counts", type=str, required=True)
    parser.add_argument("-metadata", type=str, required=True)
    args = parser.parse_args()

    scp = SalmonCountMDS(Path(args.counts), Path(args.metadata))
    scp.run()