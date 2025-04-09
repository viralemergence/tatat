from argparse import ArgumentParser
from collections import defaultdict
from csv import DictReader
import matplotlib.pyplot as plt # type: ignore
import pandas as pd # type: ignore
from pathlib import Path
from seaborn import violinplot # type: ignore
from typing import Any

class ExpressionManager:
    def __init__(self, ncbi_genes_path: Path, tissue_expression_path: Path, missing_genes_path: Path) -> None:
        self.ncbi_genes_path = ncbi_genes_path
        self.tissue_expression_path = tissue_expression_path
        self.missing_genes_path = missing_genes_path

    def run(self) -> None:
        # Prepare data per gene with highest expression level and corresponding tissue
        ncbi_genes = self.extract_ncbi_genes(self.ncbi_genes_path)
        top_tissue_per_gene = self.extract_top_tissue_per_gene(self.tissue_expression_path, ncbi_genes)

        # Add field to data indicating if gene is one of the transcriptome's missing genes
        missing_genes = self.extract_missing_genes(self.missing_genes_path)
        top_tissue_per_gene = self.append_missing_gene_field(top_tissue_per_gene, missing_genes)

        # Calculate count of missing genes per tissue and make pie chart
        top_tissue_count = self.count_top_tissues(top_tissue_per_gene)
        self.write_top_tissue_count(self.tissue_expression_path.parent, top_tissue_count)
        self.generate_pie_chart(top_tissue_count, self.tissue_expression_path.parent)

        # Generate violin plot of missing genes for testes
        self.generate_violin_plot(top_tissue_per_gene, self.tissue_expression_path.parent)

    @staticmethod
    def extract_ncbi_genes(data_path: Path) -> set[str]:
        core_genes = set()
        with data_path.open() as inhandle:
            data_reader = DictReader(inhandle, delimiter="\t")
            for data in data_reader:
                core_genes.add(data["Symbol"])
        return core_genes

    @staticmethod
    def extract_top_tissue_per_gene(data_path: Path, ncbi_genes: set[str]) -> list[dict[Any]]:
        top_tissue_per_gene = list()

        descriptive_fields = {"Name", "Description"}
        with data_path.open() as inhandle:
            data_reader = DictReader(inhandle, delimiter="\t")
            for data in data_reader:
                gene = data["Description"]
                if gene not in ncbi_genes:
                    continue
                expression_levels = {k: float(v) for k, v in data.items() if k not in descriptive_fields}
                highest_expressed_tissue = max(expression_levels, key=expression_levels.get)
                highest_expression = expression_levels[highest_expressed_tissue]
                if highest_expression == 0:
                    continue

                extracted_data = {"Gene": gene, "Tissue": highest_expressed_tissue, "Expression": highest_expression}
                top_tissue_per_gene.append(extracted_data)
        return top_tissue_per_gene

    @staticmethod
    def extract_missing_genes(missing_genes_path: Path) -> set[str]:
        missing_genes = set()
        with missing_genes_path.open() as inhandle:
            for line in inhandle:
                gene = line.strip()
                missing_genes.add(gene)
        return missing_genes

    @staticmethod
    def append_missing_gene_field(top_tissue_per_gene: list[dict[Any]], missing_genes: set[str]) -> list[dict[Any]]:
        appended_data = list()
        for data in top_tissue_per_gene:
            gene = data["Gene"]
            if gene in missing_genes:
                missing_gene = True
            else:
                missing_gene = False
            data.update({"Missing Genes": missing_gene})
            appended_data.append(data)
        return appended_data

    @staticmethod
    def count_top_tissues(top_tissue_per_gene: list[dict[Any]]) -> dict[int]:
        top_tissue_count = defaultdict(int)
        for data in top_tissue_per_gene:
            if not data["Missing Genes"]:
                continue
            tissue = data["Tissue"]
            top_tissue_count[tissue] += 1
        top_tissue_count = dict(sorted(top_tissue_count.items(), key=lambda item: item[1], reverse=True))
        return top_tissue_count

    @staticmethod
    def write_top_tissue_count(outdir: Path, top_tissue_count: dict[int]) -> None:
        out_path = outdir / "top_tissue_count.txt"
        with out_path.open("w") as outhandle:
            for tissue, count in top_tissue_count.items():
                outhandle.write(f"{tissue},{count}\n")

    @staticmethod
    def generate_pie_chart(top_tissue_count: dict[int], outdir: Path) -> None:
        tissues, counts = zip(*[(k, v) for k, v in top_tissue_count.items()])
        explode = [0.1] + [0]*(len(tissues)-1)
        labels = list(tissues[:10]) + [""]*(len(tissues)-10)
        labels = [label if len(label) < 30 else f"{label[:20]}..." for label in labels]
        autopct = lambda v: f"{v:.1f}%" if v > 2.3 else None

        plt.pie(counts, labels=labels, explode=explode, autopct=autopct, pctdistance=0.75)
        plt.title("Missing Genes By Tissue")
        out_plot = outdir / "pie.png"
        plt.savefig(out_plot, bbox_inches="tight")
        plt.close()

    @staticmethod
    def generate_violin_plot(top_tissue_per_gene: list[dict[Any]], outdir: Path) -> None:
        top_tissue_per_gene_df = pd.DataFrame(top_tissue_per_gene)

        top_tissue_per_gene_df = top_tissue_per_gene_df[top_tissue_per_gene_df["Tissue"].isin(["Testis"])]
        medians = top_tissue_per_gene_df.groupby(["Missing Genes"])["Expression"].median()
        lower_quartiles = top_tissue_per_gene_df.groupby(["Missing Genes"])["Expression"].quantile(q=0.25)
        print("Medians:")
        print(medians)
        print("\nLower Quartiles")
        print(lower_quartiles)

        vplot = violinplot(data=top_tissue_per_gene_df, x="Tissue", y="Expression", hue="Missing Genes", log_scale=True, bw_adjust=.5)
        vplot.axhline(y=medians[False], color="lightgreen", linestyle="--")
        out_plot = outdir / "violin.png"
        plt.savefig(out_plot)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-ncbi_genes_path", type=str, required=True)
    parser.add_argument("-tissue_expression_path", type=str, required=True)
    parser.add_argument("-missing_genes_path", type=str, required=True)
    args = parser.parse_args()

    em = ExpressionManager(Path(args.ncbi_genes_path), Path(args.tissue_expression_path), Path(args.missing_genes_path))
    print()
    em.run()