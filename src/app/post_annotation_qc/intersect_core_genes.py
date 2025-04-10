from argparse import ArgumentParser
from csv import DictReader
import matplotlib.pyplot as plt # type: ignore
from matplotlib_venn import venn2
from pathlib import Path

class GeneIntersector:
    def __init__(self, cds_metadata: Path, ncbi_genes_path: Path) -> None:
        self.cds_metadata = cds_metadata
        self.ncbi_genes_path = ncbi_genes_path

    def run(self) -> None:
        # Extract genes from files
        tatat_core_genes = self.extract_tatat_core_genes(self.cds_metadata)
        ncbi_genes = self.extract_ncbi_genes(self.ncbi_genes_path)

        # Remove hypothetical genes, distinguished by "LOC" prefix
        tatat_core_genes = self.remove_loc_genes(tatat_core_genes)
        ncbi_genes = self.remove_loc_genes(ncbi_genes)

        # Remove low confidence genes, distinguished by "CUNH" prefix
        tatat_core_genes = self.remove_cunh_genes(tatat_core_genes)
        ncbi_genes = self.remove_cunh_genes(ncbi_genes)

        # Generate Venn Diagram showing overlap of gene sets
        self.generate_venn_diagram(tatat_core_genes, ncbi_genes, self.ncbi_genes_path.parent)

        # Write NCBI genes missing from TATAT genes to file
        missing_ncbi_genes = ncbi_genes - tatat_core_genes
        outfile = self.ncbi_genes_path.parent / "missing_ncbi_genes.txt"
        with outfile.open("w") as outhandle:
            for gene in missing_ncbi_genes:
                outhandle.write(f"{gene}\n")

    @staticmethod
    def extract_tatat_core_genes(data_path: Path) -> set[str]:
        core_genes = set()
        with data_path.open() as inhandle:
            data_reader = DictReader(inhandle)
            for data in data_reader:
                if data["core_cds"] != "True":
                    continue
                core_genes.add(data["gene"])
        return core_genes

    @staticmethod
    def extract_ncbi_genes(data_path: Path) -> set[str]:
        core_genes = set()
        with data_path.open() as inhandle:
            data_reader = DictReader(inhandle, delimiter="\t")
            for data in data_reader:
                core_genes.add(data["Symbol"])
        return core_genes

    @staticmethod
    def remove_loc_genes(genes: set[str]) -> set[str]:
        return {gene for gene in genes if not gene.startswith("LOC")}

    @staticmethod
    def remove_cunh_genes(genes: set[str]) -> set[str]:
        return {gene for gene in genes if not gene.startswith("CUNH")}

    @staticmethod
    def generate_venn_diagram(tatat_core_genes: set[str], ncbi_genes: set[str], outdir: Path) -> None:
        venn_diagram = venn2([tatat_core_genes, ncbi_genes], ("TATAT Genes", "NCBI Genes"))
        venn_diagram.get_label_by_id("A").set_position((-.7,0))
        venn_diagram.get_label_by_id("B").set_position((.7,0))
        out_plot = outdir / "venn2.png"
        plt.savefig(out_plot, bbox_inches="tight")
        plt.close()

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-cds_metadata", type=str, required=True)
    parser.add_argument("-ncbi_genes_path", type=str, required=True)
    args = parser.parse_args()

    ga = GeneIntersector(Path(args.cds_metadata), Path(args.ncbi_genes_path))
    ga.run()