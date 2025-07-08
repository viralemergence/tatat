from argparse import ArgumentParser
import matplotlib.pyplot as plt # type: ignore
from matplotlib_venn import venn2 # type: ignore
from pathlib import Path
import sqlite3

class GeneIntersector:
    def __init__(self, sqlite_db: Path, ncbi_genes_path: Path, outdir: Path) -> None:
        self.sqlite_db = sqlite_db
        self.ncbi_genes_path = ncbi_genes_path
        self.outdir = outdir

    def run(self) -> None:
        # Extract genes from files
        tatat_core_genes = self.extract_tatat_core_genes(self.sqlite_db)
        ncbi_genes = self.extract_ncbi_genes(self.ncbi_genes_path)

        # Remove hypothetical genes, distinguished by "LOC" prefix
        tatat_core_genes = self.remove_loc_genes(tatat_core_genes)
        ncbi_genes = self.remove_loc_genes(ncbi_genes)

        # Remove low confidence genes, distinguished by "CUNH" prefix
        tatat_core_genes = self.remove_cunh_genes(tatat_core_genes)
        ncbi_genes = self.remove_cunh_genes(ncbi_genes)

        # Generate Venn Diagram showing overlap of gene sets
        self.generate_venn_diagram(tatat_core_genes, ncbi_genes, self.outdir)

        # Write NCBI genes missing from TATAT genes to file
        missing_ncbi_genes = ncbi_genes - tatat_core_genes
        outfile = self.outdir / "missing_ncbi_genes.txt"
        with outfile.open("w") as outhandle:
            for gene in missing_ncbi_genes:
                outhandle.write(f"{gene}\n")

    @staticmethod
    def extract_tatat_core_genes(sqlite_db: Path) -> set[str]:
        with sqlite3.connect(sqlite_db) as connection:
            cursor = connection.cursor()
            cursor.execute("SELECT gene_symbol FROM cds WHERE core_cds=1")
            return {row[0] for row in cursor.fetchall()}

    @staticmethod
    def extract_ncbi_genes(data_path: Path) -> set[str]:
        core_genes = set()
        with data_path.open() as inhandle:
            for line in inhandle:
                core_genes.add(line.strip())
        return core_genes

    @staticmethod
    def remove_loc_genes(genes: set[str]) -> set[str]:
        return {gene for gene in genes if not gene.startswith("LOC")}

    @staticmethod
    def remove_cunh_genes(genes: set[str]) -> set[str]:
        return {gene for gene in genes if not gene.startswith("CUNH")}

    @staticmethod
    def generate_venn_diagram(tatat_core_genes: set[str], ncbi_genes: set[str], outdir: Path) -> None:
        plt.rcParams["svg.fonttype"] = "none"
        venn_diagram = venn2([tatat_core_genes, ncbi_genes], ("TATAT Genes", "NCBI Genes"),
                             set_colors=((3/256, 10/256, 167/256), (103/256, 146/256, 103/256)), alpha=0.75)
        # Cobalt Blue (3/256, 10/256, 167/256)
        # Terre Verte (103/256, 146/256, 103/256)
        venn_diagram.get_label_by_id("A").set_position((-0.6,0.2))
        venn_diagram.get_label_by_id("B").set_position((0.6,0.2))

        venn_diagram.get_label_by_id("10").set_position((-0.7,-0.1))
        venn_diagram.get_label_by_id("01").set_position((0.7,-0.1))

        out_plot = outdir / "gene_intersection_venn.svg"
        plt.savefig(out_plot, bbox_inches="tight")
        plt.close()

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-sqlite_db", type=str, required=True)
    parser.add_argument("-ncbi_genes_path", type=str, required=True)
    parser.add_argument("-outdir", type=str, required=True)
    args = parser.parse_args()

    ga = GeneIntersector(Path(args.sqlite_db), Path(args.ncbi_genes_path), Path(args.outdir))
    ga.run()