from argparse import ArgumentParser
from csv import DictReader
from pathlib import Path

class GeneIntersector:
    def __init__(self, cds_metadata: Path, ncbi_genes_path: Path) -> None:
        self.cds_metadata = cds_metadata
        self.ncbi_genes_path = ncbi_genes_path

    def run(self) -> None:
        tatat_core_genes = self.extract_tatat_core_genes(self.cds_metadata)
        ncbi_genes = self.extract_ncbi_genes(self.ncbi_genes_path)

        print(f"TATAT genes: {len(tatat_core_genes)}")
        print(f"NCBI genes: {len(ncbi_genes)}")
        print(len(tatat_core_genes & ncbi_genes))
        print()

        tatat_core_genes = self.remove_loc_genes(tatat_core_genes)
        ncbi_genes = self.remove_loc_genes(ncbi_genes)

        tatat_core_genes = self.remove_cunh_genes(tatat_core_genes)
        ncbi_genes = self.remove_cunh_genes(ncbi_genes)
        print(f"TATAT genes: {len(tatat_core_genes)}")
        print(f"NCBI genes: {len(ncbi_genes)}")
        print(len(tatat_core_genes & ncbi_genes))
        print()

        # print(ncbi_genes - tatat_core_genes)

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

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-cds_metadata", type=str, required=True)
    parser.add_argument("-ncbi_genes_path", type=str, required=True)
    args = parser.parse_args()

    ga = GeneIntersector(Path(args.cds_metadata), Path(args.ncbi_genes_path))
    ga.run()