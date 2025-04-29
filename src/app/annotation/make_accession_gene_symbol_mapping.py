from argparse import ArgumentParser
from csv import reader
from itertools import islice
from json import loads
from pathlib import Path
import subprocess
from typing import Iterator

class AccessionGeneMapper:
    def __init__(self, blast_results: Path, datasets_mapping_path: Path) -> None:
        self.blast_results = blast_results
        self.datasets_mapping_path = datasets_mapping_path

    def run(self, upper: bool=True) -> None:
        # Extract non redundant accession numbers for submission to NCBI
        accession_numbers = self.extract_accession_numbers(self.blast_results)
        print(f"\nAll accession numbers count: {len(accession_numbers)}")

        if not self.datasets_mapping_path.exists():
            accession_numbers_gene_symbol_mapping = {}
            batch_sizes = [500, 500]
            for batch_size in batch_sizes:
                remaining_accession_numbers = {acc for acc in accession_numbers if acc not in accession_numbers_gene_symbol_mapping}
                accession_numbers_gene_symbol_mapping.update(self.batch_ncbi_datasets_accession_gene_mapping(remaining_accession_numbers, batch_size))
                print(f"Mapping keys so far: {len(accession_numbers_gene_symbol_mapping)}")

            accession_numbers_gene_symbol_mapping = self.remove_extraneous_accession_numbers(accession_numbers_gene_symbol_mapping, accession_numbers)
            if upper:
                accession_numbers_gene_symbol_mapping = self.upper_case_genes(accession_numbers_gene_symbol_mapping)

            self.write_datasets_mapping(accession_numbers_gene_symbol_mapping, self.datasets_mapping_path)
        else:
            accession_numbers_gene_symbol_mapping = self.extract_datasets_mapping(self.datasets_mapping_path)
        print(f"Datasets mapping keys count: {len(accession_numbers_gene_symbol_mapping)}")

    @staticmethod
    def extract_accession_numbers(blast_results: Path) -> set[str]:
        accession_numbers = set()
        with blast_results.open() as inhandle:
            blast_reader = reader(inhandle, delimiter="\t")
            for line in blast_reader:
                accession_number = line[1]
                accession_numbers.add(accession_number)
        return accession_numbers

    @classmethod
    def batch_ncbi_datasets_accession_gene_mapping(cls, accession_numbers: set[str], batch_size: int=500) -> dict[str]:
        print("Beginning NCBI Datasets batches")
        accession_number_gene_symbol_mapping = {}
        for i, accession_numbers_batch in enumerate(cls.chunk_set(accession_numbers, batch_size), 1):
            accession_number_gene_symbol_mapping.update(cls.submit_accession_numbers_with_ncbi_datasets(accession_numbers_batch))
            if i % 100 == 0:
                print(f"NCBI Datasets batch {i} finished")
        print(f"\nNCBI Datasets batches complete\n")
        return accession_number_gene_symbol_mapping

    @staticmethod
    def chunk_set(data: set, chunk_size: int) -> Iterator[set]:
        iterable = iter(data)
        while True:
            chunk = set(islice(iterable, chunk_size))
            if not chunk:
                break
            yield chunk

    @staticmethod
    def submit_accession_numbers_with_ncbi_datasets(accession_numbers: set[str], quiet: bool=True) -> dict[str]:
        if not quiet:
            print("Starting NCBI Datasets submission\n(This may take awhile)")
        datasets_command = ["datasets", "summary", "gene", "accession"] + list(accession_numbers)

        p1 = subprocess.Popen(datasets_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        accession_numbers_gene_mapping = {}

        while p1.poll() is None and (line := p1.stdout.readline()) != "":
            data = loads(line)
            try:
                records = data["reports"]
            except KeyError:
                continue # On iterations sometimes batches return empty jsons

            for record in records:
                try:
                    accession_number = record["query"][0]
                except KeyError:
                    continue # Can't use data without original accession number

                try:
                    chromosomes = record["gene"]["chromosomes"]
                    if len(chromosomes) == 1 and chromosomes[0] == "MT":
                        continue # Currently mitochondrial hits tend to include the whole genome and are uninformative
                except KeyError:
                    pass

                gene_symbol = record["gene"]["symbol"]

                if gene_symbol.startswith("LOC") or gene_symbol.startswith("CUN"):
                    try:
                        gene_symbol = record["gene"]["synonyms"][0]
                    except KeyError:
                        pass

                accession_numbers_gene_mapping[accession_number] = gene_symbol
        p1.wait()

        if p1.poll() != 0:
            print(p1.stderr.readlines())
            raise Exception("Datasets submission did not complete successfully")
        if not quiet:
            print("NCBI Submission complete")
        return accession_numbers_gene_mapping

    @staticmethod
    def remove_extraneous_accession_numbers(accession_numbers_gene_mapping: dict[str], accession_numbers: set[str]):
        return {k: v for k, v in accession_numbers_gene_mapping.items() if k in accession_numbers}

    @staticmethod
    def upper_case_genes(accession_numbers_gene_mapping: dict[str]):
        return {k: v.upper() for k, v in accession_numbers_gene_mapping.items()}

    @staticmethod
    def write_datasets_mapping(accession_gene_mapping: dict[str], outpath: Path) -> None:
        keys = sorted(list(accession_gene_mapping.keys()))
        with outpath.open("w") as outhandle:
            for key in keys:
                value = accession_gene_mapping[key]
                outhandle.write(f"{key},{value}\n")

    @staticmethod
    def extract_datasets_mapping(inpath: Path) -> dict[str]:
        accession_numbers_gene_symbol_mapping = {}
        with inpath.open() as inhandle:
            for line in inhandle:
                line_info = line.strip().split(",")
                accession_numbers_gene_symbol_mapping[line_info[0]] = line_info[1]
        return accession_numbers_gene_symbol_mapping

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-blast_results", type=str, required=True)
    parser.add_argument("-datasets_mapping", type=str, required=True)
    args = parser.parse_args()

    agm = AccessionGeneMapper(Path(args.blast_results), Path(args.datasets_mapping))
    agm.run()