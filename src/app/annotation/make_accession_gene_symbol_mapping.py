from argparse import ArgumentParser
from csv import reader
from itertools import islice
from json import loads
from pathlib import Path
import sqlite3
import subprocess
from typing import Iterator

class AccessionGeneMapper:
    def __init__(self, blast_results: Path, sqlite_db: Path, table_name: str) -> None:
        self.blast_results = blast_results
        self.sqlite_db = sqlite_db
        self.table_name = table_name

    def run(self, rna_type: str, upper: bool=True) -> None:
        # Extract non redundant accession numbers for submission to NCBI
        accession_numbers = self.extract_accession_numbers(self.blast_results)
        print(f"\nAll accession numbers count: {len(accession_numbers)}")

        accession_numbers_gene_symbol_mapping = {}
        batch_sizes = [500, 500]
        for batch_size in batch_sizes:
            remaining_accession_numbers = {acc for acc in accession_numbers if acc not in accession_numbers_gene_symbol_mapping}
            accession_numbers_gene_symbol_mapping.update(
                self.batch_ncbi_datasets_accession_gene_mapping(remaining_accession_numbers, rna_type, batch_size)
                )
            print(f"Mapping keys so far: {len(accession_numbers_gene_symbol_mapping)}")

        accession_numbers_gene_symbol_mapping = self.remove_extraneous_accession_numbers(accession_numbers_gene_symbol_mapping, accession_numbers)
        if upper:
            accession_numbers_gene_symbol_mapping = self.upper_case_genes(accession_numbers_gene_symbol_mapping)

        with sqlite3.connect(self.sqlite_db) as connection:
            self.insert_accession_gene_mapping_into_table(connection, accession_numbers_gene_symbol_mapping, self.table_name)

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
    def batch_ncbi_datasets_accession_gene_mapping(cls, accession_numbers: set[str], rna_type, batch_size: int=500) -> dict[str]:
        print("Beginning NCBI Datasets batches")
        accession_number_gene_symbol_mapping = {}
        for i, accession_numbers_batch in enumerate(cls.chunk_set(accession_numbers, batch_size), 1):
            accession_number_gene_symbol_mapping.update(
                cls.submit_accession_numbers_with_ncbi_datasets(accession_numbers_batch, rna_type)
                )
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
    def submit_accession_numbers_with_ncbi_datasets(accession_numbers: set[str], rna_type: str, quiet: bool=True) -> dict[str]:
        if rna_type == "coding":
            acceptable_gene_types = ["PROTEIN_CODING"]
        elif rna_type == "non_coding":
            acceptable_gene_types = ["ncRNA", "rRNA", "snRNA", "snoRNA", "PSEUDO"]
        else:
            print(f"rna_type '{rna_type}' not recognized. Terminating script")
            quit(1)

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
                        gene_symbol = record["gene"]["synonyms"][0] # Better to get a real gene symbol if possible
                    except KeyError:
                        pass

                try:
                    gene_type = record["gene"]["type"]
                except KeyError:
                    continue # Some annotations lack this information, making it unreliable

                if gene_type not in acceptable_gene_types:
                    continue

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
    def insert_accession_gene_mapping_into_table(connection: sqlite3.Connection,
                                                 accession_numbers_gene_symbol_mapping: dict[str],
                                                 table_name: str) -> None:
        cursor = connection.cursor()
        values = [(acc, gene) for acc, gene in accession_numbers_gene_symbol_mapping.items()]
        sql_statement = f"INSERT INTO {table_name} VALUES (?, ?)"
        cursor.executemany(sql_statement, values)
        connection.commit()

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-blast_results", type=str, required=True)
    parser.add_argument("-sqlite_db", type=str, required=True)
    parser.add_argument("-table_name", type=str, required=True)
    parser.add_argument("-rna_type", type=str, required=True)
    args = parser.parse_args()

    agm = AccessionGeneMapper(Path(args.blast_results), Path(args.sqlite_db), args.table_name)
    agm.run(args.rna_type)