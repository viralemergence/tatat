from argparse import ArgumentParser
from Bio import Entrez # type: ignore
from csv import reader
from itertools import islice
from json import loads
from os import getenv
from pathlib import Path
import subprocess
from typing import Iterator

class AccessionGeneMapper:
    def __init__(self, blast_results: Path, datasets_mapping_path: Path,
                 accepted_accession_numbers_path: Path,
                 accession_gene_id_mapping_path: Path) -> None:
        self.blast_results = blast_results
        self.datasets_mapping_path = datasets_mapping_path

        self.accepted_accession_numbers_path = accepted_accession_numbers_path
        self.accession_gene_id_mapping_path = accession_gene_id_mapping_path

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

        remaining_accession_numbers = {acc for acc in accession_numbers if acc not in accession_numbers_gene_symbol_mapping}
        print(len(remaining_accession_numbers))
        for acc, gene_symbol in accession_numbers_gene_symbol_mapping.items():
            if gene_symbol.startswith("LOC") or gene_symbol.startswith("CUN"):
                remaining_accession_numbers.add(acc)
        print(len(remaining_accession_numbers))


        Entrez.email = "alexander.brown@wsu.edu"
        Entrez.local_cache = self.accepted_accession_numbers_path.parent
        Entrez.api_key = getenv("NCBI_API_KEY")

        acc_test = list(remaining_accession_numbers)[1]
        acc_test = "XM_016136009.2"
        acc_test2 = "XP_015991495.1"
        # with Entrez.efetch(db="nucleotide", id=acc_test, rettype="gene_table", retmode="text") as handle:
        #     print(handle.read())
        with Entrez.efetch(db="protein", id=acc_test2, rettype="gb", retmode="text") as handle:
            print(handle.read())

        # print(Entrez.read(Entrez.einfo()))
        # with Entrez.elink(dbfrom="nucleotide", id=acc_test, db="protein",
        #                   linkname="nucleotide_protein", idtype="acc") as handle:
        #     for record in Entrez.parse(handle):
        #         print(record)
        return

        Entrez.email = "alexander.brown@wsu.edu"
        Entrez.local_cache = self.accepted_accession_numbers_path.parent
        Entrez.api_key = getenv("NCBI_API_KEY")

        potential_uniprot_accession_numbers = set()
        for acc in accession_numbers:
            if len(acc.split(".")[0]) in [6, 10]:
                if acc[-2] == "_" or acc[0] == "K":
                    continue
                potential_uniprot_accession_numbers.add(acc)
                print(acc)
        print(len(potential_uniprot_accession_numbers))
        return

        if self.accepted_accession_numbers_path.exists():
            print("Accepted accession numbers file detected\nLoading file")
            accepted_accession_numbers = self.extract_accepted_accession_numbers(self.accepted_accession_numbers_path)
        else:
            accepted_accession_numbers = self.batch_efetch_accepted_accession_numbers(accession_numbers)
            self.write_accepted_accession_numbers(accepted_accession_numbers, self.accepted_accession_numbers_path)
        print(f"Accepted accession numbers count: {len(accepted_accession_numbers)}")

        batch_test = list(accepted_accession_numbers)[:2_000]

        #accession_number_entrez_id_mapping = self.batch_elink_accession_gene_id_mapping(batch_test)
        accession_number_entrez_id_mapping = self.batch_elink_accession_gene_id_mapping(accepted_accession_numbers)
        self.write_accession_gene_id_mapping(accession_number_entrez_id_mapping, self.accession_gene_id_mapping_path)

        gene_ids = list(set(accession_number_entrez_id_mapping.values()))
        
        print(len(accession_number_entrez_id_mapping))
        print(len(gene_ids))
        return

        accession_number_entrez_id_mapping = {}
        with Entrez.elink(dbfrom="protein", id=batch_test, db="gene", linkname="protein_gene", idtype="acc") as handle:
            for record in Entrez.parse(handle):
                try:
                    gene_id = record["LinkSetDb"][0]["Link"][0]["Id"]
                except:
                    continue
                accession_numbers = record["IdList"]
                for accession_number in accession_numbers:
                    accession_number_entrez_id_mapping[accession_number] = gene_id

        gene_ids = list(set(accession_number_entrez_id_mapping.values()))
        #print(accession_number_entrez_id_mapping)
        print(len(accession_number_entrez_id_mapping))
        print(len(gene_ids))

        return
        records = Entrez.read(handle)
        handle.close()

        accession_number_entrez_id_mapping = {}
        for record in records:
            print(record)
            quit()
            try:
                gene_id = record["LinkSetDb"][0]["Link"][0]["Id"]
            except:
                continue
            accession_numbers = record["IdList"]
            for accession_number in accession_numbers:
                accession_number_entrez_id_mapping[accession_number] = gene_id

        gene_ids = list(set(accession_number_entrez_id_mapping.values()))
        print(accession_number_entrez_id_mapping)
        print(gene_ids)
        return

        handle = Entrez.efetch(db="protein", id=batch_test, rettype="acc")
        ids = [id.strip() for id in handle]
#        records = Entrez.read(handle)
        handle.close()
        print(ids[:10])
        print(len(ids))
        quit()

        accession_number_query_list = [f"({acc}[ACCN] AND {feature}[FKEY])" for acc in accession_number_set]
        accession_number_query = " OR ".join(accession_number_query_list)

        with Entrez.esearch(db="nucleotide", term=accession_number_query, idtype="acc", retmax=10_000) as handle:
            record = Entrez.read(handle)
        return set(record["IdList"])

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

    @staticmethod
    def extract_accepted_accession_numbers(file_path: Path) -> set[str]:
        accession_numbers = set()
        with file_path.open() as inhandle:
            for line in inhandle:
                accession_number = line.strip()
                accession_numbers.add(accession_number)
        return accession_numbers

    @classmethod
    def batch_efetch_accepted_accession_numbers(cls, accession_numbers: set[str]) -> set[str]:
        accepted_accession_numbers = set()
        chunk_size = 5_000
        for i, accession_numbers_batch in enumerate(cls.chunk_set(accession_numbers, chunk_size), 1):
            accepted_accession_numbers_batch = cls.efetch_accepted_accession_numbers(accession_numbers_batch)
            accepted_accession_numbers.update(accepted_accession_numbers_batch)
            print(f"\033[KEfetch batch {i} finished", end="\r", flush=True)
        print(f"\nEfetch batches complete\nFiltered {len(accession_numbers)} to {len(accepted_accession_numbers)} accesion numbers")
        return accepted_accession_numbers

    @staticmethod
    def efetch_accepted_accession_numbers(accession_numbers: set[str]) -> set[str]:
        with Entrez.efetch(db="protein", id=accession_numbers, rettype="acc") as handle:
            return {acc.strip() for acc in handle if acc.strip()}

    @staticmethod
    def write_accepted_accession_numbers(accession_numbers: set[str], outpath: Path) -> None:
        accession_numbers = sorted(list(accession_numbers))
        with outpath.open("w") as outhandle:
            for acc in accession_numbers:
                outhandle.write(f"{acc}\n")

    @classmethod
    def batch_elink_accession_gene_id_mapping(cls, accession_numbers: set[str]) -> dict[str]:
        print("Beginning elink batches")
        accession_number_entrez_id_mapping = {}
        chunk_size = 2_000
        print(chunk_size)
        for i, accession_numbers_batch in enumerate(cls.chunk_set(accession_numbers, chunk_size), 1):
            accession_number_entrez_id_mapping_batch = cls.elink_accession_gene_id_mapping(accession_numbers_batch)
            accession_number_entrez_id_mapping.update(accession_number_entrez_id_mapping_batch)
            print(f"Elink batch {i} finished")
            #print(f"\033[KElink batch {i} finished", end="\r", flush=True)
        print(f"\nElink batches complete\n")
        return accession_number_entrez_id_mapping

    @staticmethod
    def elink_accession_gene_id_mapping(accession_numbers: set[str]) -> dict[str]:
        accession_number_entrez_id_mapping = {}
        with Entrez.elink(dbfrom="protein", id=accession_numbers, db="gene",
                          linkname="protein_gene", idtype="acc") as handle:
            for record in Entrez.parse(handle):
                try:
                    gene_id = record["LinkSetDb"][0]["Link"][0]["Id"]
                except:
                    continue
                accession_numbers = record["IdList"]
                for accession_number in accession_numbers:
                    accession_number_entrez_id_mapping[accession_number] = gene_id
        return accession_number_entrez_id_mapping

    @staticmethod
    def write_accession_gene_id_mapping(accession_gene_id_mapping: dict[str], outpath: Path) -> None:
        keys = sorted(list(accession_gene_id_mapping.keys()))
        with outpath.open("w") as outhandle:
            for key in keys:
                value = accession_gene_id_mapping[key]
                outhandle.write(f"{key},{value}\n")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-blast_results", type=str, required=True)
    parser.add_argument("-datasets_mapping", type=str, required=True)
    parser.add_argument("-accepted_accession_numbers", type=str, required=True)
    parser.add_argument("-accession_gene_id_mapping", type=str, required=True)
    args = parser.parse_args()

    agm = AccessionGeneMapper(Path(args.blast_results), Path(args.datasets_mapping),
                              Path(args.accepted_accession_numbers),
                              Path(args.accession_gene_id_mapping))
    agm.run()