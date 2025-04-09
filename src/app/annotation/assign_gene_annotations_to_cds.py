from argparse import ArgumentParser
from csv import DictReader, DictWriter, reader
from pathlib import Path
from shutil import copyfileobj
import subprocess
from tempfile import NamedTemporaryFile
from typing import Any

class GeneAssigner:
    def __init__(self, diamond_results: Path, accession_numbers_out: Path, cds_metadata: Path) -> None:
        self.diamond_results = diamond_results
        self.accession_numbers_out = accession_numbers_out
        self.cds_metadata = cds_metadata

    def run(self, upper: bool) -> None:
        # Extract non redundant accession numbers for submission to NCBI
        cds_id_accession_numbers_mapping = self.extract_cds_id_accession_numbers_mapping(self.diamond_results)
        accession_numbers = self.extract_accession_numbers(cds_id_accession_numbers_mapping)
        self.write_accession_numbers(self.accession_numbers_out, accession_numbers)

        # Submit accession numbers to NCBI and make mapping to gene symbols
        accession_numbers_gene_mapping = self.submit_accession_numbers_with_ncbi_datasets(self.accession_numbers_out)
        accession_numbers_gene_mapping = self.remove_extraneous_accession_numbers(accession_numbers_gene_mapping, accession_numbers)
        if upper:
            accession_numbers_gene_mapping = self.upper_case_genes(accession_numbers_gene_mapping)

        # Calculate "core" CDS ids as the longest CDS representing a given gene
        core_cds_ids = self.calculate_core_cds_ids(self.diamond_results,
                                                   accession_numbers_gene_mapping, cds_id_accession_numbers_mapping)

        # Collate metadata of interest calculated so far, in preparation for adding to CDS metadata file
        diamond_results_metadata_fields = ["accession_number", "gene", "core_cds"]
        diamond_results_metadata = self.collate_diamond_results_metadata(cds_id_accession_numbers_mapping,
                                                                         accession_numbers_gene_mapping, core_cds_ids)
        
        # Append new metadata to CDS metadata file
        tmpfile_path = self.write_appended_metadata_to_tempfile(self.cds_metadata, self.diamond_results.parent,
                                                                diamond_results_metadata, diamond_results_metadata_fields,
                                                                "cds_id")
        self.copy_file(tmpfile_path, self.cds_metadata)
        tmpfile_path.unlink()

    @staticmethod
    def extract_cds_id_accession_numbers_mapping(diamond_results: Path) -> dict[str]:
        cds_id_accession_numbers_mapping = dict()
        with diamond_results.open() as inhandle:
            diamond_reader = reader(inhandle, delimiter="\t")
            for line in diamond_reader:
                cds_id = line[0]
                accession_number = line[1]
                cds_id_accession_numbers_mapping[cds_id] = accession_number
        return cds_id_accession_numbers_mapping

    @staticmethod
    def extract_accession_numbers(cds_id_accession_numbers_mapping: dict[str]) -> set[str]:
        return set(cds_id_accession_numbers_mapping.values())

    @staticmethod
    def write_accession_numbers(accession_numbers_out: Path, accession_numbers: set[str]) -> None:
        with accession_numbers_out.open("w") as outhandle:
            for accession_number in accession_numbers:
                outhandle.write(f"{accession_number}\n")

    @staticmethod
    def submit_accession_numbers_with_ncbi_datasets(accession_numbers_out: Path) -> dict[str]:
        print("Starting NCBI Datasets submission\n(This may take awhile)")
        datasets_command = ["datasets", "summary", "gene", "accession", "--inputfile", f"{accession_numbers_out}",
                            "--report", "product", "--as-json-lines"]
        
        dataformat_command = ["dataformat", "tsv", "gene-product", "--elide-header", "--fields", "transcript-protein-accession,symbol"]

        p1 = subprocess.Popen(datasets_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2 = subprocess.Popen(dataformat_command, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        p1.stdout.close()

        accession_numbers_gene_mapping = {}

        while p2.poll() is None and (line := p2.stdout.readline()) != "":
            gene_info = line.strip().split()
            if len(gene_info) == 2:
                accession_numbers_gene_mapping[gene_info[0]] = gene_info[1]
            else:
                pass
        p2.wait()

        if p2.poll() != 0:
            print(p2.stderr.readlines())
            raise Exception("Datasets submission did not complete successfully")
        print("NCBI Submission complete")
        return accession_numbers_gene_mapping

    @staticmethod
    def remove_extraneous_accession_numbers(accession_numbers_gene_mapping: dict[str], accession_numbers: set[str]):
        return {k: v for k, v in accession_numbers_gene_mapping.items() if k in accession_numbers}

    @staticmethod
    def upper_case_genes(accession_numbers_gene_mapping: dict[str]):
        return {k: v.upper() for k, v in accession_numbers_gene_mapping.items()}

    @staticmethod
    def extract_cds_id_len_mapping(diamond_results: Path) -> dict[str]:
        cds_id_len_mapping = dict()
        with diamond_results.open() as inhandle:
            diamond_reader = reader(inhandle, delimiter="\t")
            for line in diamond_reader:
                cds_id = line[0]
                cds_len = int(line[-1])
                cds_id_len_mapping[cds_id] = cds_len
        return cds_id_len_mapping

    @classmethod
    def calculate_core_cds_ids(cls, diamond_results: Path, accession_numbers_gene_mapping: dict[str],
                               cds_id_accession_numbers_mapping: dict[str]) -> set[str]:
        cds_id_len_mapping = cls.extract_cds_id_len_mapping(diamond_results)

        longest_cds = {gene: {"cds_id": "", "cds_len": 0} for gene in set(accession_numbers_gene_mapping.values())}
        for cds_id, cds_len in cds_id_len_mapping.items():
            accession_number = cds_id_accession_numbers_mapping[cds_id]
            try:
                gene = accession_numbers_gene_mapping[accession_number]
            except KeyError:
                continue
            if cds_len > longest_cds[gene]["cds_len"]:
                longest_cds[gene]["cds_id"] = cds_id
                longest_cds[gene]["cds_len"] = cds_len

        core_cds_ids = {cds_info["cds_id"] for cds_info in longest_cds.values()}
        return core_cds_ids

    @staticmethod
    def collate_diamond_results_metadata(cds_id_accession_numbers_mapping: dict[str],
                                       accession_numbers_gene_mapping: dict[str],
                                       core_cds_ids: set[str]) -> dict[Any]:
        diamond_results_metadata = {}
        for cds_id, accession_number in cds_id_accession_numbers_mapping.items():
            try:
                gene = accession_numbers_gene_mapping[accession_number]
            except KeyError:
                gene = ""
            core_cds = True if cds_id in core_cds_ids else ""
            metadata = {"accession_number": accession_number, "gene": gene, "core_cds": core_cds}
            diamond_results_metadata[cds_id] = metadata
        return diamond_results_metadata

    @staticmethod
    def write_appended_metadata_to_tempfile(metadata_path: Path, outdir: Path,
                                            new_metadata: dict[Any], new_fields: list[str],
                                            primary_key_column: str) -> Path:
        print("Starting to write to tmpfile\n(This may take a moment)")

        empty_metadata_singlet = {field: "" for field in new_fields}

        with metadata_path.open() as inhandle, NamedTemporaryFile(dir=outdir, mode="w", delete=False) as tmpfile:
            reader = DictReader(inhandle)
            write_field_names = reader.fieldnames + [field for field in new_fields if field not in reader.fieldnames]

            writer = DictWriter(tmpfile, fieldnames=write_field_names)
            writer.writeheader()
            for data in reader:
                primary_key = data[primary_key_column]
                try:
                    new_metadata_singlet = new_metadata[primary_key]
                    data.update(new_metadata_singlet)
                except KeyError:
                    data.update(empty_metadata_singlet)
                writer.writerow(data)
        return Path(tmpfile.name)

    @staticmethod
    def copy_file(infile_path: Path, outfile_path: Path) -> None:
        print("Starting to copy tmpfile to metadata file")
        with infile_path.open("rb") as inhandle, outfile_path.open("wb") as outhandle:
            copyfileobj(inhandle, outhandle)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-diamond_results", type=str, required=True)
    parser.add_argument("-accesion_numbers_out", type=str, required=True)
    parser.add_argument("-cds_metadata", type=str, required=True)
    parser.add_argument("-upper", type=bool, required=False, default=True)
    args = parser.parse_args()

    ga = GeneAssigner(Path(args.diamond_results), Path(args.accesion_numbers_out), Path(args.cds_metadata))
    ga.run(args.upper)