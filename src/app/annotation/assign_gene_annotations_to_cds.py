from argparse import ArgumentParser
from collections import defaultdict
from csv import DictReader, DictWriter, reader
from pathlib import Path
from shutil import copyfileobj
from tempfile import NamedTemporaryFile
from typing import Any

class GeneAssigner:
    def __init__(self, blast_results: Path, accession_gene_mapping: Path, cds_metadata: Path) -> None:
        self.blast_results = blast_results
        self.accession_gene_mapping_path = accession_gene_mapping
        self.cds_metadata = cds_metadata

    def run(self) -> None:
        # Extract cds id to accession number mapping
        cds_id_accession_numbers_mapping = self.extract_cds_id_accession_numbers_mapping(self.blast_results)

        # Extract accession numbers to gene symbols mapping
        accession_numbers_gene_mapping = self.extract_accession_number_gene_symbol_mapping(self.accession_gene_mapping_path)
        print(len(accession_numbers_gene_mapping))

        # Assign best "gene" hit to cds id
        cds_id_best_gene_mapping = self.assign_best_gene_to_cds_ids(accession_numbers_gene_mapping, cds_id_accession_numbers_mapping)
        print(len(cds_id_best_gene_mapping))

        # Calculate "core" CDS ids as the longest CDS representing a given gene
        core_cds_ids = self.calculate_core_cds_ids(self.blast_results, cds_id_best_gene_mapping)
        print(f"Core CDS count: {len(core_cds_ids)}")

        # Collate metadata of interest calculated so far, in preparation for adding to CDS metadata file
        blast_results_metadata_fields = ["accession_number", "gene", "unambiguous_gene", "core_cds"]
        blast_results_metadata = self.collate_blast_results_metadata(cds_id_accession_numbers_mapping,
                                                                         cds_id_best_gene_mapping, core_cds_ids)
        
        # Append new metadata to CDS metadata file
        tmpfile_path = self.write_appended_metadata_to_tempfile(self.cds_metadata, self.blast_results.parent,
                                                                blast_results_metadata, blast_results_metadata_fields,
                                                                "cds_id")
        self.copy_file(tmpfile_path, self.cds_metadata)
        tmpfile_path.unlink()

    @staticmethod
    def extract_cds_id_accession_numbers_mapping(blast_results: Path) -> dict[list[str]]:
        cds_id_accession_numbers_mapping = defaultdict(list)
        with blast_results.open() as inhandle:
            blast_reader = reader(inhandle, delimiter="\t")
            for line in blast_reader:
                cds_id = line[0]
                accession_number = line[1]
                cds_id_accession_numbers_mapping[cds_id].append(accession_number)
        return dict(cds_id_accession_numbers_mapping)

    @staticmethod
    def extract_accession_number_gene_symbol_mapping(accession_gene_mapping_path: Path) -> dict[str]:
        accession_numbers_gene_mapping = {}
        with accession_gene_mapping_path.open() as inhandle:
            for line in inhandle:
                line = line.strip().split(",")
                accession_numbers_gene_mapping[line[0]] = line[1]
        return accession_numbers_gene_mapping

    @staticmethod
    def assign_best_gene_to_cds_ids(accession_numbers_gene_mapping: dict[str],
                               cds_id_accession_numbers_mapping: dict[list[str]]) -> dict[dict[str]]:

        cds_id_best_gene_mapping = {}

        for cds_id, accession_numbers in cds_id_accession_numbers_mapping.items():
            first_hit = None
            gene_hit = None
            for accession_number in accession_numbers:
                try:
                    gene = accession_numbers_gene_mapping[accession_number]
                except KeyError:
                    continue
                
                if first_hit is None:
                    first_hit = {"accession": accession_number, "gene": gene}

                if gene.startswith("LOC"):
                    continue
                else:
                    gene_hit = {"accession": accession_number, "gene": gene}
                    break
            
            if first_hit is None and gene_hit is None:
                continue
            elif first_hit and gene_hit is None:
                best_hit = first_hit
            elif first_hit and gene_hit:
                best_hit = gene_hit
            
            cds_id_best_gene_mapping[cds_id] = best_hit

        return cds_id_best_gene_mapping

    @classmethod
    def calculate_core_cds_ids(cls, blast_results: Path, cds_id_best_gene_mapping: dict[dict[str]]) -> set[str]:
        cds_id_len_mapping = cls.extract_cds_id_len_mapping(blast_results)

        longest_cds = {gene: {"cds_id": "", "cds_len": 0} for gene in set(cds_info["gene"] for cds_info in cds_id_best_gene_mapping.values())}
        for cds_id, cds_len in cds_id_len_mapping.items():
            try:
                gene = cds_id_best_gene_mapping[cds_id]["gene"]
            except KeyError:
                continue

            if cds_len > longest_cds[gene]["cds_len"]:
                longest_cds[gene]["cds_id"] = cds_id
                longest_cds[gene]["cds_len"] = cds_len

        core_cds_ids = {cds_info["cds_id"] for cds_info in longest_cds.values()}
        return core_cds_ids

    @staticmethod
    def extract_cds_id_len_mapping(blast_results: Path) -> dict[str]:
        cds_id_len_mapping = dict()
        with blast_results.open() as inhandle:
            blast_reader = reader(inhandle, delimiter="\t")
            for line in blast_reader:
                cds_id = line[0]
                cds_len = int(line[-1])
                cds_id_len_mapping[cds_id] = cds_len
        return cds_id_len_mapping

    @staticmethod
    def collate_blast_results_metadata(cds_id_accession_numbers_mapping: dict[str],
                                       cds_id_best_gene_mapping: dict[dict[str]],
                                       core_cds_ids: set[str]) -> dict[Any]:
        blast_results_metadata = {}
        for cds_id, accession_numbers in cds_id_accession_numbers_mapping.items():
            try:
                accession_number = cds_id_best_gene_mapping[cds_id]["accession"]
                gene = cds_id_best_gene_mapping[cds_id]["gene"]
                if gene.startswith("LOC") or gene.startswith("CUN"):
                    unambiguous = False
                else:
                    unambiguous = True
            except KeyError:
                accession_number = accession_numbers[0]
                gene = ""
                unambiguous = ""
            core_cds = True if cds_id in core_cds_ids else ""

            metadata = {"accession_number": accession_number, "gene": gene, "unambiguous_gene": unambiguous, "core_cds": core_cds}
            blast_results_metadata[cds_id] = metadata
        return blast_results_metadata

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
    parser.add_argument("-blast_results", type=str, required=True)
    parser.add_argument("-accession_gene_mapping", type=str, required=True)
    parser.add_argument("-cds_metadata", type=str, required=True)
    args = parser.parse_args()

    ga = GeneAssigner(Path(args.blast_results), Path(args.accession_gene_mapping), Path(args.cds_metadata))
    ga.run()