from argparse import ArgumentParser
from collections import defaultdict
from csv import reader
from pathlib import Path
import sqlite3
from typing import Any

class GeneAssigner:
    def __init__(self, blast_results: Path, sqlite_db: Path) -> None:
        self.blast_results = blast_results
        self.sqlite_db = sqlite_db

    def run(self, transcriptome: str) -> None:
        # Extract cds id to accession number mapping, then filter by transcriptome
        cds_id_accession_numbers_mapping = self.extract_cds_id_accession_numbers_mapping(self.blast_results)
        cds_id_accession_numbers_mapping = (
            self.filter_cds_id_accession_numbers_mapping_by_transcriptome(cds_id_accession_numbers_mapping,
                                                                          self.sqlite_db, transcriptome)
                                            )

        # Extract accession numbers to gene symbols mapping
        accession_numbers_gene_mapping = self.extract_accession_number_gene_symbol_mapping(self.sqlite_db)
        print(f"Accession number count: {len(accession_numbers_gene_mapping)}")

        # Assign best "gene" hit to cds id
        cds_id_best_gene_mapping = self.assign_best_gene_to_cds_ids(accession_numbers_gene_mapping, cds_id_accession_numbers_mapping)
        print(len(cds_id_best_gene_mapping))

        # Calculate "core" CDS ids as the longest CDS representing a given gene
        core_cds_ids = self.calculate_core_cds_ids(self.blast_results, cds_id_best_gene_mapping)
        print(f"Core CDS count: {len(core_cds_ids)}")

        # Collate metadata of interest calculated so far, in preparation for adding to CDS metadata file
        blast_results_metadata = self.collate_blast_results_metadata(cds_id_accession_numbers_mapping,
                                                                         cds_id_best_gene_mapping, core_cds_ids)
        
        # Insert new metadata to CDS metadata table
        with sqlite3.connect(self.sqlite_db) as connection:
            self.update_cds_table(connection, blast_results_metadata)

    @staticmethod
    def extract_cds_id_accession_numbers_mapping(blast_results: Path) -> dict[list[str]]:
        cds_id_accession_numbers_mapping = defaultdict(list)
        with blast_results.open() as inhandle:
            blast_reader = reader(inhandle, delimiter="\t")
            for line in blast_reader:
                cds_id = int(line[0])
                accession_number = line[1]
                cds_id_accession_numbers_mapping[cds_id].append(accession_number)
        return dict(cds_id_accession_numbers_mapping)

    @classmethod
    def filter_cds_id_accession_numbers_mapping_by_transcriptome(cls, cds_id_accession_numbers_mapping: dict[list[str]],
                                                                 sqlite_db: Path, transcriptome: str) -> dict[list[str]]:
        transcriptome_filtered_cds_ids = cls.extract_transcriptome_filtered_cds_ids(sqlite_db, transcriptome)
        return {k: v for k, v in cds_id_accession_numbers_mapping.items() if k in transcriptome_filtered_cds_ids}

    @staticmethod
    def extract_transcriptome_filtered_cds_ids(sqlite_db: Path, transcriptome: str) -> set[int]:
        print(f"\nExtracting cds ids that belong to transcriptome: {transcriptome}")
        with sqlite3.connect(sqlite_db) as connection:
            cursor = connection.cursor()
            sql_query = ("SELECT c.uid "
                         "FROM cds c "
                         "LEFT OUTER JOIN transcripts t ON c.transcript_uid = t.uid "
                         "LEFT OUTER JOIN samples s ON t.sample_uid = s.uid "
                         f"WHERE s.transcriptome = '{transcriptome}'")
            cursor.execute(sql_query)
            return {row[0] for row in cursor.fetchall()}

    @staticmethod
    def extract_accession_number_gene_symbol_mapping(sqlite_db: Path) -> dict[str]:
        with sqlite3.connect(sqlite_db) as connection:
            cursor = connection.cursor()
            sql_query = "SELECT accession_number,gene_symbol FROM accession_numbers"
            cursor.execute(sql_query)
            return {row[0]: row[1] for row in cursor.fetchall()}

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
    def calculate_core_cds_ids(cls, blast_results: Path, cds_id_best_gene_mapping: dict[dict[str]]) -> set[int]:
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
    def extract_cds_id_len_mapping(blast_results: Path) -> dict[int]:
        cds_id_len_mapping = dict()
        with blast_results.open() as inhandle:
            blast_reader = reader(inhandle, delimiter="\t")
            for line in blast_reader:
                cds_id = int(line[0])
                cds_len = int(line[-1])
                cds_id_len_mapping[cds_id] = cds_len
        return cds_id_len_mapping

    @staticmethod
    def collate_blast_results_metadata(cds_id_accession_numbers_mapping: dict[str],
                                       cds_id_best_gene_mapping: dict[dict[str]],
                                       core_cds_ids: set[int]) -> dict[Any]:
        blast_results_metadata = {}
        for cds_id, accession_numbers in cds_id_accession_numbers_mapping.items():
            try:
                accession_number = cds_id_best_gene_mapping[cds_id]["accession"]
                gene = cds_id_best_gene_mapping[cds_id]["gene"]
                if gene.startswith("LOC") or gene.startswith("CUN"):
                    unambiguous = 0
                else:
                    unambiguous = 1
            except KeyError:
                accession_number = accession_numbers[0]
                gene = None
                unambiguous = None
            if gene is None:
                core_cds = None
            else:
                core_cds = 1 if cds_id in core_cds_ids else 0

            metadata = {"accession_number": accession_number, "gene_symbol": gene, "unambiguous_gene": unambiguous, "core_cds": core_cds}
            blast_results_metadata[cds_id] = metadata
        return blast_results_metadata

    @staticmethod
    def update_cds_table(connection: sqlite3.Connection, blast_results_metadata: dict[Any]) -> None:
        cursor = connection.cursor()
        values = [(data["accession_number"], data["gene_symbol"],
                    data["unambiguous_gene"], data["core_cds"], uid)
                    for uid, data in blast_results_metadata.items()]
        
        sql_statement = ("UPDATE cds "
                            "SET accession_number = ?, gene_symbol = ?, "
                            "unambiguous_gene = ?, core_cds = ? "
                            "WHERE uid = ?")
        cursor.executemany(sql_statement, values)
        connection.commit()

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-blast_results", type=str, required=True)
    parser.add_argument("-sqlite_db", type=str, required=True)
    parser.add_argument("-transcriptome", type=str, required=True)
    args = parser.parse_args()

    ga = GeneAssigner(Path(args.blast_results), Path(args.sqlite_db))
    ga.run(args.transcriptome)