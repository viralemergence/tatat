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
        # Extract ncrna id to accession number mapping, then filter by transcriptome
        ncrna_id_accession_numbers_mapping = self.extract_ncrna_id_accession_numbers_mapping(self.blast_results)
        ncrna_id_accession_numbers_mapping = (
            self.filter_ncrna_id_accession_numbers_mapping_by_transcriptome(ncrna_id_accession_numbers_mapping,
                                                                          self.sqlite_db, transcriptome)
                                            )

        # Extract accession numbers to gene symbols mapping
        accession_numbers_gene_mapping = self.extract_accession_number_gene_symbol_mapping(self.sqlite_db)
        print(f"Accession number count: {len(accession_numbers_gene_mapping)}")

        # Assign best "gene" hit to ncrna id
        ncrna_id_best_gene_mapping = self.assign_best_gene_to_ncrna_ids(accession_numbers_gene_mapping, ncrna_id_accession_numbers_mapping)
        print(len(ncrna_id_best_gene_mapping))

        # Calculate "core" ncRNA ids as the longest sequence representing a given "gene"
        core_ncrna_ids = self.calculate_core_ncrna_ids(self.blast_results, ncrna_id_best_gene_mapping)
        print(f"Core ncRNA count: {len(core_ncrna_ids)}")

        # Collate metadata of interest calculated so far, in preparation for adding to ncRNA metadata table
        blast_results_metadata = self.collate_blast_results_metadata(ncrna_id_accession_numbers_mapping,
                                                                         ncrna_id_best_gene_mapping, core_ncrna_ids)
        
        # Insert new metadata to ncRNA metadata table
        with sqlite3.connect(self.sqlite_db) as connection:
            self.update_ncrna_table(connection, blast_results_metadata)

    @staticmethod
    def extract_ncrna_id_accession_numbers_mapping(blast_results: Path) -> dict[list[str]]:
        ncrna_id_accession_numbers_mapping = defaultdict(list)
        with blast_results.open() as inhandle:
            blast_reader = reader(inhandle, delimiter="\t")
            for line in blast_reader:
                ncrna_id = int(line[0])
                accession_number = line[1]
                ncrna_id_accession_numbers_mapping[ncrna_id].append(accession_number)
        return dict(ncrna_id_accession_numbers_mapping)

    @classmethod
    def filter_ncrna_id_accession_numbers_mapping_by_transcriptome(cls, ncrna_id_accession_numbers_mapping: dict[list[str]],
                                                                 sqlite_db: Path, transcriptome: str) -> dict[list[str]]:
        transcriptome_filtered_ncrna_ids = cls.extract_transcriptome_filtered_ncrna_ids(sqlite_db, transcriptome)
        return {k: v for k, v in ncrna_id_accession_numbers_mapping.items() if k in transcriptome_filtered_ncrna_ids}

    @staticmethod
    def extract_transcriptome_filtered_ncrna_ids(sqlite_db: Path, transcriptome: str) -> set[int]:
        print(f"\nExtracting ncrna ids that belong to transcriptome: {transcriptome}")
        with sqlite3.connect(sqlite_db) as connection:
            cursor = connection.cursor()
            sql_query = ("SELECT n.uid "
                         "FROM ncrna n "
                         "LEFT OUTER JOIN transcripts t ON n.uid = t.uid "
                         "LEFT OUTER JOIN samples s ON t.sample_uid = s.uid "
                         f"WHERE s.transcriptome = '{transcriptome}'")
            cursor.execute(sql_query)
            return {row[0] for row in cursor.fetchall()}

    @staticmethod
    def extract_accession_number_gene_symbol_mapping(sqlite_db: Path) -> dict[str]:
        with sqlite3.connect(sqlite_db) as connection:
            cursor = connection.cursor()
            sql_query = "SELECT accession_number,gene_symbol FROM nc_accession_numbers"
            cursor.execute(sql_query)
            return {row[0]: row[1] for row in cursor.fetchall()}

    @staticmethod
    def assign_best_gene_to_ncrna_ids(accession_numbers_gene_mapping: dict[str],
                               ncrna_id_accession_numbers_mapping: dict[list[str]]) -> dict[dict[str]]:

        ncrna_id_best_gene_mapping = {}

        for ncrna_id, accession_numbers in ncrna_id_accession_numbers_mapping.items():
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
            
            ncrna_id_best_gene_mapping[ncrna_id] = best_hit

        return ncrna_id_best_gene_mapping

    @classmethod
    def calculate_core_ncrna_ids(cls, blast_results: Path, ncrna_id_best_gene_mapping: dict[dict[str]]) -> set[int]:
        ncrna_id_len_mapping = cls.extract_ncrna_id_len_mapping(blast_results)

        longest_ncrna = {gene: {"ncrna_id": "", "ncrna_len": 0} for gene in set(ncrna_info["gene"] for ncrna_info in ncrna_id_best_gene_mapping.values())}
        for ncrna_id, ncrna_len in ncrna_id_len_mapping.items():
            try:
                gene = ncrna_id_best_gene_mapping[ncrna_id]["gene"]
            except KeyError:
                continue

            if ncrna_len > longest_ncrna[gene]["ncrna_len"]:
                longest_ncrna[gene]["ncrna_id"] = ncrna_id
                longest_ncrna[gene]["ncrna_len"] = ncrna_len

        core_ncrna_ids = {ncrna_info["ncrna_id"] for ncrna_info in longest_ncrna.values()}
        return core_ncrna_ids

    @staticmethod
    def extract_ncrna_id_len_mapping(blast_results: Path) -> dict[int]:
        ncrna_id_len_mapping = dict()
        with blast_results.open() as inhandle:
            blast_reader = reader(inhandle, delimiter="\t")
            for line in blast_reader:
                ncrna_id = int(line[0])
                ncrna_len = int(line[-1])
                ncrna_id_len_mapping[ncrna_id] = ncrna_len
        return ncrna_id_len_mapping

    @staticmethod
    def collate_blast_results_metadata(ncrna_id_accession_numbers_mapping: dict[str],
                                       ncrna_id_best_gene_mapping: dict[dict[str]],
                                       core_ncrna_ids: set[int]) -> dict[Any]:
        blast_results_metadata = {}
        for ncrna_id, accession_numbers in ncrna_id_accession_numbers_mapping.items():
            try:
                accession_number = ncrna_id_best_gene_mapping[ncrna_id]["accession"]
                gene = ncrna_id_best_gene_mapping[ncrna_id]["gene"]
                if gene.startswith("LOC") or gene.startswith("CUN"):
                    unambiguous = 0
                else:
                    unambiguous = 1
            except KeyError:
                accession_number = accession_numbers[0]
                gene = None
                unambiguous = None
            if gene is None:
                core_ncrna = None
            else:
                core_ncrna = 1 if ncrna_id in core_ncrna_ids else 0

            metadata = {"accession_number": accession_number, "gene_symbol": gene, "unambiguous_gene": unambiguous, "core_ncrna": core_ncrna}
            blast_results_metadata[ncrna_id] = metadata
        return blast_results_metadata

    @staticmethod
    def update_ncrna_table(connection: sqlite3.Connection, blast_results_metadata: dict[Any]) -> None:
        cursor = connection.cursor()
        values = [(data["accession_number"], data["gene_symbol"],
                    data["core_ncrna"], uid)
                    for uid, data in blast_results_metadata.items()]
        
        sql_statement = ("UPDATE ncrna "
                            "SET accession_number = ?, gene_symbol = ?, "
                            "core_ncrna = ? "
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