from argparse import ArgumentParser
from pathlib import Path
import sqlite3
from typing import Iterator, Union

class NcrnaFastaManager:
    def __init__(self, assembly_fasta: Path, sqlite_db: Path, ncrna_fasta: Path) -> None:
        self.assembly_fasta = assembly_fasta
        self.sqlite_db = sqlite_db
        self.ncrna_fasta = ncrna_fasta

    def run(self, transcriptome: str, add_gene_name: True=None) -> None:
        core_ncrna_ids = self.extract_transcriptome_filtered_core_ncrna_ids(self.sqlite_db, transcriptome)
        if add_gene_name:
            ncrna_id_gene_mapping = self.extract_ncrna_id_gene_mapping(self.sqlite_db)
        else:
            ncrna_id_gene_mapping = None

        self.extract_and_write_ncrna(self.assembly_fasta, self.ncrna_fasta, core_ncrna_ids, ncrna_id_gene_mapping)

    @staticmethod
    def extract_transcriptome_filtered_core_ncrna_ids(sqlite_db: Path, transcriptome: str) -> set[int]:
        print(f"\nExtracting ncRNA ids that belong to transcriptome: {transcriptome}")
        with sqlite3.connect(sqlite_db) as connection:
            cursor = connection.cursor()
            sql_query = ("SELECT n.uid "
                         "FROM ncrna n "
                         "LEFT OUTER JOIN transcripts t ON n.uid = t.uid "
                         "LEFT OUTER JOIN samples s ON t.sample_uid = s.uid "
                         f"WHERE s.transcriptome = '{transcriptome}' "
                         "AND n.core_ncrna = 1")
            cursor.execute(sql_query)
            return {row[0] for row in cursor.fetchall()}

    @staticmethod
    def extract_ncrna_id_gene_mapping(sqlite_db: Path) -> dict[str]:
        print("Extracting ncRNA id gene mapping")
        with sqlite3.connect(sqlite_db) as connection:
            cursor = connection.cursor()
            sql_query = ("SELECT uid, gene_symbol "
                         "FROM ncrna")
            cursor.execute(sql_query)
            return {row[0]: row[1] for row in cursor.fetchall()}

    @classmethod
    def extract_and_write_ncrna(cls, assembly_fasta: Path, ncrna_fasta: Path,
                                core_ncrna_ids: set[int], ncrna_gene_mapping: Union[None, dict[str]]) -> None:
        print("Starting ncRNA extraction and writing\n(This may take awhile)")
        with ncrna_fasta.open("w") as ncrna_outhandle:
            for fasta_seq in cls.fasta_chunker(assembly_fasta):
                transcript_id = int(fasta_seq[0][1:])
                if transcript_id not in core_ncrna_ids:
                    continue

                transcript_seq = "".join(fasta_seq[1:])
                if ncrna_gene_mapping:
                    gene = ncrna_gene_mapping[transcript_id]
                    ncrna_outhandle.write(f">{transcript_id};{gene}\n")
                else:
                    ncrna_outhandle.write(f">{transcript_id}\n")
                ncrna_outhandle.write(f"{transcript_seq}\n")

    @staticmethod
    def fasta_chunker(fasta_path: Path) -> Iterator[list[str]]:
        fasta_seq = []
        first_chunk = True
        with fasta_path.open() as inhandle:
            for line in inhandle:
                line = line.strip()
                if not line.startswith(">"):
                    fasta_seq.append(line)
                else:
                    if first_chunk:
                        fasta_seq.append(line)
                        first_chunk = False
                        continue
                    yield fasta_seq
                    fasta_seq = [line]
            if fasta_seq:
                yield fasta_seq

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-assembly_fasta", type=str, required=True)
    parser.add_argument("-sqlite_db", type=str, required=True)
    parser.add_argument("-ncrna_fasta", type=str, required=True)
    parser.add_argument("-add_gene_name", action="store_true", required=False)
    parser.add_argument("-transcriptome", type=str, required=True)
    args = parser.parse_args()

    nfm = NcrnaFastaManager(Path(args.assembly_fasta), Path(args.sqlite_db), Path(args.ncrna_fasta))
    nfm.run(args.transcriptome, args.add_gene_name)
    print("\nFinished")