from argparse import ArgumentParser
from pathlib import Path
import sqlite3
from typing import Iterator, Union

class NcrnaFastaManager:
    def __init__(self, assembly_fasta: Path, sqlite_db: Path, ncrna_fasta: Path) -> None:
        self.assembly_fasta = assembly_fasta
        self.sqlite_db = sqlite_db
        self.ncrna_fasta = ncrna_fasta

    def run(self) -> None:
        ncrna_ids = self.extract_cd_hit_filtered_ncrna_ids(self.sqlite_db)
        self.extract_and_write_ncrna(self.assembly_fasta, self.ncrna_fasta, ncrna_ids)

    @staticmethod
    def extract_cd_hit_filtered_ncrna_ids(sqlite_db: Path) -> set[int]:
        print(f"\nExtracting ncRNA ids for blast")
        with sqlite3.connect(sqlite_db) as connection:
            cursor = connection.cursor()
            sql_query = ("SELECT uid FROM ncrna "
                         "WHERE cd_hit_pass = 1")
            cursor.execute(sql_query)
            return {row[0] for row in cursor.fetchall()}

    @classmethod
    def extract_and_write_ncrna(cls, assembly_fasta: Path, ncrna_fasta: Path,
                                ncrna_ids: set[int]) -> None:
        print("Starting ncRNA extraction and writing\n(This may take awhile)")
        with ncrna_fasta.open("w") as ncrna_outhandle:
            for fasta_seq in cls.fasta_chunker(assembly_fasta):
                transcript_id = int(fasta_seq[0][1:])
                if transcript_id not in ncrna_ids:
                    continue

                for line in fasta_seq:
                    ncrna_outhandle.write(f"{line}\n")

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
    args = parser.parse_args()

    nfm = NcrnaFastaManager(Path(args.assembly_fasta), Path(args.sqlite_db), Path(args.ncrna_fasta))
    nfm.run()
    print("\nFinished")