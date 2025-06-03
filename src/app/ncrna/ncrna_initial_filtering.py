from argparse import ArgumentParser
from hashlib import sha256
from pathlib import Path
import sqlite3
from typing import Iterator

class NcrnaInitialManager:
    def __init__(self, sqlite_db: Path, transcripts_fasta: Path) -> None:
        self.sqlite_db = sqlite_db
        self.transcripts_fasta = transcripts_fasta

    def run(self, transcriptome: str) -> None:
        transcript_ids = self.extract_length_and_gene_filtered_transcript_ids(self.sqlite_db, transcriptome)
        transcript_ids = self.remove_ids_of_duplicate_sequences(self.transcripts_fasta, transcript_ids)
        print(len(transcript_ids))

        values = [(id,) for id in transcript_ids]
        with sqlite3.connect(self.sqlite_db, timeout=600) as connection:
            cursor = connection.cursor()
            sql_statement = ("INSERT INTO ncrna (uid) VALUES (?)")
            cursor.executemany(sql_statement, values)
            connection.commit()

    @staticmethod
    def extract_length_and_gene_filtered_transcript_ids(sqlite_db: Path, transcriptome: str) -> set[int]:
        print("\nExtracting filtered transcript ids")
        with sqlite3.connect(sqlite_db, timeout=600) as connection:
            cursor = connection.cursor()
            sql_query = ("SELECT t.uid "
                         "FROM transcripts t "
                         "LEFT OUTER JOIN cds c ON t.uid = c.transcript_uid "
                         "LEFT OUTER JOIN samples s ON t.sample_uid = s.uid "
                         "WHERE t.length < 5000 AND c.gene_symbol IS NULL "
                         f"AND s.transcriptome = '{transcriptome}'")
            cursor.execute(sql_query)
            return {row[0] for row in cursor.fetchall()}

    @classmethod
    def remove_ids_of_duplicate_sequences(cls, transcripts_fasta: Path, transcript_ids: set[int]) -> set[int]:
        print("Removing transcript ids of duplicate sequences")
        non_duplicate_ids = set()
        sha256_hashes = set()
        for fasta_seq in cls.fasta_chunker(transcripts_fasta):
            transcript_id = int(fasta_seq[0][1:])

            if transcript_id not in transcript_ids:
                continue

            seq = "".join(fasta_seq[1:])
            seq_hash = sha256(seq.encode("utf-8")).digest()
            if seq_hash in sha256_hashes:
                continue

            rev_seq = cls.reverse_translate_dna(seq)
            rev_seq_hash = sha256(rev_seq.encode("utf-8")).digest()
            if rev_seq_hash in sha256_hashes:
                continue
            
            sha256_hashes.add(seq_hash)
            sha256_hashes.add(rev_seq_hash)
            non_duplicate_ids.add(transcript_id)
        return non_duplicate_ids

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

    @staticmethod
    def reverse_translate_dna(dna_sequence: str) -> str:
        translation_mapping = str.maketrans("ATCG", "TAGC")
        return dna_sequence.translate(translation_mapping)[::-1]

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-sqlite_db", type=str, required=True)
    parser.add_argument("-transcripts_fasta", type=str, required=True)
    parser.add_argument("-transcriptome", type=str, required=True)
    args = parser.parse_args()

    nim = NcrnaInitialManager(Path(args.sqlite_db), Path(args.transcripts_fasta))
    nim.run(args.transcriptome)