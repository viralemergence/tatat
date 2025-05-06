from argparse import ArgumentParser
from pathlib import Path
import sqlite3
from typing import Any, Iterator, TextIO

class DeNovoAssemblyManager:
    def __init__(self, assembly_fasta_dir: Path, merged_path: Path, sqlite_db: Path) -> None:
        self.assembly_fasta_dir = assembly_fasta_dir
        self.merged_path = merged_path
        self.sqlite_db = sqlite_db

    def run(self) -> None:
        transcript_metadata = self.merge_fastas_and_extract_metadata(self.assembly_fasta_dir, self.merged_path)

        with sqlite3.connect(self.sqlite_db) as connection:
            self.create_transcripts_table(connection)
            self.insert_into_transcripts_table(connection, transcript_metadata)

    @classmethod
    def merge_fastas_and_extract_metadata(cls, assembly_fasta_dir: Path, merged_path: Path) -> list[tuple[Any]]:
        fasta_files = cls.get_file_list(assembly_fasta_dir)
        print(f"A total of {len(fasta_files)} files detected")

        transcript_metadata = list()
        with merged_path.open("w") as merged_outhandle:
            transcript_id = 0
            for fasta_file in fasta_files:
                sample_uid = fasta_file.stem
                print(f"Starting on sample: {sample_uid}")
                for fasta_seq in cls.fasta_chunker(fasta_file):
                    transcript_id += 1
                    cls.write_renamed_fasta_seq(transcript_id, merged_outhandle, fasta_seq[1:])

                    transcript_len = len("".join(fasta_seq[1:]))
                    transcript_metadata.append((transcript_id, sample_uid, transcript_len))
        return transcript_metadata

    @staticmethod
    def get_file_list(directory: Path) -> list[Path]:
        return sorted([file for file in directory.iterdir()])

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
    def write_renamed_fasta_seq(seq_id: int, outhandle: TextIO, fasta_seq: list[str]) -> None:
        # NOTE: fasta_seq should JUST be sequence, not including header
        # Otherwise, new and old header will both be written
        new_header = f">{seq_id}"
        outhandle.write(f"{new_header}\n")
        for line in fasta_seq:
            outhandle.write(f"{line}\n")

    @staticmethod
    def create_transcripts_table(connection: sqlite3.Connection) -> None:
            cursor = connection.cursor()
            cursor.execute("DROP TABLE IF EXISTS transcripts")
            cursor.execute('''CREATE TABLE transcripts
                           (uid INTEGER NOT NULL PRIMARY KEY,
                           sample_uid TEXT NOT NULL,
                           length INTEGER NOT NULL
                           )''')
            connection.commit()

    @staticmethod
    def insert_into_transcripts_table(connection: sqlite3.Connection, transcript_metadata: list[tuple[Any]]) -> None:
        cursor = connection.cursor()
        sql_statement = '''INSERT INTO transcripts VALUES (?,?,?)'''
        cursor.executemany(sql_statement, transcript_metadata)
        connection.commit()

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-assembly_fasta_dir", type=str, required=True)
    parser.add_argument("-merged_path", type=str, required=True)
    parser.add_argument("-sqlite_db", type=str, required=True)
    args = parser.parse_args()

    dnam = DeNovoAssemblyManager(Path(args.assembly_fasta_dir), Path(args.merged_path), Path(args.sqlite_db))
    dnam.run()