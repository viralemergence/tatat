from argparse import ArgumentParser
import pandas as pd # type: ignore
from pathlib import Path
import sqlite3

class SqliteDbManager:
    def __init__(self, sqlite_db_dir: Path) -> None:
        self.sqlite_db = sqlite_db_dir / "tatat.db"

    def create_and_insert_samples_table(self, sample_metadata: Path) -> None:
        sample_metadata = self.extract_sample_metadata(sample_metadata)
        with sqlite3.connect(self.sqlite_db) as connection:
            self.create_samples_table(connection, sample_metadata)
            sample_metadata.to_sql("samples", connection, if_exists="append", index=False)

    @staticmethod
    def extract_sample_metadata(sample_metadata_path: Path) -> pd.DataFrame:
        return pd.read_csv(sample_metadata_path)

    @staticmethod
    def create_samples_table(connection: sqlite3.Connection, sample_metadata: pd.DataFrame) -> None:
        cursor = connection.cursor()
        cursor.execute("DROP TABLE IF EXISTS samples")
        cursor.execute('''CREATE TABLE samples
                        (uid TEXT NOT NULL PRIMARY KEY,
                        transcriptome TEXT NOT NULL)''')

        for column in sample_metadata.columns:
            if column not in ["uid", "transcriptome"]:
                cursor.execute(f"ALTER TABLE samples ADD COLUMN {column} TEXT")
        connection.commit()

    def create_transcripts_table(self) -> None:
        with sqlite3.connect(self.sqlite_db) as connection:
            cursor = connection.cursor()
            cursor.execute("DROP TABLE IF EXISTS transcripts")
            cursor.execute('''CREATE TABLE transcripts
                           (uid INTEGER NOT NULL PRIMARY KEY,
                           sample_uid TEXT NOT NULL,
                           length INTEGER NOT NULL,
                           transcript_class TEXT,
                           evigene_pass INTEGER,
                           cds_ids TEXT)''')
            connection.commit()

    def create_cds_table(self) -> None:
        with sqlite3.connect(self.sqlite_db) as connection:
            cursor = connection.cursor()
            cursor.execute("DROP TABLE IF EXISTS cds")
            cursor.execute('''CREATE TABLE cds
                           (uid INTEGER NOT NULL PRIMARY KEY,
                           transcript_uid INTEGER NOT NULL,
                           evigene_class TEXT NOT NULL,
                           strand TEXT NOT NULL,
                           start INTEGER NOT NULL,
                           end INTEGER NOT NULL,
                           length INTEGER NOT NULL,
                           accession_number TEXT,
                           gene_symbol TEXT,
                           unambiguous_gene TEXT,
                           core_cds INTEGER)''')
            connection.commit()

    def create_accession_numbers_table(self) -> None:
        with sqlite3.connect(self.sqlite_db) as connection:
            cursor = connection.cursor()
            cursor.execute("DROP TABLE IF EXISTS accession_numbers")
            cursor.execute('''CREATE TABLE accession_numbers
                           (accession_number TEXT NOT NULL PRIMARY KEY,
                           gene_symbol TEXT NOT NULL)''')
            connection.commit()

    def create_ncrna_table(self) -> None:
        with sqlite3.connect(self.sqlite_db) as connection:
            cursor = connection.cursor()
            cursor.execute("DROP TABLE IF EXISTS ncrna")
            cursor.execute('''CREATE TABLE ncrna
                           (uid INTEGER NOT NULL PRIMARY KEY,
                           accession_number TEXT,
                           gene_symbol TEXT,
                           core_ncrna INTEGER)''')
            connection.commit()

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-sqlite_db_dir", type=str, required=True)
    parser.add_argument("-sample_metadata", type=str, required=False)
    parser.add_argument("-create_transcripts_table", action="store_true", required=False)
    parser.add_argument("-create_cds_table", action="store_true", required=False)
    parser.add_argument("-create_acc_num_table", action="store_true", required=False)
    parser.add_argument("-create_ncrna_table", action="store_true", required=False)
    args = parser.parse_args()

    sdm = SqliteDbManager(Path(args.sqlite_db_dir))
    if args.sample_metadata:
        sdm.create_and_insert_samples_table(Path(args.sample_metadata))
    if args.create_transcripts_table:
        sdm.create_transcripts_table()
    if args.create_cds_table:
        sdm.create_cds_table()
    if args.create_acc_num_table:
        sdm.create_accession_numbers_table()
    if args.create_ncrna_table:
        sdm.create_ncrna_table()