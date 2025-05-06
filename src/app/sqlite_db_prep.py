from argparse import ArgumentParser
import pandas as pd # type: ignore
from pathlib import Path
import sqlite3

class SampleMetadataManager:
    def __init__(self, sample_metadata: Path, sqlite_db_dir: Path) -> None:
        self.sample_metadata = self.extract_sample_metadata(sample_metadata)
        self.sqlite_db = sqlite_db_dir / "tatat.db"

    @staticmethod
    def extract_sample_metadata(sample_metadata_path: Path) -> pd.DataFrame:
        return pd.read_csv(sample_metadata_path)

    def run(self) -> None:
        with sqlite3.connect(self.sqlite_db) as connection:
            self.create_samples_table(connection, self.sample_metadata)
            self.sample_metadata.to_sql("samples", connection, if_exists="append", index=False)

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

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-sample_metadata", type=str, required=True)
    parser.add_argument("-sqlite_db_dir", type=str, required=True)
    args = parser.parse_args()

    smm = SampleMetadataManager(Path(args.sample_metadata), Path(args.sqlite_db_dir))
    smm.run()