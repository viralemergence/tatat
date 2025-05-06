from argparse import ArgumentParser
from pathlib import Path
import sqlite3
from typing import Any

class SampleMetadataManager:
    def __init__(self, sqlite_db: Path) -> None:
        self.sample_metadata = self.extract_sample_metadata(sqlite_db)

    @staticmethod
    def extract_sample_metadata(sqlite_db: Path) -> list[dict[Any]]:
        with sqlite3.connect(sqlite_db) as connection:
            connection.row_factory = sqlite3.Row
            cursor = connection.cursor()
            cursor.execute("SELECT * from samples ORDER BY uid")
            return [dict(row) for row in cursor.fetchall()]

    def print_uid(self, array_index: int) -> None:
        print(self.sample_metadata[array_index]["uid"], flush=True)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-sqlite_db", type=str, required=True)
    parser.add_argument("-array_index", type=int, required=True)
    args = parser.parse_args()

    smm = SampleMetadataManager(Path(args.sqlite_db))
    smm.print_uid(args.array_index)