from argparse import ArgumentParser
from pathlib import Path
import sqlite3

class SampleMetadataManager:
    def __init__(self, sqlite_db: Path) -> None:
        self.sqlite_db = sqlite_db

    def print_uid(self, array_index: int) -> None:
        with sqlite3.connect(self.sqlite_db) as connection:
            cursor = connection.cursor()
            sql_query = "SELECT uid FROM samples ORDER BY uid"
            cursor.execute(sql_query)
            uids = [row[0] for row in cursor.fetchall()]
        print(uids[array_index], flush=True)

    def print_transcriptome(self, array_index: int) -> None:
        with sqlite3.connect(self.sqlite_db) as connection:
            cursor = connection.cursor()
            sql_query = "SELECT DISTINCT transcriptome FROM samples ORDER BY transcriptome"
            cursor.execute(sql_query)
            transcriptomes= [row[0] for row in cursor.fetchall()]
        print(transcriptomes[array_index], flush=True)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-sqlite_db", type=str, required=True)
    parser.add_argument("-array_index", type=int, required=True)
    parser.add_argument("-return_uid", action="store_true", required=False)
    parser.add_argument("-return_transcriptome", action="store_true", required=False)
    args = parser.parse_args()

    smm = SampleMetadataManager(Path(args.sqlite_db))
    if args.return_uid:
        smm.print_uid(args.array_index)
    if args.return_transcriptome:
        smm.print_transcriptome(args.array_index)