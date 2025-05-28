from argparse import ArgumentParser
from pathlib import Path
import sqlite3
import subprocess

class SraReadDownload:
    def __init__(self, accession_number: str, download_dir: Path, collate_dir: Path, sqlite_db: Path, testing: bool) -> None:
        self.accession_number = accession_number
        self.download_dir = download_dir
        self.collate_dir = collate_dir
        self.sqlite_db = sqlite_db
        self.testing = testing

    def run(self) -> None:
        fastq_files = self.download_reads(self.accession_number, self.download_dir, self.testing)
        self.move_read_files(fastq_files, self.collate_dir)
        self.update_samples_table_read_file_names(self.sqlite_db, self.accession_number, fastq_files)

    @staticmethod
    def download_reads(accession_number: str, outdir: Path, testing: bool) -> list[Path]:
        print("Starting prefetch")
        prefetch_command = ["prefetch", accession_number, "-O", outdir]
        result = subprocess.run(prefetch_command, capture_output=True, text=True)

        print("Starting decompress")
        fastq_outdir = outdir / Path(accession_number)
        if testing:
            decompress_command = ["fastq-dump", fastq_outdir, "--outdir", fastq_outdir, "--split-3", "-X", "25000"]
        else:
            decompress_command = ["fasterq-dump", fastq_outdir, "--temp", fastq_outdir, "--outdir", fastq_outdir]
        result = subprocess.run(decompress_command, capture_output=True, text=True)
        
        acceptable_fastq_files = [f"{accession_number}.fastq", f"{accession_number}_1.fastq", f"{accession_number}_2.fastq"]
        fastq_files = [file for file in fastq_outdir.iterdir() if file.name in acceptable_fastq_files]
        return sorted(fastq_files, reverse=False)

    @staticmethod
    def move_read_files(fastq_files: list[Path], outdir: Path) -> None:
        print("Starting fastq move to collated directory")
        for fastq in fastq_files:
            mv_command = ["mv", fastq, outdir]
            result = subprocess.run(mv_command, capture_output=True, text=True)

    @staticmethod
    def update_samples_table_read_file_names(sqlite_db: Path, accession_number: str,
                                             fastq_files: list[Path]) -> None:
        with sqlite3.connect(sqlite_db) as connection:
            cursor = connection.cursor()
            for i, fastq_file in enumerate(fastq_files, 1):
                sql_statement = f"UPDATE samples SET r{i}_reads = ? WHERE uid = ?"
                cursor.execute(sql_statement, (fastq_file.name, accession_number))
                connection.commit()

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-sra_number", type=str, required=True)
    parser.add_argument("-download_dir", type=str, required=True)
    parser.add_argument("-collate_dir", type=str, required=True)
    parser.add_argument("-sqlite_db", type=str, required=True)
    parser.add_argument("-testing", action="store_true", required=False)
    args = parser.parse_args()

    srd = SraReadDownload(args.sra_number, Path(args.download_dir),
                          Path(args.collate_dir), Path(args.sqlite_db),
                          args.testing)
    srd.run()