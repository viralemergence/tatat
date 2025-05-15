from argparse import ArgumentParser
from pathlib import Path
import sqlite3
import subprocess
from typing import Iterator

class CdHitManager:
    def __init__(self, sqlite_db: Path, transcripts_fasta: Path, ncrna_dir: Path, cds_fasta: Path) -> None:
        self.sqlite_db = sqlite_db
        self.transcripts_fasta = transcripts_fasta
        self.ncrna_dir = ncrna_dir
        self.cds_fasta = cds_fasta

        self.temp_ncrna_fasta = ncrna_dir / "temp_ncrna.fna"
        self.ncrna_cd_hit_est_2d_fasta = ncrna_dir / "ncrna_cd_hit_est_2d.fna"
        self.ncrna_cd_hit_est_fasta = ncrna_dir / "ncrna_cd_hit_est.fna"

    def run(self, cpus: int, memory: int) -> None:
        ncrna_ids = self.extract_ncrna_ids(self.sqlite_db)
        self.write_temporary_ncrna_fasta(ncrna_ids, self.transcripts_fasta, self.temp_ncrna_fasta)
        self.run_cd_hit_est_2d(self.cds_fasta, self.temp_ncrna_fasta, self.ncrna_cd_hit_est_2d_fasta, cpus, memory)
        self.run_cd_hit_est(self.ncrna_cd_hit_est_2d_fasta, self.ncrna_cd_hit_est_fasta, cpus, memory)

    @staticmethod
    def extract_ncrna_ids(sqlite_db: Path) -> set[int]:
        print("\nExtracting ncrna ids")
        with sqlite3.connect(sqlite_db) as connection:
            cursor = connection.cursor()
            sql_query = ("SELECT uid FROM ncrna")
            cursor.execute(sql_query)
            return {row[0] for row in cursor.fetchall()}

    @classmethod
    def write_temporary_ncrna_fasta(cls, ncrna_ids: set[int], transcripts_fasta: Path, temp_ncrna_fasta: Path) -> None:
        print("Writing temporary ncrna fasta\n(This may take a while)")

        with temp_ncrna_fasta.open("w") as outhandle:
            for fasta_seq in cls.fasta_chunker(transcripts_fasta):
                transcript_id = int(fasta_seq[0][1:])
                if transcript_id not in ncrna_ids:
                    continue
                for line in fasta_seq:
                    outhandle.write(f"{line}\n")

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
    def run_cd_hit_est_2d(cds_fasta: Path, temp_ncrna_fasta: Path, ncrna_cd_hit_est_2d_fasta: Path,
                          cpus: int=1, memory: int=1_000) -> None:
        print("\nStarting cd-hit-est-2d\n(This may take even longer)")
        cd_hit_command = ["cd-hit-est-2d",
                          "-i", f"{cds_fasta}",
                          "-i2", f"{temp_ncrna_fasta}",
                          "-o", f"{ncrna_cd_hit_est_2d_fasta}",
                          "-c", "0.99",
                          "-T", f"{cpus}",
                          "-M", f"{memory}",
                          "-s2", "0",
                          "-S2", "999999"]

        p = subprocess.Popen(cd_hit_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        while p.poll() is None and (line := p.stdout.readline()) != "":
            line = line.strip()
            try:
                progress = float(line[:-1])
                if progress % 10 == 0:
                    print(line)
            except ValueError:
                print(line)
        p.wait()
        print(f"Exit code: {p.poll()}")
        for line in p.stderr:
            print(line)

        if p.poll() != 0:
            raise Exception("cd-hit-est-2d did not complete successfully")

    @staticmethod
    def run_cd_hit_est(ncrna_cd_hit_est_2d_fasta: Path, ncrna_cd_hit_est_fasta: Path,
                          cpus: int=1, memory: int=1_000) -> None:
        print("\nStarting cd-hit-est\n(This may take a while)")
        cd_hit_command = ["cd-hit-est",
                          "-i", f"{ncrna_cd_hit_est_2d_fasta}",
                          "-o", f"{ncrna_cd_hit_est_fasta}",
                          "-c", "0.90",
                          "-T", f"{cpus}",
                          "-M", f"{memory}"]

        p = subprocess.Popen(cd_hit_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        while p.poll() is None and (line := p.stdout.readline()) != "":
            print(line.strip())
        p.wait()
        print(f"Exit code: {p.poll()}")
        for line in p.stderr:
            print(line)

        if p.poll() != 0:
            raise Exception("cd-hit-est did not complete successfully")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-sqlite_db", type=str, required=True)
    parser.add_argument("-transcripts_fasta", type=str, required=True)
    parser.add_argument("-ncrna_dir", type=str, required=True)
    parser.add_argument("-cds_fasta", type=str, required=True)
    parser.add_argument("-cpus", type=int, default=1, required=False)
    parser.add_argument("-memory", type=int, default=1_000, required=False)
    args = parser.parse_args()

    chm = CdHitManager(Path(args.sqlite_db), Path(args.transcripts_fasta), Path(args.ncrna_dir), Path(args.cds_fasta))
    chm.run(args.cpus, args.memory)