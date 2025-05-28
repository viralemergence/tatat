from argparse import ArgumentParser
from pathlib import Path
import sqlite3
import subprocess

class FastqPathManager:
    def __init__(self, fastq_dir: Path, sqlite_db : Path, uid: str, outdir: Path) -> None:
        fastq_file_names = self.extract_file_names(sqlite_db, uid)
        self.uid_fastq_paths = self.construct_uid_fastq_paths(fastq_dir, fastq_file_names)
        self.output_fastq_paths = self.set_output_paths(self.uid_fastq_paths, outdir, uid)

    @staticmethod
    def extract_file_names(sqlite_db: Path, uid: str) -> list[str]:
        with sqlite3.connect(sqlite_db) as connection:
            cursor = connection.cursor()
            sql_query = ("SELECT uid,r1_reads,r2_reads "
                         "FROM samples "
                         f"WHERE uid = '{uid}'")
            cursor.execute(sql_query)
            return {row[0]: [row[1], row[2]] for row in cursor.fetchall()}[uid]

    @staticmethod
    def construct_uid_fastq_paths(fastq_dir: Path, fastq_file_names: list[str]) -> list[Path]:
        uid_fastq_paths = []
        for file_name in fastq_file_names:
            if file_name is None:
                continue
            uid_fastq_path = fastq_dir / file_name
            uid_fastq_paths.append(uid_fastq_path)
        return uid_fastq_paths

    @staticmethod
    def set_output_paths(input_fastqs: list[Path], outdir: Path, uid: str) -> list[Path]:
        output_fastqs = []
        for i, fastq in enumerate(input_fastqs, 1):
            out_file_name = f"{uid}_r{i}_fastp.fastq.gz"
            outpath = outdir / out_file_name
            output_fastqs.append(outpath)
        return output_fastqs

class FastpManager:
    def __init__(self, input_fastqs: list[Path], output_fastqs: list[Path]) -> None:
        self.input_fastqs = input_fastqs
        self.output_fastqs = output_fastqs

    def run_fastp(self, r1_adapter: str, r2_adapter: str, cpus: int) -> None:
        if len(self.input_fastqs) == 2:
            fastp_command = ["fastp", "-V",
                            "-i", self.input_fastqs[0], "-I", self.input_fastqs[1],
                            "-o", self.output_fastqs[0], "-O", self.output_fastqs[1],
                            "--adapter_sequence", r1_adapter,
                            "--adapter_sequence_r2", r2_adapter,
                            "--thread", f"{cpus}"]
        if len(self.input_fastqs) == 1:
            fastp_command = ["fastp", "-V",
                            "-i", self.input_fastqs[0],
                            "-o", self.output_fastqs[0],
                            "--adapter_sequence", r1_adapter,
                            "--thread", f"{cpus}"]

        p = subprocess.Popen(fastp_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        while p.poll() is None and (line := p.stderr.readline()) != "":
            print(line.strip())
        p.wait()
        print(f"Exit code: {p.poll()}")

        if p.poll() != 0:
            raise Exception("Fastp did not complete successfully")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-fastq_dir", type=str, required=True)
    parser.add_argument("-outdir", type=str, required=True)
    parser.add_argument("-sqlite_db", type=str, required=True)
    parser.add_argument("-uid", type=str, required=True)
    parser.add_argument("-r1_adapter", type=str, required=True)
    parser.add_argument("-r2_adapter", type=str, required=True)
    parser.add_argument("-cpus", type=int, default=3, required=False)
    args = parser.parse_args()

    fpm = FastqPathManager(Path(args.fastq_dir), Path(args.sqlite_db), args.uid, Path(args.outdir))

    fm = FastpManager(fpm.uid_fastq_paths, fpm.output_fastq_paths)
    fm.run_fastp(args.r1_adapter, args.r2_adapter, args.cpus)