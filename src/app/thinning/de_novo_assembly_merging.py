from argparse import ArgumentParser
from pathlib import Path
from shutil import copyfileobj
import subprocess
from typing import Iterator

class DeNovoAssemblyMerger:
    def __init__(self, assembly_fasta_dir: Path, outdir: Path, cpus: int, memory: int, concat_size: int) -> None:
        self.assembly_fasta_dir = assembly_fasta_dir
        self.outdir = outdir
        self.cpus = cpus
        self.memory = memory
        self.concat_size = concat_size

    def run(self) -> None:
        fasta_files = self.get_file_list(self.assembly_fasta_dir)

        base_file = fasta_files[0]
        i = 1
        for fasta_files_chunk in self.chunk_concat_files(fasta_files[1:], self.concat_size-1):
            i += (self.concat_size-1)
            concat_file_path = self.set_concat_path(self.outdir, "concat", i, base_file.suffix)

            fasta_files = [base_file] + fasta_files_chunk
            self.concat_files(fasta_files, concat_file_path)

            merged_assembly_file_path =self.set_merged_assembly_path(concat_file_path)
            self.cd_hit_est_merge(concat_file_path, merged_assembly_file_path, self.memory, self.cpus)

            concat_file_path.unlink()
            if i != self.concat_size:
                base_file.unlink()
            base_file = merged_assembly_file_path

    @staticmethod
    def get_file_list(directory: Path) -> list[Path]:
        return sorted([file for file in directory.iterdir()])
    
    @staticmethod
    def chunk_concat_files(files: list[Path], chunk_size: int) -> Iterator[list[Path]]:
        for i in range(0, len(files), chunk_size):
            yield files[i:i+chunk_size]

    @staticmethod
    def set_concat_path(outdir: Path, file_prefix: str, i: int, file_suffix: str) -> Path:
        return outdir / f"{file_prefix}_{i}{file_suffix}"

    @classmethod
    def concat_files(cls, file_paths: list[Path], outfile_path: Path) -> None:
        with outfile_path.open("wb") as outhandle:
            for file in file_paths:
                with file.open("rb") as inhandle:
                    copyfileobj(inhandle, outhandle)

    @staticmethod
    def set_merged_assembly_path(concat_file_path: Path) -> Path:
        return concat_file_path.parent / f"{concat_file_path.stem}_merged{concat_file_path.suffix}"

    @classmethod
    def cd_hit_est_merge(cls, concat_fasta_path: Path, merged_fasta_path: Path, memory: int, cpus: int) -> None:
        cd_hit_est_command = ["cd-hit-est", "-i", f"{concat_fasta_path}", "-o", f"{merged_fasta_path}",
                              "-c", "0.99", "-M", f"{memory}", "-T", f"{cpus}"]
        p = subprocess.Popen(cd_hit_est_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        while p.poll() is None and (line := p.stdout.readline()) != "":
            print(line.strip())
        p.wait()

        if p.poll() != 0:
            raise Exception("CD-HIT-EST did not complete successfully")
        
        cls.remove_cluster_file(merged_fasta_path)
        
    @staticmethod
    def remove_cluster_file(merged_fasta_path: Path) -> None:
        cluster_path = merged_fasta_path.parent / f"{merged_fasta_path.name}.clstr"
        cluster_path.unlink()

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-assembly_fasta_dir", type=str, required=True)
    parser.add_argument("-outdir", type=str, required=True)
    parser.add_argument("-cpus", type=int, required=False, default=1)
    parser.add_argument("-mem", type=int, required=False, default=1_000)
    parser.add_argument("-concat_size", type=int, required=False, default=2)
    args = parser.parse_args()

    dnam = DeNovoAssemblyMerger(Path(args.assembly_fasta_dir), Path(args.outdir), args.cpus, args.mem, args.concat_size)
    dnam.run()
