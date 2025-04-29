from argparse import ArgumentParser
from csv import DictWriter
from pathlib import Path
from typing import Any, Iterator, TextIO

class DeNovoAssemblyManager:
    def __init__(self, assembly_fasta_dir: Path, merged_path: Path, metadata_path: Path) -> None:
        self.assembly_fasta_dir = assembly_fasta_dir
        self.merged_path = merged_path
        self.metadata_path = metadata_path

    def run(self) -> None:
        fasta_files = self.get_file_list(self.assembly_fasta_dir)
        print(f"A total of {len(fasta_files)} files detected")

        with self.merged_path.open("w") as merged_outhandle, \
            self.metadata_path.open("w") as metadata_outhandle:

            metadata_fields = ["sequence_id", "sample_uid", "length", "kmer_coverage"]
            metadata_writer = DictWriter(metadata_outhandle, fieldnames=metadata_fields)
            metadata_writer.writeheader()
            
            sequence_id = 0
            for fasta_file in fasta_files:
                sample_uid = fasta_file.stem
                print(f"Starting on sample: {sample_uid}")
                for fasta_seq in self.fasta_chunker(fasta_file):
                    sequence_id += 1
                    
                    self.write_renamed_fasta_seq(sequence_id, merged_outhandle, fasta_seq[1:])

                    fasta_header = fasta_seq[0].replace(">", "")
                    self.write_assembly_metadata(fasta_header, sequence_id, sample_uid, metadata_fields, metadata_writer)

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

    @classmethod
    def write_assembly_metadata(cls, fasta_header: str, seq_id: int, sample_uid: str, metadata_fields: list[str], metadata_writer: DictWriter) -> None:
        assembly_info = cls.extract_assembly_info(fasta_header)
        metadata_values = [seq_id, sample_uid, assembly_info["length"], assembly_info["cov"]]
        metadata = {key: value for key, value in zip(metadata_fields, metadata_values)}
        metadata_writer.writerow(metadata)

    @staticmethod
    def extract_assembly_info(fasta_header: str) -> dict[Any]:
        assembly_info = fasta_header.split("_")
        assembly_dict = dict(zip(assembly_info[::2], assembly_info[1::2]))
        assembly_dict["length"] = int(assembly_dict["length"])
        assembly_dict["cov"] = round(float(assembly_dict["cov"]), 2)
        return assembly_dict

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-assembly_fasta_dir", type=str, required=True)
    parser.add_argument("-merged_path", type=str, required=True)
    parser.add_argument("-metadata_path", type=str, required=True)
    args = parser.parse_args()

    dnam = DeNovoAssemblyManager(Path(args.assembly_fasta_dir), Path(args.merged_path), Path(args.metadata_path))
    dnam.run()