from argparse import ArgumentParser
from csv import DictWriter, reader
from pathlib import Path
from typing import Any, Iterator, TextIO

class DeNovoAssemblyManager:
    def __init__(self, assembly_fasta_dir: Path, outdir: Path, merged_name: str, metadata_name: str) -> None:
        self.assembly_fasta_dir = assembly_fasta_dir
        self.merged_path = outdir / f"{merged_name}"
        self.metadata_path = outdir / f"{metadata_name}"

    def run(self) -> None:
        fasta_files = self.get_file_list(self.assembly_fasta_dir)
        print(f"A total of {len(fasta_files)} files detected")

        with self.merged_path.open("w") as merged_outhandle, \
            self.metadata_path.open("w") as metadata_outhandle:

            metadata_fields = ["sequence_id", "sample_id", "length", "kmer_coverage"]
            metadata_writer = DictWriter(metadata_outhandle, fieldnames=metadata_fields)
            metadata_writer.writeheader()
            
            i = 0
            for fasta_file in fasta_files:
                sample_id = fasta_file.stem
                print(f"Starting on sample: {sample_id}")
                for fasta_seq in self.fasta_chunker(fasta_file):
                    i += 1
                    seq_id = f"seq_{i}"
                    
                    self.write_renamed_fasta_seq(seq_id, merged_outhandle, fasta_seq[1:])

                    fasta_header = fasta_seq[0].replace(">", "")
                    self.write_assembly_metadata(fasta_header, seq_id, sample_id, metadata_fields, metadata_writer)

    @staticmethod
    def get_file_list(directory: Path) -> list[Path]:
        return sorted([file for file in directory.iterdir()])

    @staticmethod
    def fasta_chunker(fasta_path: Path) -> Iterator[list[str]]:
        fasta_seq = []
        first_chunk = True
        with fasta_path.open() as inhandle:
            reader_iterator = reader(inhandle)
            for line in reader_iterator:
                line = line[0]
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
    def write_renamed_fasta_seq(seq_id: str, outhandle: TextIO, fasta_seq: list[str]) -> None:
        # NOTE: fasta_seq should JUST be sequence, not including header
        # Otherwise, new and old header will both be written
        new_header = f">{seq_id}"
        outhandle.write(new_header + "\n")
        for line in fasta_seq:
            outhandle.write(line + "\n")

    @classmethod
    def write_assembly_metadata(cls, fasta_header: str, seq_id: str, sample_id: str, metadata_fields: list[str], metadata_writer: DictWriter) -> None:
        assembly_info = cls.extract_assembly_info(fasta_header)
        metadata_values = [seq_id, sample_id, assembly_info["length"], assembly_info["cov"]]
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
    parser.add_argument("-outdir", type=str, required=True)
    parser.add_argument("-merged_name", type=str, required=True)
    parser.add_argument("-metadata_name", type=str, required=True)
    args = parser.parse_args()

    dnam = DeNovoAssemblyManager(Path(args.assembly_fasta_dir), Path(args.outdir),
                                 args.merged_name, args.metadata_name)
    dnam.run()