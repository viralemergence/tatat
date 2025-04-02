from argparse import ArgumentParser
from csv import DictReader, reader
from pathlib import Path
from typing import Any, Iterator

class CdsAaFastaManager:
    def __init__(self, assembly_fasta: Path, transcript_metadata: Path, cds_metadata: Path, cds_fasta: Path) -> None:
        self.assembly_fasta = assembly_fasta
        self.transcript_metadata = transcript_metadata
        self.cds_metadata = cds_metadata
        self.cds_fasta = cds_fasta

    def run(self) -> None:
        transcript_cds_id_mapping = self.extract_transcript_cds_id_mapping(self.transcript_metadata)
        print(len(transcript_cds_id_mapping))
        cds_positions = self.extract_cds_positions(self.cds_metadata)
        print(len(cds_positions))

        if self.cds_fasta:
            self.extract_and_write_cds(transcript_cds_id_mapping, cds_positions, self.assembly_fasta, self.cds_fasta)

    @staticmethod
    def extract_transcript_cds_id_mapping(transcript_metadata: Path) -> dict[list[str]]:
        print("Starting transcript id to cds id mapping")
        transcript_cds_id_mapping = dict()

        with transcript_metadata.open() as inhandle:
            reader = DictReader(inhandle)
            for data in reader:
                evigene_pass = data["evigene_pass"]
                if evigene_pass != "True":
                    continue

                cds_ids = data["cds_ids"].split(";")
                if cds_ids[0] == "":
                    continue

                transcript_id = data["sequence_id"]
                transcript_cds_id_mapping[transcript_id] = cds_ids
        return transcript_cds_id_mapping
    
    @staticmethod
    def extract_cds_positions(cds_metadata: Path) -> dict[Any]:
        print("Starting cds positions extraction")
        cds_positions = dict()

        with cds_metadata.open() as inhandle:
            reader = DictReader(inhandle)
            for data in reader:
                cds_id = data["cds_id"]
                start = int(data["start"]) - 1
                end = int(data["end"])
                strand = data["strand"]

                cds_data = {"start": start, "end": end, "strand": strand}

                cds_positions[cds_id] = cds_data
        return cds_positions

    @classmethod
    def extract_and_write_cds(cls, transcript_cds_id_mapping: dict[list[str]], cds_positions: dict[Any],
                              assembly_fasta:Path, cds_fasta: Path) -> None:
        print("Starting cds extraction and writing\n(This may take awhile)")
        with cds_fasta.open("w") as cds_outhandle:
            for fasta_seq in cls.fasta_chunker(assembly_fasta):
                transcript_id = fasta_seq[0][1:]
                try:
                    cds_ids = transcript_cds_id_mapping[transcript_id]
                except KeyError:
                    continue

                transcript_seq = "".join(fasta_seq[1:])

                for cds_id in cds_ids:
                    cds_seq = cls.extract_cds_sequence(cds_positions[cds_id], transcript_seq)
                    cds_outhandle.write(f">{cds_id}\n")
                    cds_outhandle.write(f"{cds_seq}\n")

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

    @classmethod
    def extract_cds_sequence(cls, cds_position: dict[Any], transcript_seq: str) -> str:
        cds_start = cds_position["start"]
        cds_end = cds_position["end"]
        cds_seq = transcript_seq[cds_start: cds_end]

        if cds_position["strand"] == "-":
            cds_seq = cls.reverse_translate_dna(cds_seq)
        return cds_seq

    @staticmethod
    def reverse_translate_dna(dna_sequence: str) -> str:
        translation_mapping = str.maketrans("ATCG", "TAGC")
        return dna_sequence.translate(translation_mapping)[::-1]

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-assembly_fasta", type=str, required=True)
    parser.add_argument("-transcript_metadata", type=str, required=True)
    parser.add_argument("-cds_metadata", type=str, required=True)
    parser.add_argument("-cds_fasta", type=str, required=True)
    args = parser.parse_args()

    cafm = CdsAaFastaManager(Path(args.assembly_fasta), Path(args.transcript_metadata), Path(args.cds_metadata), Path(args.cds_fasta))
    print("\nStarting CDS/AA Fasta Manager")
    cafm.run()
    print("\nFinished")