from argparse import ArgumentParser
from constants import CODON_TO_AMINO_ACID
from csv import DictReader, reader
from json import load
from pathlib import Path
from typing import Any, Iterator, Union

class CdsAaFastaManager:
    def __init__(self, assembly_fasta: Path, transcript_metadata: Path, cds_metadata: Path,
                 extraction_fields: Path, cds_fasta: Union[None, Path], aa_fasta: Union[None, Path]) -> None:
        self.assembly_fasta = assembly_fasta
        self.transcript_metadata = transcript_metadata
        self.cds_metadata = cds_metadata
        self.extraction_fields = self.extract_extraction_fields(extraction_fields)
        self.cds_fasta = cds_fasta
        self.aa_fasta = aa_fasta

    @staticmethod
    def extract_extraction_fields(json_path: Path) -> dict[dict[str]]:
        with json_path.open() as inhandle:
            return load(inhandle)

    def run(self, codon_to_aa: dict[str]=None, add_gene_name: True=None) -> None:
        # Extract a dictionary of transcript ids mapping to cds ids (sometimes more than one cds id per transcript id)
        transcript_cds_id_mapping = self.extract_transcript_cds_id_mapping(self.transcript_metadata, self.extraction_fields)
        print(len(transcript_cds_id_mapping))

        # If cds metadata extraction fields are present, extract cds ids that meet the criteria,
        # then remove any transcript ids whose cds ids do not meet that criteria
        if self.extraction_fields["cds_metadata"] is not None:
            filtered_cds_ids = self.extract_filtered_cds_ids(self.cds_metadata, self.extraction_fields)
            transcript_cds_id_mapping = self.filter_transcript_cds_id_mapping(transcript_cds_id_mapping, filtered_cds_ids)
            print(len(transcript_cds_id_mapping))

        if add_gene_name:
            cds_gene_mapping = self.extract_cds_gene_mapping(self.cds_metadata)
        else:
            cds_gene_mapping = None

        # Extract cds_position info
        cds_positions = self.extract_cds_positions(self.cds_metadata)
        print(len(cds_positions))

        # Using the remaining transcript ids write cds and/or aa fasta files
        if (self.cds_fasta) and (not self.aa_fasta):
            self.extract_and_write_cds(transcript_cds_id_mapping, cds_positions, self.assembly_fasta, self.cds_fasta, cds_gene_mapping)
        elif (not self.cds_fasta) and (self.aa_fasta):
            self.extract_and_write_aa(transcript_cds_id_mapping, cds_positions, self.assembly_fasta, self.aa_fasta, codon_to_aa, cds_gene_mapping)
        elif self.cds_fasta and self.aa_fasta:
            self.extract_and_write_cds_aa(transcript_cds_id_mapping, cds_positions, self.assembly_fasta,
                                          self.cds_fasta, self.aa_fasta, codon_to_aa, cds_gene_mapping)
        else:
            raise Exception("Neither output cds nor aa fasta paths detected\nNo calculation performed")

    @staticmethod
    def extract_transcript_cds_id_mapping(transcript_metadata: Path, extraction_fields: dict[dict[str]]) -> dict[list[str]]:
        transcript_extraction_fields = extraction_fields["transcriptome_metadata"]
        print("Starting transcript id to cds id mapping")
        transcript_cds_id_mapping = dict()

        with transcript_metadata.open() as inhandle:
            reader = DictReader(inhandle)
            for data in reader:
                skip = False
                for field, value in transcript_extraction_fields.items():
                    if data[field] != value:
                        skip = True
                        break
                if skip:
                    continue

                cds_ids = data["cds_ids"].split(";")
                if cds_ids[0] == "":
                    continue

                transcript_id = data["sequence_id"]
                transcript_cds_id_mapping[transcript_id] = cds_ids
        return transcript_cds_id_mapping

    @staticmethod
    def extract_filtered_cds_ids(cds_metadata: Path, extraction_fields: dict[dict[str]]) -> set[str]:
        cds_extraction_fields = extraction_fields["cds_metadata"]
        print("Starting filtered cds id mapping")
        filtered_cds_ids = set()

        with cds_metadata.open() as inhandle:
            reader = DictReader(inhandle)
            for data in reader:
                skip = False
                for field, value in cds_extraction_fields.items():
                    if data[field] != value:
                        skip = True
                        break
                if skip:
                    continue

                cds_id = data["cds_id"]
                filtered_cds_ids.add(cds_id)
        return filtered_cds_ids

    @staticmethod
    def filter_transcript_cds_id_mapping(transcript_cds_id_mapping: dict[list[str]], filtered_cds_ids: set[str]) -> dict[list[str]]:
        filtered_transcript_cds_id_mapping = dict()

        for transcript_id, cds_ids in transcript_cds_id_mapping.items():
            cds_ids_keepers = []
            for cds_id in cds_ids:
                if cds_id in filtered_cds_ids:
                    cds_ids_keepers.append(cds_id)
            if len(cds_ids_keepers) < 1:
                continue
            filtered_transcript_cds_id_mapping[transcript_id] = cds_ids_keepers
        return filtered_transcript_cds_id_mapping

    @staticmethod
    def extract_cds_gene_mapping(cds_metadata: Path) -> dict[str]:
        print("Starting cds to gene name mapping")
        cds_gene_mapping = dict()

        with cds_metadata.open() as inhandle:
            reader = DictReader(inhandle)
            for data in reader:
                cds_id = data["cds_id"]
                gene = data["gene"]
                cds_gene_mapping[cds_id] = gene
        return cds_gene_mapping

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
                              assembly_fasta:Path, cds_fasta: Path, cds_gene_mapping: Union[None, dict[str]]) -> None:
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
                    if cds_gene_mapping:
                        gene = cds_gene_mapping[cds_id]
                        cds_outhandle.write(f">{cds_id};{gene}\n")
                    else:
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

    @classmethod
    def extract_and_write_aa(cls, transcript_cds_id_mapping: dict[list[str]], cds_positions: dict[Any],
                              assembly_fasta:Path, aa_fasta: Path, codon_to_aa: dict[str],
                              cds_gene_mapping: Union[None, dict[str]]) -> None:
        print("Starting cds extraction, aa translation, and aa writing\n(This may take awhile)")

        with aa_fasta.open("w") as aa_outhandle:
            for fasta_seq in cls.fasta_chunker(assembly_fasta):
                transcript_id = fasta_seq[0][1:]
                try:
                    cds_ids = transcript_cds_id_mapping[transcript_id]
                except KeyError:
                    continue

                transcript_seq = "".join(fasta_seq[1:])

                for cds_id in cds_ids:
                    cds_seq = cls.extract_cds_sequence(cds_positions[cds_id], transcript_seq)
                    try:
                        aa_seq = cls.translate_cds_to_aa(cds_seq, codon_to_aa)
                    except KeyError:
                        continue
                    if cds_gene_mapping:
                        gene = cds_gene_mapping[cds_id]
                        aa_outhandle.write(f">{cds_id};{gene}\n")
                    else:
                        aa_outhandle.write(f">{cds_id}\n")
                    aa_outhandle.write(f"{aa_seq}\n")

    @classmethod
    def translate_cds_to_aa(cls, dna_sequence: str, codon_to_aa: dict[str]) -> str:
        return "".join([codon_to_aa[codon] for codon in cls.codon_chunker(dna_sequence)])

    @staticmethod
    def codon_chunker(dna_sequence: str) -> Iterator[str]:
        codon_size = 3
        for i in range(0, len(dna_sequence), codon_size):
            yield dna_sequence[i:i+codon_size]

    @classmethod
    def extract_and_write_cds_aa(cls, transcript_cds_id_mapping: dict[list[str]], cds_positions: dict[Any],
                              assembly_fasta:Path, cds_fasta:Path, aa_fasta: Path, codon_to_aa: dict[str],
                              cds_gene_mapping: Union[None, dict[str]]) -> None:
        print("Starting cds extraction, aa translation, and cds and aa writing\n(This may take awhile)")

        with cds_fasta.open("w") as cds_outhandle, aa_fasta.open("w") as aa_outhandle:
            for fasta_seq in cls.fasta_chunker(assembly_fasta):
                transcript_id = fasta_seq[0][1:]
                try:
                    cds_ids = transcript_cds_id_mapping[transcript_id]
                except KeyError:
                    continue

                transcript_seq = "".join(fasta_seq[1:])

                for cds_id in cds_ids:
                    cds_seq = cls.extract_cds_sequence(cds_positions[cds_id], transcript_seq)
                    if cds_gene_mapping:
                        gene = cds_gene_mapping[cds_id]
                        cds_outhandle.write(f">{cds_id};{gene}\n")
                    else:
                        cds_outhandle.write(f">{cds_id}\n")
                    cds_outhandle.write(f"{cds_seq}\n")
                    try:
                        aa_seq = cls.translate_cds_to_aa(cds_seq, codon_to_aa)
                    except KeyError:
                        continue
                    if cds_gene_mapping:
                        gene = cds_gene_mapping[cds_id]
                        aa_outhandle.write(f">{cds_id};{gene}\n")
                    else:
                        aa_outhandle.write(f">{cds_id}\n")
                    aa_outhandle.write(f"{aa_seq}\n")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-assembly_fasta", type=str, required=True)
    parser.add_argument("-transcript_metadata", type=str, required=True)
    parser.add_argument("-cds_metadata", type=str, required=True)
    parser.add_argument("-extraction_fields", type=str, required=True)
    parser.add_argument("-cds_fasta", type=str, required=False)
    parser.add_argument("-aa_fasta", type=str, required=False)
    parser.add_argument("-add_gene_name", action="store_true", required=False)
    args = parser.parse_args()

    cds_fasta_path = Path(args.cds_fasta) if args.cds_fasta else None
    aa_fasta_path = Path(args.aa_fasta) if args.aa_fasta else None

    cafm = CdsAaFastaManager(Path(args.assembly_fasta), Path(args.transcript_metadata), Path(args.cds_metadata),
                             Path(args.extraction_fields), cds_fasta_path, aa_fasta_path)
    print("\nStarting CDS/AA Fasta Manager")
    cafm.run(CODON_TO_AMINO_ACID, args.add_gene_name)
    print("\nFinished")