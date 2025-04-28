from argparse import ArgumentParser
from csv import DictReader
from pathlib import Path
from typing import Any

class SampleMetadataManager:
    def __init__(self, sample_metadata: Path, encoding: str) -> None:
        self.sample_metadata = self.extract_sample_metadata(sample_metadata, encoding)

    @staticmethod
    def extract_sample_metadata(sample_metadata_path: Path, encoding: str) -> list[dict[Any]]:
        with sample_metadata_path.open("r", encoding=encoding) as inhandle:
            reader = DictReader(inhandle)
            return [data for data in reader]

    def print_uid(self, array_index: int) -> None:
        print(self.sample_metadata[array_index]["uid"], flush=True)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-sample_metadata", type=str, required=True)
    parser.add_argument("-array_index", type=int, required=True)
    parser.add_argument("-encoding", type=str, required=False, default="utf-8")
    args = parser.parse_args()

    smm = SampleMetadataManager(Path(args.sample_metadata), args.encoding)
    smm.print_uid(args.array_index)