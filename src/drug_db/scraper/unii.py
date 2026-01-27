import pickle
import re
from datetime import datetime
from pathlib import Path

import pandas as pd


class UNIIResolver:
    def __init__(self):
        self.mapping = {}
        self.file_path = self._find_latest_file()

        if self.file_path:
            self.cache_path = self.file_path.with_suffix(".pkl")

            # 3. Load Cache if valid, else Process CSV
            if self._is_cache_valid():
                print(f"Loading cached UNII data from: {self.cache_path.name}")
                self._load_from_cache()
            else:
                print(f"Parsing raw CSV: {self.file_path.name}")
                self._load_data()
        else:
            print("ERROR: No UNII_Records_*.txt file found in data/ directory")

    def _find_latest_file(self) -> Path | None:
        repo_root = Path.cwd()
        data_dir = repo_root / "data"

        if not data_dir.exists():
            print(f"ERROR: Could not locate data directory at {data_dir}")
            return None

        files = list(data_dir.glob("UNII_Records_*.txt"))
        if not files:
            return None

        def parse_file_date(file_path: Path):
            match = re.search(r"UNII_Records_(.+)\.txt", file_path.name)
            if match:
                date_str = match.group(1)
                try:
                    return datetime.strptime(date_str, "%d%b%Y")
                except ValueError:
                    print(f"Skipping file with unparseable date: {file_path.name}")
                    return datetime.min
            return datetime.min

        latest_file = max(files, key=parse_file_date)
        return latest_file

    def _is_cache_valid(self) -> bool:
        return (
            self.cache_path.exists()
            and self.cache_path.stat().st_mtime >= self.file_path.stat().st_mtime
        )

    def _load_from_cache(self):
        try:
            with open(self.cache_path, "rb") as f:
                self.mapping = pickle.load(f)
        except Exception as e:
            print(f"Cache load failed ({e}), falling back to CSV.")
            self._load_data()

    def _load_data(self):
        try:
            cols = ["UNII", "INCHIKEY", "SMILES", "MF"]

            df = pd.read_csv(
                self.file_path, sep="\t", usecols=cols, dtype=str, on_bad_lines="skip"
            )

            df = df.rename(
                columns={
                    "UNII": "unii",
                    "INCHIKEY": "inchi_key",
                    "SMILES": "smiles",
                    "MF": "mol_formula",
                }
            )

            df = df.dropna(subset=["unii"])
            df = df.dropna(subset=["inchi_key", "smiles", "mol_formula"], how="all")
            df = df.where(pd.notnull(df), None)

            self.mapping = df.set_index("unii").to_dict(orient="index")
            with open(self.cache_path, "wb") as f:
                pickle.dump(self.mapping, f, protocol=pickle.HIGHEST_PROTOCOL)

        except Exception as e:
            print(f"ERROR: cannot load UNII file: {e}")

    def resolve(self, unii_id: str) -> dict | None:
        return self.mapping.get(unii_id)


unii_resolver = UNIIResolver()
