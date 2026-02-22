"""
JSON SCHEMA:
ames_mutagenic: bool
"""

from multiprocessing import Pool
from typing import Iterator

import numpy as np
import pandas as pd
from sqlalchemy.orm import Session

from drug_db.models.schemas import DrugSchema, RecordSchema, ScraperPayload
from drug_db.scraper.adapters import BaseAdapter
from drug_db.scraper.base import BaseFileScraper


def _process_worker(raw_data: dict) -> dict:
    smiles = raw_data.get("Canonical_Smiles")
    activity = raw_data.get("Activity")

    if smiles is None or activity is None:
        return {
            "status": "skipped",
            "raw_data": raw_data,
            "error": "Missing SMILES or Activity",
        }

    try:
        drug = DrugSchema(smiles=smiles, generic_name="Unknown")
        is_mutagenic = bool(int(activity))

        return {
            "status": "success",
            "drug": drug.model_dump(),
            "record_json": {"ames_mutagenic": is_mutagenic},
        }
    except Exception as e:
        return {"status": "error", "raw_data": raw_data, "error": str(e)}


class AmesAdapter(BaseAdapter):
    def __init__(self):
        self.skipped = []

    def standardise(self, row_data: dict) -> ScraperPayload | None:
        status = row_data.get("status")

        if status == "skipped":
            self.skipped.append(
                {"raw_data": row_data.get("raw_data"), "reason": row_data.get("error")}
            )
            return None

        if status == "success":
            return ScraperPayload(
                drug=DrugSchema(**row_data["drug"]),
                record=RecordSchema(data_json=row_data["record_json"]),
            )

        if status == "error":
            self.skipped.append(
                {"raw_data": row_data.get("raw_data"), "error": row_data.get("error")}
            )
            return None

        return None


class AmesScraper(BaseFileScraper):
    def __init__(
        self,
        session: Session,
        file_path: str,
        version: str,
        run_skip: bool = False,
    ):
        super().__init__(session, "Ames Mutagenicity", version, file_path)
        self.run_skip = run_skip

    def _get_adapter(self) -> BaseAdapter:
        return AmesAdapter()

    def _extract(self) -> Iterator[dict]:
        reader = pd.read_csv(self.file_path, chunksize=1000)

        def row_stream(reader_obj):
            for chunk in reader_obj:
                chunk = chunk.replace({np.nan: None})
                for row in chunk.to_dict(orient="records"):
                    yield row

        with Pool() as pool:
            results = pool.imap_unordered(
                _process_worker, row_stream(reader), chunksize=20
            )

            for result in results:
                if result:
                    yield result


if __name__ == "__main__":
    import sys
    from pathlib import Path

    from drug_db.database.manager import db

    BASE_DIR = Path(__file__).parent.parent.parent.resolve()
    DATA_DIR = BASE_DIR / ".." / "data"
    data_path = DATA_DIR / "ames_muta_v2.csv"

    db.create_tables()
    session = db.get_session()
    try:
        with session:
            scraper = AmesScraper(db.get_session(), data_path.resolve(), "v2")
            scraper.run()
    except KeyboardInterrupt:
        sys.exit(0)
    finally:
        session.close()
