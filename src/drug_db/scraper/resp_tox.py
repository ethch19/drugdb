"""
JSON SCHEMA:
respiratory_toxicity: bool
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
    smiles = raw_data.get("SMILES")
    label = raw_data.get("Label")

    if smiles is None or label is None:
        return {
            "status": "skipped",
            "raw_data": raw_data,
            "error": "Missing SMILES or Label",
        }

    try:
        drug = DrugSchema(smiles=smiles, generic_name=None)
        is_toxic = bool(int(label))

        return {
            "status": "success",
            "drug": drug.model_dump(),
            "record_json": {"respiratory_toxicity": is_toxic},
        }
    except Exception as e:
        return {"status": "error", "raw_data": raw_data, "error": str(e)}


class RespToxAdapter(BaseAdapter):
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


class RespToxScraper(BaseFileScraper):
    def __init__(
        self,
        session: Session,
        file_path: str,
        version: str = "1.0",
        run_skip: bool = False,
    ):
        super().__init__(session, "Respiratory Toxicity", version, file_path)
        self.run_skip = run_skip

    def _get_adapter(self) -> BaseAdapter:
        return RespToxAdapter()

    def _extract(self) -> Iterator[dict]:
        df = pd.read_excel(self.file_path, sheet_name="Table S1", skiprows=1)

        if pd.isna(df.iloc[0]["SMILES"]):
            df = df.iloc[1:].copy()

        df = df.replace({np.nan: None})
        records = df.to_dict(orient="records")

        with Pool() as pool:
            results = pool.imap_unordered(_process_worker, records, chunksize=20)

            for result in results:
                if result:
                    yield result


if __name__ == "__main__":
    import sys
    from pathlib import Path

    from drug_db.database.manager import db

    BASE_DIR = Path(__file__).parent.parent.parent.resolve()
    DATA_DIR = BASE_DIR / ".." / "data"
    data_path = DATA_DIR / "resp_tox_data.xlsx"

    db.create_tables()
    session = db.get_session()
    try:
        with session:
            scraper = RespToxScraper(session, data_path.resolve())
            scraper.run()
    except KeyboardInterrupt:
        sys.exit(0)
    finally:
        session.close()
