"""
JSON SCHEMA:
withdrawn_id: str
atc_code: str | None
pubchem_cid: str | None
ctd_id: str | None
inchikey_jchem: str | None


ld50: float | None
protox_toxclass: int | None

first_withdrawn_year: int | None
last_withdrawn_year: int | None
first_approval_year: int | None
withdrawal_countries: List[str] | None
toxicity_types: List[str] | None
"""

from multiprocessing import Pool
from typing import Iterator, List, NamedTuple

import numpy as np
import pandas as pd
from pydantic import BaseModel, Field, ValidationError, field_validator

from drug_db.database.manager import Session
from drug_db.models.schemas import DrugSchema, RecordSchema, ScraperPayload
from drug_db.scraper.base import BaseAdapter, BaseFileScraper

dtype_spec = {
    "id": "str",
    "withdrawnID": "str",
    "drugname": "str",
    "atc": "str",
    "cid": "str",  # PubChem compound id
    "chemblid": "str",
    "inchi": "str",
    "inchikey": "str",
    "smiles": "str",
    "ld50": "Float64",
    "toxclass": "Int8",  # protox-ii tox class (I to VI, according GHS)
    "firstwithdrawn": "Int64",
    "lastwithdrawn": "Int64",
    "firstapproval": "Int64",
    "firstdeath": "Int64",
    "toxtype": "str",  # comma-separated list (e.g: "neurological, dermatological")
    "reference": "str",
    "hba": "Int8",
    "hbd": "Int8",
    "nrbonds": "Int8",
    "formula": "str",
    "molwt": "Float64",
    "tpsa": "Float64",
    "maccsfp": "str",
    "dataset": "str",
    "inchikey_jchem": "str",
    "ctdid": "str",  # Comparative Toxicogenomics Database
    "withdrawalCountries": "str",  # semi-colon separated list (e.g: "EU;JPN;PER")
    "molfile": "str",
    "DrugbankID": "str",
    "fp24": "str",
}


class WithdrawnSourceData(BaseModel):
    # unique
    withdrawn_id: str = Field(validation_alias="withdrawnID")

    # untested identifiers
    atc_code: str | None = Field(None, validation_alias="atc")
    pubchem_cid: str | None = Field(None, validation_alias="cid")
    ctd_id: str | None = Field(None, validation_alias="ctdid")
    inchikey_jchem: str | None = Field(None, validation_alias="inchikey_jchem")

    # standardised
    ld50: float | None = Field(None, validation_alias="ld50")
    protox_toxclass: int | None = Field(None, validation_alias="toxclass")

    # to be standardised
    first_withdrawn_year: int | None = Field(None, validation_alias="firstwithdrawn")
    last_withdrawn_year: int | None = Field(None, validation_alias="lastwithdrawn")
    first_approval_year: int | None = Field(None, validation_alias="firstapproval")
    withdrawal_countries: List[str] | None = Field(
        None, validation_alias="withdrawalCountries"
    )
    toxicity_types: List[str] | None = Field(None, validation_alias="toxtype")

    @field_validator("toxicity_types", mode="before")
    @classmethod
    def parse_comma_separated(cls, value):
        if not value or not isinstance(value, str):
            return None
        return [x.strip() for x in value.split(",") if x.strip()]

    @field_validator("withdrawal_countries", mode="before")
    @classmethod
    def parse_semicolon_separated(cls, value):
        if not value or not isinstance(value, str):
            return None
        return [x.strip() for x in value.split(";") if x.strip()]


def _process_row(raw_data: dict) -> dict:
    try:
        source_model = WithdrawnSourceData(**raw_data)

        raw_key = raw_data.get("inchikey")
        if raw_key and isinstance(raw_key, str):
            raw_key = raw_key.strip()
            if raw_key == "":
                raw_key = None

        if not raw_key:
            return {
                "status": "skipped",
                "name": raw_data.get("drugname", "Unknown"),
                "raw_data": raw_data,
            }

        try:
            drug_entry = DrugSchema(
                inchi_key=raw_key,
                generic_name=raw_data.get("drugname"),
                inchi=raw_data.get("inchi"),
                smiles=raw_data.get("smiles"),
                chem_formula=raw_data.get("formula"),
                mol_weight=raw_data.get("molwt"),
                drugbank_id=raw_data.get("DrugbankID"),
                chembl_id=raw_data.get("chemblid"),
            )
        except ValidationError:
            return {
                "status": "skipped",
                "name": raw_data.get("drugname", "Unknown"),
                "raw_data": raw_data,
            }

        return {
            "status": "success",
            "drug": drug_entry.model_dump(),
            "record_json": source_model.model_dump(exclude_none=True),
        }

    except Exception as e:
        return {"status": "error", "error": str(e), "raw_data": raw_data}


class WithdrawnDbScraper(BaseFileScraper):
    def __init__(self, session: Session, version: str, file_path: str):
        super().__init__(session, "Withdrawn", version, file_path)

    def _extract(self) -> Iterator[pd.Series]:
        reader = pd.read_csv(
            self.file_path,
            chunksize=1000,
            dtype=dtype_spec,
            parse_dates=False,
            na_values=["NULL", "N/A", "null", "n/a", "", "..."],
            keep_default_na=True,
            engine="c",
        )

        def row_stream(reader):
            for chunk in reader:
                chunk = chunk.replace({np.nan: None})
                for row in chunk.to_dict(orient="records"):
                    yield row

        with Pool() as pool:
            print(f"Workers: {pool._processes}")
            results = pool.imap_unordered(
                _process_row, row_stream(reader), chunksize=20
            )

            for result in results:
                if result:
                    yield result

    def _get_adapter(self) -> WithdrawnDbAdapter:
        return WithdrawnDbAdapter()


class WithdrawnDbAdapter(BaseAdapter):
    def __init__(self):
        self.skipped = []

    def standardise(self, row_data: NamedTuple) -> ScraperPayload | None:
        status = row_data.get("status")

        if status == "skipped":
            self.skipped.append(
                {"name": row_data.get("name"), "raw_data": row_data.get("raw_data")}
            )
            return None

        if status == "success":
            source_json = row_data["record_json"]
            if "withdrawn_id" in source_json:
                del source_json["withdrawn_id"]

            return ScraperPayload(
                drug=DrugSchema(**row_data["drug"]),
                record=RecordSchema(data_json=source_json),
            )

        if status == "error":
            print(f"[WORKER ERROR] {row_data.get('error')}")
            return None

        return None


if __name__ == "__main__":
    import sys
    from pathlib import Path

    from drug_db.database.manager import db

    BASE_DIR = Path(__file__).parent.parent.parent.resolve()
    data_path = BASE_DIR / ".." / "data" / "withdrawn_2_dataset.csv"

    db.create_tables()
    print(data_path.resolve())
    session = db.get_session()
    try:
        with session:
            scraper = WithdrawnDbScraper(session, "2.0", data_path.resolve())

            scraper.run()
    except KeyboardInterrupt:
        sys.exit(0)
    finally:
        session.close()
