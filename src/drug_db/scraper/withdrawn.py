from typing import Iterator, List, NamedTuple

import numpy as np
import pandas as pd
from pydantic import BaseModel, Field, field_validator

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
    "inchi_key": "str",
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


class WithdrawnDbScraper(BaseFileScraper):
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

        for chunk in reader:
            chunk = chunk.replace({np.nan: None})
            for row in chunk.itertuples(index=False, name="Row"):
                yield row

    def _get_adapter(self) -> WithdrawnDbAdapter:
        return WithdrawnDbAdapter()


class WithdrawnDbAdapter(BaseAdapter):
    def standardise(self, row_data: NamedTuple) -> ScraperPayload | None:
        raw_data = row_data._asdict()
        try:
            source_json = WithdrawnSourceData(**raw_data)

            valid_data = source_json.model_dump(exclude_none=True)

            if "withdrawn_id" in valid_data:
                del valid_data["withdrawn_id"]

            if not valid_data:
                return None

            drug_data = DrugSchema(
                inchi_key=raw_data.get("inchi_key"),
                generic_name=raw_data.get("drugname"),
                inchi=raw_data.get("inchi"),
                smiles=raw_data.get("smiles"),
                chem_formula=raw_data.get("formula"),
                mol_weight=raw_data.get("molwt"),
                drugbank_id=raw_data.get("DrugbankID"),
                chembl_id=raw_data.get("chemblid"),
            )

            return ScraperPayload(
                drug=drug_data,
                record=RecordSchema(
                    data_json=source_json.model_dump(exclude_none=True)
                ),
            )
        except ValueError as e:
            print(f"ERROR: invalid data cannot be mapped: {e}")
            return None
        except Exception as e:
            print(f"ERROR: {e}")
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
            scraper = WithdrawnDbScraper(
                db.get_session(), "Withdrawn", "2.0", data_path.resolve()
            )

            scraper.run()
    except KeyboardInterrupt:
        sys.exit(0)
    finally:
        session.close()
