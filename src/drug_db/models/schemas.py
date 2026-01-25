from typing import Any

import pubchempy as pcp
from pydantic import BaseModel, Field, field_validator, model_validator
from rdkit import Chem


class DrugSchema(BaseModel):
    inchi_key: str | None = Field(
        default=None, pattern=r"^[A-Z]{14}-[A-Z]{8}[SN]A-[A-Z]$"
    )
    generic_name: str | None
    inchi: str | None
    chem_formula: str | None
    mol_weight: float | None
    smiles: str | None = None
    drugbank_id: str | None = None
    chembl_id: str | None = None
    synonyms: list[str] = Field(default_factory=list)

    @model_validator(mode="before")
    @classmethod
    def missing_ids(cls, data: dict):
        if not isinstance(data, dict):
            return data

        inchi_key = data.get("inchi_key")
        inchi = data.get("inchi")
        smiles = data.get("smiles")
        name = data.get("generic_name")

        if not inchi_key and smiles:
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    if not inchi:
                        inchi = Chem.MolToInchi(mol)
                        data["inchi"] = inchi
                    new_key = Chem.InchiToInchiKey(inchi)
                    data["inchi_key"] = new_key
            except Exception:
                pass

        missing_ids = not all(
            [data.get("inchi_key"), data.get("drugbank_id"), data.get("chembl_id")]
        )

        if missing_ids and name:
            try:
                results = pcp.get_compounds(name, "name")

                if results:
                    top_result = results[0]
                    # perhaps require human intervention??? entirely relying on search algo + name being correct

                    if not data.get("inchi_key"):
                        data["inchi_key"] = top_result.inchikey
                    if not data.get("inchi"):
                        data["inchi"] = top_result.inchi
                    if not data.get("smiles"):
                        data["smiles"] = (
                            top_result.isomeric_smiles
                        )  # Fixed deprecation warning
                    if not data.get("mol_weight"):
                        data["mol_weight"] = top_result.molecular_weight
                    if not data.get("chem_formula"):
                        data["chem_formula"] = top_result.molecular_formula

                    # query another db for missing chembl_id or drugbank_id
            except Exception as e:
                print(f"ERROR: pubchem lookup failed: {e}")
        return data

    @model_validator(mode="after")
    def verify_id(self):
        if not self.inchi_key:
            raise ValueError(f"ERROR: Could not get InChIKey for {self.generic_name}")
        return self

    @field_validator("synonyms", mode="before")
    @classmethod
    def convert_none_to_list(cls, value: any) -> list:
        if value is None:
            return []
        return value


class RecordSchema(BaseModel):
    data_json: dict[str, Any] = Field(default_factory=dict)


class ScraperPayload(BaseModel):
    drug: DrugSchema
    record: RecordSchema
