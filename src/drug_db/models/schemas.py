from typing import Any

import pubchempy as pcp
from pydantic import BaseModel, Field, field_validator, model_validator
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt, CalcMolFormula


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
    def prop_filler(cls, data: dict):
        if not isinstance(data, dict):
            return data

        inchi_key = data.get("inchi_key")
        inchi = data.get("inchi")
        name = data.get("generic_name")
        mol_weight = data.get("mol_weight")
        chem_formula = data.get("chem_formula")
        smiles = data.get("smiles")

        checklist = [inchi, inchi_key, mol_weight, chem_formula]
        empty_fields = any(
            v is None or (isinstance(v, str) and not v.strip()) for v in checklist
        )

        # local compute from smiles or inchi/inchi_key
        if empty_fields:
            mol = None
            try:
                if smiles and isinstance(smiles, str) and smiles.strip():
                    mol = Chem.MolFromSmiles(smiles)
                if not mol and inchi and isinstance(inchi, str) and inchi.strip():
                    mol = Chem.MolFromInchi(inchi)

                if mol:
                    calculated = {
                        "inchi": Chem.MolToInchi(mol)[0],
                        "inchi_key": Chem.MolToInchiKey(mol),
                        "mol_weight": CalcExactMolWt(mol),
                        "chem_formula": CalcMolFormula(mol),
                        "smiles": Chem.MolToSmiles(mol),
                    }

                    for field, new_val in calculated.items():
                        cur_val = data.get(field)
                        if cur_val is None or (
                            isinstance(cur_val, str) and not cur_val.strip()
                        ):
                            data[field] = new_val
                    return data
            except Exception:
                pass

        # api calls with no smiles + only name
        if not data.get("inchi_key") and name:
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
                        data["smiles"] = top_result.smiles
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

    @field_validator("smiles", mode="before")
    @classmethod
    def canonical_smiles(cls, value: any) -> str:
        if value is None or not str(value).strip():
            return ""

        mol = Chem.MolFromSmiles(value)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)


class RecordSchema(BaseModel):
    data_json: dict[str, Any] = Field(default_factory=dict)


class ScraperPayload(BaseModel):
    drug: DrugSchema
    record: RecordSchema
