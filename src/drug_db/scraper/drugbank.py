from __future__ import annotations

import os
import time
from multiprocessing import Pool
from typing import Any, Dict, Iterator, List

from lxml import etree
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt, CalcMolFormula

from drug_db.database.manager import Session
from drug_db.models.schemas import DrugSchema, RecordSchema, ScraperPayload
from drug_db.scraper.base import BaseAdapter, BaseFileScraper
from drug_db.scraper.unii import unii_resolver


class Property(BaseModel):
    kind: str
    value: str
    source: str | None = None


class ExternalID(BaseModel):
    resource: str
    identifier: str


class DrugbankData(BaseModel):
    model_config = ConfigDict(validate_by_name=True)

    drugbank_id: str = Field(alias="drugbank-id")
    drugname: str = Field(alias="name")
    unii: str | None = None

    # optional primitive fields
    description: str | None = None
    toxicity: str | None = None
    moa: str | None = Field(default=None, alias="mechanism-of-action")

    # nested lists
    synonyms: List[str] = Field(default_factory=list)
    groups: List[str] = Field(default_factory=list)
    categories: List[Dict[str, Any]] = Field(default_factory=list)
    atc_codes: List[str] = Field(default_factory=list, alias="atc-codes")
    calculated_props: List[Property] = Field(
        default_factory=list, alias="calculated-properties"
    )
    experimental_props: List[Property] = Field(
        default_factory=list, alias="experimental-properties"
    )
    external_ids: List[ExternalID] = Field(
        default_factory=list, alias="external-identifiers"
    )

    # derived props
    inchi_key: str | None = None
    inchi: str | None = None
    smiles: str | None = None
    formula: str | None = None
    mol_weight: str | None = None
    chembl_id: str | None = None

    @field_validator("drugbank_id", mode="before")
    @classmethod
    def extract_id(cls, value: any) -> str:
        if isinstance(value, list):
            for vid in value:
                if isinstance(vid, dict) and vid.get("primary") == "true":
                    return vid.get("#text")
            first = value[0]
            return first.get("#text") if isinstance(first, dict) else str(first)

        if isinstance(value, dict):
            return value.get("#text")

        return str(value)

    @field_validator(
        "synonyms",
        "groups",
        "categories",
        "calculated_props",
        "experimental_props",
        "external_ids",
        mode="before",
    )
    @classmethod
    def unwrap_list(cls, value: any) -> list[any]:
        if value == "" or value is None:
            return []

        if isinstance(value, dict):
            found_list = False
            for k, v in value.items():
                if isinstance(v, list):
                    value = v
                    found_list = True
                    break

            if not found_list:
                found_child = False
                for k, v in value.items():
                    if isinstance(v, dict):
                        value = [v]
                        found_child = True
                        break

                if not found_child:
                    return []

        if not isinstance(value, list):
            value = [value]

        clean_list = []
        for item in value:
            if isinstance(item, dict):
                if "#text" in item and "kind" not in item:
                    clean_list.append(item["#text"])
                else:
                    clean_list.append(item)
            else:
                clean_list.append(item)
        return clean_list

    @field_validator("atc_codes", mode="before")
    @classmethod
    def decompose_atc_codes(cls, value: any) -> list[str]:
        if not value:
            return []

        if isinstance(value, dict):
            value = value.get("atc-code", [])

        if isinstance(value, dict):
            value = [value]

        if not isinstance(value, list):
            return []

        extracted_codes = []
        for item in value:
            if isinstance(item, dict):
                code = item.get("code")
                if code:
                    extracted_codes.append(code)
            elif isinstance(item, str):
                extracted_codes.append(item)

        return extracted_codes

    @model_validator(mode="after")
    def decompose_props(self) -> DrugbankData:
        prop_map = {
            "InChIKey": "inchi_key",
            "InChI": "inchi",
            "SMILES": "smiles",
            "Molecular Formula": "formula",
            "Molecular Weight": "mol_weight",
        }

        for prop in self.calculated_props:
            if field_name := prop_map.get(prop.kind):
                setattr(self, field_name, prop.value)

        for vid in self.external_ids:
            if vid.resource == "ChEMBL":
                self.chembl_id = vid.identifier
                break

        try:
            self.synonyms.remove(self.drugname)
        except ValueError:
            pass

        if self.inchi and not self.inchi_key:
            try:
                self.inchi_key = Chem.InchiToInchiKey(self.inchi)
                return self
            except Exception:
                pass

        if self.unii and (not self.inchi_key or not self.smiles):
            cal_data = unii_resolver.resolve(self.unii)
            if cal_data:
                self.inchi_key = self.inchi_key or cal_data.get("inchi_key")
                self.smiles = self.smiles or cal_data.get("smiles")
                self.formula = self.formula or cal_data.get("mol_formula")

        if self.smiles and not self.inchi_key:
            try:
                mol = Chem.MolFromSmiles(self.smiles)
                if mol:
                    self.inchi_key = Chem.MolToInchiKey(mol)
                    self.inchi = self.inchi or Chem.MolToInchi(mol)[0]
                    self.formula = self.formula or CalcMolFormula(mol)
                    self.mol_weight = self.mol_weight or CalcExactMolWt(mol)
            except Exception:
                pass
        return self


def _drug_processing_thread(raw_xml: bytes) -> dict | None:
    timer_start = time.perf_counter()

    try:
        import sys

        sys.setrecursionlimit(5000)

        root = etree.fromstring(raw_xml)

        drug_data = xml_to_dict(root)

        drugdata_model = DrugbankData(**drug_data)
        if not drugdata_model.inchi_key:
            clean_return = drugdata_model.model_dump(exclude_none=True)
            return {
                "pre_processed": True,
                "status": "skipped",
                "name": drugdata_model.drugname,
                "raw_data": clean_return,
                "time_taken": time.perf_counter() - timer_start,
            }

        exclude_keys = {
            "drugname",
            "inchi_key",
            "inchi",
            "smiles",
            "formula",
            "mol_weight",
            "drugbank_id",
            "chembl_id",
            "calculated_props",
            "experimental_props",
            "external_ids",
        }

        drug_cols = {
            "inchi_key": drugdata_model.inchi_key,
            "generic_name": drugdata_model.drugname,
            "inchi": drugdata_model.inchi,
            "smiles": drugdata_model.smiles,
            "chem_formula": drugdata_model.formula,
            "mol_weight": drugdata_model.mol_weight,
            "drugbank_id": drugdata_model.drugbank_id,
            "chembl_id": drugdata_model.chembl_id,
            "synonyms": drugdata_model.synonyms,
        }

        record_json = drugdata_model.model_dump(
            mode="json", exclude_none=True, exclude=exclude_keys
        )

        return {
            "pre_processed": True,
            "status": "success",
            "drug_cols": drug_cols,
            "record_json": record_json,
            "raw_backup": drug_data,
            "time_taken": time.perf_counter() - timer_start,
        }

    except Exception as e:
        return {
            "pre_processed": False,
            "error": str(e),
            "raw_data": {"bytes": len(raw_xml)},
            "time_taken": time.perf_counter() - timer_start,
        }


def xml_to_dict(element):
    if len(element) == 0 and not element.attrib:
        return element.text.strip() if element.text else ""

    result = {}
    if element.attrib:
        for k, v in element.attrib.items():
            attr_name = k.split("}", 1)[1] if "}" in k else k
            result[attr_name] = v

    for child in element:
        tag_name = child.tag
        if isinstance(tag_name, str) and "}" in tag_name:
            tag_name = tag_name.split("}", 1)[1]
        child_data = xml_to_dict(child)
        if tag_name in result:
            if isinstance(result[tag_name], list):
                result[tag_name].append(child_data)
            else:
                result[tag_name] = [result[tag_name], child_data]
        else:
            result[tag_name] = child_data

    text = element.text.strip() if element.text else ""
    if text:
        if not result:
            return text
        result["#text"] = text

    if len(result) == 1 and "#text" in result:
        return result["#text"]

    return result


class DrugbankScraper(BaseFileScraper):
    def __init__(
        self,
        session: Session,
        xml_file_path: str,
        xsd_file_path: str,
        version: str = None,
    ):
        if not os.path.exists(xml_file_path):
            raise FileNotFoundError(f"ERROR: source file not found: {xml_file_path}")

        if not os.path.exists(xsd_file_path):
            raise FileNotFoundError(f"ERROR: source file not found: {xsd_file_path}")

        if not version:
            try:
                content = etree.iterparse(xml_file_path, events=("start",))
                for _, element in content:
                    if element.tag.endswith("drugbank"):
                        version = element.get("version")
                        if not version:
                            raise ValueError("ERROR: <drugbank> tag missing")
                        break
            except etree.XMLSyntaxError as e:
                raise ValueError(f"Invalid XML file: {e}")

        self.schema_file_path = xsd_file_path
        super().__init__(session, "Drugbank", version, xml_file_path)

    def _extract(self) -> Iterator[dict]:
        ns = "{http://www.drugbank.ca}"

        context = etree.iterparse(
            self.file_path, events=("end",), tag=f"{ns}drug", huge_tree=True
        )

        xml_stats = {"count": 0, "total_time": 0.0}
        worker_stats = {"count": 0, "total_time": 0.0}

        def data_stream(ctx):
            t0 = time.perf_counter()
            for _, elem in ctx:
                if elem.get("type") == "small molecule":
                    xml_bytes = etree.tostring(elem)
                    t1 = time.perf_counter()
                    xml_stats["total_time"] += t1 - t0
                    xml_stats["count"] += 1
                    yield xml_bytes
                elem.clear()
                parent = elem.getparent()
                if parent is not None:
                    while elem.getprevious() is not None:
                        del parent[0]
                t0 = time.perf_counter()

        pool = Pool()
        try:
            print(f"\nNumber of worker processes: {pool._processes}")
            num_processes = int(pool._processes)

            results_stream = pool.imap_unordered(
                _drug_processing_thread, data_stream(context), chunksize=1
            )

            for result in results_stream:
                if result:
                    if "time_taken" in result:
                        worker_stats["count"] += 1
                        worker_stats["total_time"] += result["time_taken"]
                    yield result

        except GeneratorExit:
            pool.terminate()
            raise

        except Exception:
            pool.terminate()
            raise

        finally:
            pool.close()
            pool.join()
            print("\n" + "=" * 40)
            print("EXTRACTION PERF REPORT")
            print("=" * 40)
            if xml_stats["count"] > 0:
                avg_xml = xml_stats["total_time"] / xml_stats["count"]
                xml_rate = 1 / avg_xml if avg_xml > 0 else 0
                print("XML Parsing")
                print(f"  - Avg time per record: {avg_xml * 1000:.4f} ms")
                print(f"  - Max Throughput:      {xml_rate:.1f} records/sec")
            if worker_stats["count"] > 0:
                avg_worker = worker_stats["total_time"] / worker_stats["count"]
                print("\nWorker Validation")
                print(f"  - Avg time per record: {avg_worker * 1000:.4f} ms")
                print(
                    f"  - Theoretical Max:     {(1 / avg_worker) * num_processes:.1f} records/sec (across {num_processes} processses/cores)"
                )
            print("=" * 40 + "\n")
            del context

    def _element_to_dict(self, elem) -> dict:
        text = elem.text.strip() if elem.text else ""
        attribs = dict(elem.attrib)

        if len(elem) == 0:
            if attribs:
                attribs["#text"] = text
                return attribs
            return text

        result = {}

        if attribs:
            result.update(attribs)

        for child in elem:
            tag = child.tag.split("}")[-1] if "}" in child.tag else child.tag
            child_value = self._element_to_dict(child)

            if tag in result:
                if isinstance(result[tag], list):
                    result[tag].append(child_value)
                else:
                    result[tag] = [result[tag], child_value]
            else:
                result[tag] = child_value

        return result

    def _get_adapter(self) -> DrugbankAdapter:
        return DrugbankAdapter()


class DrugbankAdapter(BaseAdapter):
    def __init__(self):
        self.skipped = []

    def standardise(self, drug_data: Dict[str, any]) -> ScraperPayload | None:
        if "error" in drug_data:
            print(f"\n[WORKER ERROR] {drug_data.get('error')}")
            return None

        if drug_data.get("pre_processed"):
            status = drug_data.get("status")
            if status == "skipped":
                self.skipped.append(
                    {"name": drug_data["name"], "raw_data": drug_data["raw_data"]}
                )
                return None

            if status == "success":
                return ScraperPayload(
                    drug=DrugSchema(**drug_data["drug_cols"]),
                    record=RecordSchema(data_json=drug_data["record_json"]),
                )
        return None


if __name__ == "__main__":
    import sys
    from pathlib import Path

    from drug_db.database.manager import db

    BASE_DIR = Path(__file__).parent.parent.parent.resolve()
    data_path = BASE_DIR / ".." / "data" / "drugbank.xml"
    schema_path = BASE_DIR / ".." / "data" / "drugbank.xsd"

    db.create_tables()
    session = db.get_session()
    try:
        with session:
            scraper = DrugbankScraper(
                db.get_session(), data_path.resolve(), schema_path.resolve()
            )

            scraper.run()
    except KeyboardInterrupt:
        sys.exit(0)
    finally:
        session.close()
