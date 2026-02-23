import collections
import sqlite3
from typing import Iterator

import pandas as pd
from sqlalchemy.orm import Session

from drug_db.models.schemas import DrugSchema, RecordSchema, ScraperPayload
from drug_db.scraper.adapters import BaseAdapter
from drug_db.scraper.base import BaseScraper


class ChemblAdapter(BaseAdapter):
    def standardise(self, raw_data: dict) -> ScraperPayload:
        drug = DrugSchema(
            inchi_key=raw_data.get("standard_inchi_key"),
            smiles=raw_data.get("canonical_smiles"),
            generic_name=raw_data.get("generic_name"),
        )

        record = RecordSchema(
            data_json={"chembl_pchembl_array": raw_data.get("pchembl_array")}
        )
        return ScraperPayload(drug=drug, record=record)


class ChemblScraper(BaseScraper):
    def __init__(self, session: Session, db_path: str, workers: int = 4):
        self.db_path, self.workers = db_path, workers
        super().__init__(session, source_name="ChEMBL", version="36.0")
        self.target_vocabulary = self._build_target_vocab()

    def _build_target_vocab(self, min_assay_count: int = 100) -> list:
        conn = sqlite3.connect(self.db_path)

        query = """
            SELECT td.pref_name as target_name, COUNT(*) as assay_count
            FROM activities a
            JOIN assays ass ON a.assay_id = ass.assay_id
            JOIN target_dictionary td ON ass.tid = td.tid
            WHERE a.pchembl_value IS NOT NULL
              AND td.organism = 'Homo sapiens'
              AND ass.assay_type IN ('B', 'F')
              AND ass.confidence_score >= 8
            GROUP BY td.pref_name
            HAVING COUNT(*) >= ?
            ORDER BY td.pref_name ASC
        """
        df = pd.read_sql_query(query, conn, params=(min_assay_count,))
        conn.close()

        dynamic_vocab = df["target_name"].tolist()
        print(
            f"Dynamically discovered {len(dynamic_vocab)} targets with at least {min_assay_count} assays."
        )

        return dynamic_vocab

    def _get_adapter(self) -> BaseAdapter:
        return ChemblAdapter()

    def _extract(self) -> Iterator[dict]:
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute("SELECT MAX(molregno) FROM activities")
        max_id = cursor.fetchone()[0] or 0
        conn.close()

        chunk_size = 2500

        escaped_vocab = [t.replace("'", "''") for t in self.target_vocabulary]
        targets_sql = "','".join(escaped_vocab)

        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row

        for start_id in range(1, max_id + 1, chunk_size):
            end_id = start_id + chunk_size - 1

            query = f"""
                SELECT 
                    a.molregno,
                    cs.standard_inchi_key, 
                    cs.canonical_smiles, 
                    md.pref_name, 
                    td.pref_name as target, 
                    a.pchembl_value as val
                FROM activities a
                JOIN molecule_dictionary md ON a.molregno = md.molregno
                LEFT JOIN compound_structures cs ON a.molregno = cs.molregno
                JOIN assays ass ON a.assay_id = ass.assay_id
                JOIN target_dictionary td ON ass.tid = td.tid
                WHERE a.molregno BETWEEN {start_id} AND {end_id}
                  AND td.organism = 'Homo sapiens' 
                  AND ass.assay_type IN ('B', 'F') 
                  AND td.pref_name IN ('{targets_sql}')
                  AND ass.confidence_score >= 8
            """

            molecule_data = collections.defaultdict(lambda: {"assays": {}})

            cursor = conn.execute(query)
            for row in cursor:
                m_id = row["molregno"]
                m_info = molecule_data[m_id]

                if "inchi" not in m_info:
                    m_info["inchi"] = row["standard_inchi_key"]
                    m_info["smiles"] = row["canonical_smiles"]
                    m_info["name"] = row["pref_name"]

                target = row["target"]
                raw_val = row["val"]
                val = float(raw_val) if raw_val is not None else 0.0

                if target not in m_info["assays"] or val > m_info["assays"][target]:
                    m_info["assays"][target] = val

            for d in molecule_data.values():
                if d["inchi"] and d["smiles"]:
                    yield {
                        "standard_inchi_key": d["inchi"],
                        "canonical_smiles": d["smiles"],
                        "generic_name": d["name"],
                        "pchembl_array": [
                            d["assays"].get(t, 0.0) for t in self.target_vocabulary
                        ],
                    }

            del molecule_data

        conn.close()


if __name__ == "__main__":
    import sys
    from pathlib import Path

    from drug_db.database.manager import db
    from drug_db.models.source import Source, SourceRecord

    BASE_DIR = Path(__file__).parent.parent.parent.resolve()
    DATA_DIR = BASE_DIR / ".." / "data"
    data_path = DATA_DIR / "chembl_36.db"

    db.create_tables()
    session = db.get_session()
    try:
        with session:
            print("Checking for existing ChEMBL data...")
            existing_source = (
                session.query(Source).filter(Source.name == "ChEMBL").first()
            )

            if existing_source:
                print(
                    "Found existing source. Purging old records to prevent duplication..."
                )
                session.query(SourceRecord).filter(
                    SourceRecord.source_id == existing_source.id
                ).delete()
                session.delete(existing_source)
                session.commit()
                print("Cleanup complete. Starting fresh!")

            scraper = ChemblScraper(session, str(data_path.resolve()))
            scraper.run()

    except KeyboardInterrupt:
        sys.exit(0)
    finally:
        session.close()
