import json
import random
import time
from datetime import datetime
from pathlib import Path

from sqlalchemy import and_, func, not_, or_, select
from sqlalchemy.orm import Session
from tqdm import tqdm

from drug_db.database.manager import db
from drug_db.inferences.protox import ProToxBuilder, ProToxModel
from drug_db.models.drug import Drug
from drug_db.models.source import Source, SourceRecord

CUR_DIR = Path(__file__).parent.resolve()
ROOT_DIR = CUR_DIR.parent.parent.parent
DATASTORE_DIR = ROOT_DIR / "datastore"
OUTPUT_FILE = DATASTORE_DIR / "withdrawn_protox_results.jsonl"


def get_withdrawn_smiles(session: Session):
    drugbank_criteria = and_(
        Source.name.ilike("%DrugBank%"),
        func.json_extract(SourceRecord.data_json, "$.groups").like('%"withdrawn"%'),
    )

    withdrawn_db_criteria = and_(
        Source.name.ilike("%Withdrawn DB%"),
        or_(
            func.json_extract(SourceRecord.data_json, "$.first_withdrawn_year").isnot(
                None
            ),
            and_(
                func.json_extract(
                    SourceRecord.data_json, "$.withdrawal_countries"
                ).isnot(None),
                func.json_extract(SourceRecord.data_json, "$.withdrawal_countries")
                != "[]",
            ),
        ),
    )

    stmt = (
        select(Drug.smiles, Drug.generic_name, Drug.inchi_key)
        .join(SourceRecord, Drug.inchi_key == SourceRecord.drug_inchi_key)
        .join(Source, SourceRecord.source_id == Source.id)
        .where(or_(drugbank_criteria, withdrawn_db_criteria))
        .distinct()
    )

    results = session.execute(stmt).all()

    clean_results = [
        {"smiles": row.smiles, "name": row.generic_name, "id": row.inchi_key}
        for row in results
        if row.smiles
    ]

    print(f"Found {len(clean_results)} withdrawn drugs with SMILES.")
    return clean_results


def get_investigation_smiles(session):
    print(
        "Querying database for 'Approved' drugs (excluding experimental/withdrawn)..."
    )

    groups_json = func.json_extract(SourceRecord.data_json, "$.groups")

    criteria = and_(
        Source.name.ilike("%DrugBank%"),
        groups_json.like('%"approved"%'),  # MUST have approved
        not_(groups_json.like('%"experimental"%')),  # MUST NOT have experimental
        not_(groups_json.like('%"withdrawn"%')),  # MUST NOT have withdrawn
    )

    stmt = (
        select(Drug.smiles, Drug.generic_name, Drug.inchi_key)
        .join(SourceRecord, Drug.inchi_key == SourceRecord.drug_inchi_key)
        .join(Source, SourceRecord.source_id == Source.id)
        .where(criteria)
        .distinct()
    )

    results = session.execute(stmt).all()
    clean_results = [
        {"smiles": row.smiles, "name": row.generic_name, "id": row.inchi_key}
        for row in results
        if row.smiles
    ]
    print(f"Found {len(clean_results)} drugs meeting criteria.")
    return clean_results


def run_batch_inference():
    DATASTORE_DIR.mkdir(parents=True, exist_ok=True)

    session = db.get_session()
    model = ProToxModel()

    config = ProToxBuilder().with_ld50().with_tox_class().with_pred_acc().build()

    try:
        drugs = get_withdrawn_smiles(session)
    finally:
        session.close()

    processed_ids = set()
    if OUTPUT_FILE.exists():
        with open(OUTPUT_FILE, "r") as f:
            for line in f:
                try:
                    processed_ids.add(json.loads(line)["id"])
                except:
                    pass
        print(f"Resuming... {len(processed_ids)} drugs already processed.")

    print(f"Saving results to: {OUTPUT_FILE}")

    with open(OUTPUT_FILE, "a") as f_out:
        for drug in tqdm(drugs, desc="Running ProTox Inference"):
            if drug["id"] in processed_ids:
                continue

            try:
                time.sleep(random.uniform(5, 10))

                result = model.execute(drug["smiles"], config)

                if result:
                    output_record = {
                        "id": drug["id"],
                        "name": drug["name"],
                        "smiles": drug["smiles"],
                        "timestamp": datetime.now().isoformat(),
                        "inference": result,
                    }
                    f_out.write(json.dumps(output_record) + "\n")
                    f_out.flush()
                else:
                    error_rec = {"id": drug["id"], "error": "No result returned"}
                    f_out.write(json.dumps(error_rec) + "\n")
                    f_out.flush()

            except Exception as e:
                print(f"\nError processing {drug['name']}: {e}")


def get_average_tox_class(file_path: str | Path) -> float:
    total_score = 0
    count = 0

    path = Path(file_path)
    if not path.exists():
        print(f"File not found: {path}")
        return 0.0

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            try:
                if not line.strip():
                    continue

                record = json.loads(line)

                inference = record.get("inference")
                if not isinstance(inference, dict):
                    continue

                tox_class_str = inference.get("tox_class")

                if tox_class_str and tox_class_str.isdigit():
                    total_score += int(tox_class_str)
                    count += 1

            except (json.JSONDecodeError, ValueError):
                continue

    if count == 0:
        return 0.0

    average = total_score / count
    print(f"Processed {count} records. Average Tox Class: {average:.2f}")
    return average


if __name__ == "__main__":
    # run_batch_inference()
    avg = get_average_tox_class(OUTPUT_FILE.resolve())
