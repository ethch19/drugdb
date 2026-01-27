import multiprocessing
from pathlib import Path

import pandas as pd
from sqlalchemy import create_engine, inspect, text
from tqdm import tqdm

from drug_db.pubchem.base import Identifier, Structure, Synonym

CUR_DIR = Path(__file__).parent
ROOT_DIR = CUR_DIR.parent.parent.parent
DATA_DIR = ROOT_DIR / "data"
DATASTORE_DIR = ROOT_DIR / "datastore"

DATA_DIR.mkdir(parents=True, exist_ok=True)
DATASTORE_DIR.mkdir(parents=True, exist_ok=True)

TASKS = [
    {
        "id": 0,
        "source": "CID-InChI-Key.gz",
        "db": "pubchem_identifiers.db",
        "model": Identifier,
        "cols": ["cid", "inchi", "inchi_key"],
    },
    {
        "id": 1,
        "source": "CID-SMILES.gz",
        "db": "pubchem_structures.db",
        "model": Structure,
        "cols": ["cid", "smiles"],
    },
    {
        "id": 2,
        "source": "CID-Synonym-filtered.gz",
        "db": "pubchem_synonyms.db",
        "model": Synonym,
        "cols": ["cid", "synonym"],
    },
]


def get_existing_count(engine, table_name):
    try:
        insp = inspect(engine)
        if not insp.has_table(table_name):
            return 0

        with engine.connect() as conn:
            result = conn.execute(text(f"SELECT COUNT(*) FROM {table_name}"))
            return result.scalar()
    except Exception:
        return 0


def process_task(task):
    print(f"Processing {task['db']}...")
    db_path = DATASTORE_DIR / f"{task['db']}"
    engine = create_engine(f"sqlite:///{db_path.absolute()}", echo=False)
    table_name = task["model"].__tablename__

    task["model"].metadata.create_all(engine)
    existing_rows = get_existing_count(engine, table_name)

    start_msg = f"Starting {task['db']}"
    if existing_rows > 0:
        start_msg = f"Resuming {task['db']} from row {existing_rows:,}"

    progress_bar = tqdm(
        desc=f"Building {task['db']}", position=task["id"], unit="rows", leave=True
    )
    progress_bar.write(start_msg)

    chunk_size = 100_000

    source_file = DATA_DIR / f"{task['source']}"

    try:
        with pd.read_csv(
            source_file.absolute(),
            sep="\t",
            names=task["cols"],
            chunksize=chunk_size,
            compression="gzip",
            quoting=3,
            dtype=str,
            skiprows=existing_rows,
        ) as reader:
            for chunk in reader:
                if "cid" in chunk.columns:
                    chunk["cid"] = pd.to_numeric(chunk["cid"], errors="coerce")
                chunk.to_sql(
                    task["model"].__tablename__, engine, if_exists="append", index=False
                )
                progress_bar.update(len(chunk))
    except Exception as e:
        progress_bar.write(f"ERROR: in {task['db']}: {e}")
    finally:
        progress_bar.close()


if __name__ == "__main__":
    try:
        with multiprocessing.Pool(len(TASKS)) as p:
            result = p.map_async(process_task, TASKS)
            result.get(timeout=999999)
        print("\nAll PubChem databases built successfully")
    except KeyboardInterrupt:
        pass
