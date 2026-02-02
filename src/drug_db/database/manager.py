import os
from pathlib import Path

from sqlalchemy import Engine, create_engine, event, text
from sqlalchemy.orm import Session, sessionmaker

from drug_db.models.base import Base

CUR_DIR = Path(__file__).parent
ROOT_DIR = CUR_DIR.parent.parent.parent
DATASTORE_DIR = ROOT_DIR / "datastore"


@event.listens_for(Engine, "connect")
def set_sqlite_pragma(dbapi_connection, connection_record):
    ac = dbapi_connection.autocommit
    dbapi_connection.autocommit = True

    cursor = dbapi_connection.cursor()
    cursor.execute("PRAGMA foreign_keys=ON")
    mode = os.environ.get("DB_MODE", "WAL")
    cursor.execute(f"PRAGMA journal_mode={mode}")
    cursor.execute("PRAGMA busy_timeout=5000")
    cursor.close()

    dbapi_connection.autocommit = ac


class DbManager:
    def __init__(self, db_name: str = "drugs"):
        self.db_path = DATASTORE_DIR / f"{db_name}.db"
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self.engine = create_engine(f"sqlite:////{self.db_path.absolute()}", echo=False)
        self.session = sessionmaker(bind=self.engine)

    def create_tables(self):
        Base.metadata.create_all(self.engine)

    def get_session(self) -> Session:
        return self.session()


db = DbManager()

if __name__ == "__main__":
    os.environ["DB_MODE"] = "DELETE"

    try:
        with db.engine.connect() as conn:
            print("   > Merging WAL data...")
            conn.execute(text("VACUUM;"))

        print("Files merged.")
    except Exception as e:
        print(f"Error: {e}")
