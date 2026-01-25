from pathlib import Path

from sqlalchemy import Engine, create_engine, event
from sqlalchemy.orm import Session, sessionmaker

from drug_db.models.base import Base

CUR_DIR = Path(__file__).parent
ROOT_DIR = CUR_DIR.parent.parent.parent


@event.listens_for(Engine, "connect")
def set_sqlite_pragma(dbapi_connection, connection_record):
    ac = dbapi_connection.autocommit
    dbapi_connection.autocommit = True

    cursor = dbapi_connection.cursor()
    cursor.execute("PRAGMA foreign_keys=ON")
    cursor.close()

    dbapi_connection.autocommit = ac


class DbManager:
    def __init__(self, db_name: str = "drugs"):
        self.db_path = ROOT_DIR / "data" / f"{db_name}.db"
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self.engine = create_engine(f"sqlite:////{self.db_path.absolute()}", echo=False)
        self.session = sessionmaker(bind=self.engine)

    def create_tables(self):
        Base.metadata.create_all(self.engine)

    def get_session(self) -> Session:
        return self.session()


db = DbManager()
