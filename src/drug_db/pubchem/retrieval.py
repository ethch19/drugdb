from pathlib import Path

from sqlalchemy import (
    Column,
    Integer,
    MetaData,
    String,
    Table,
    create_engine,
    select,
    text,
)

CUR_DIR = Path(__file__).parent
ROOT_DIR = CUR_DIR.parent.parent.parent
DATASTORE_DIR = ROOT_DIR / "datastore"


class PubChemRetriever:
    def __init__(self):
        self.engine = create_engine(
            f"sqlite:///{DATASTORE_DIR / 'pubchem_identifiers.db'}", echo=False
        )
        self.conn = self.engine.connect()
        self.conn.execute(
            text(
                f"ATTACH DATABASE '{(DATASTORE_DIR / 'pubchem_structures.db').absolute()}' AS struct_db"
            )
        )
        self.conn.execute(
            text(
                f"ATTACH DATABASE '{(DATASTORE_DIR / 'pubchem_synonyms.db').absolute()}' AS syn_db"
            )
        )

        self.metadata = MetaData()
        self.t_identifiers = Table(
            "identifiers",
            self.metadata,
            Column("cid", Integer, primary_key=True),
            Column("inchi", String),
            Column("inchi_key", String),
        )
        self.t_structures = Table(
            "structures",
            self.metadata,
            Column("cid", Integer),
            Column("smiles", String),
            schema="struct_db",
        )
        self.t_synonyms = Table(
            "synonyms",
            self.metadata,
            Column("cid", Integer),
            Column("synonym", String),
            schema="syn_db",
        )

    def get_by_name(self, name: str):
        query = (
            select(self.t_synonyms.c.cid)
            .where(self.t_synonyms.c.synonym == name)
            .limit(1)
        )

        result = self.conn.execute(query).fetchone()

        if result:
            return self.get_by_cid(result.cid)
        return None

    def get_by_inchi_key(self, inchi_key: str):
        query = (
            select(
                self.t_identifiers.c.cid,
                self.t_identifiers.c.inchi,
                self.t_identifiers.c.inchi_key,
                self.t_structures.c.smiles,
                self.t_synonyms.c.synonym,
            )
            .select_from(
                self.t_identifiers.join(
                    self.t_structures,
                    self.t_identifiers.c.cid == self.t_structures.c.cid,
                    isouter=True,
                ).join(
                    self.t_synonyms,
                    self.t_identifiers.c.cid == self.t_synonyms.c.cid,
                    isouter=True,
                )
            )
            .where(self.t_identifiers.c.inchi_key == inchi_key)
            .limit(1)
        )

        row = self.conn.execute(query).mappings().fetchone()
        return dict(row) if row else None

    def get_by_cid(self, cid: int):
        query = (
            select(
                self.t_identifiers.c.cid,
                self.t_identifiers.c.inchi,
                self.t_identifiers.c.inchi_key,
                self.t_structures.c.smiles,
                self.t_synonyms.c.synonym,
            )
            .select_from(
                self.t_identifiers.join(
                    self.t_structures,
                    self.t_identifiers.c.cid == self.t_structures.c.cid,
                    isouter=True,
                ).join(
                    self.t_synonyms,
                    self.t_identifiers.c.cid == self.t_synonyms.c.cid,
                    isouter=True,
                )
            )
            .where(self.t_identifiers.c.cid == cid)
            .limit(1)
        )

        row = self.conn.execute(query).mappings().fetchone()
        return dict(row) if row else None

    def close(self):
        self.conn.close()
