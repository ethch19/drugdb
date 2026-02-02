from pathlib import Path

from sqlalchemy import (
    Column,
    Integer,
    MetaData,
    NullPool,
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
            f"sqlite:///{DATASTORE_DIR / 'pubchem_identifiers.db'}",
            echo=False,
            poolclass=NullPool,
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

    def _get_base_query(self, include_synonyms: bool):
        col = [
            self.t_identifiers.c.cid,
            self.t_identifiers.c.inchi,
            self.t_identifiers.c.inchi_key,
            self.t_structures.c.smiles,
        ]

        j = self.t_identifiers.join(
            self.t_structures,
            self.t_identifiers.c.cid == self.t_structures.c.cid,
            isouter=True,
        )

        if include_synonyms:
            col.append(self.t_synonyms.c.synonym)
            j = j.join(
                self.t_synonyms,
                self.t_identifiers.c.cid == self.t_synonyms.c.cid,
                isouter=True,
            )

        return select(*col).select_from(j)

    def _aggregate_rows(self, rows):
        if not rows:
            return None

        first_row = rows[0]
        result = dict(first_row)
        if "synonym" not in result:
            return result

        result["synonyms"] = []
        result.pop("synonym", None)

        seen_synonyms = set()
        for row in rows:
            syn = row["synonym"]
            if syn and syn not in seen_synonyms:
                result["synonyms"].append(syn)
                seen_synonyms.add(syn)
        return result

    def get_by_name(self, name: str, include_synonyms: bool = False):
        subquery = (
            select(self.t_synonyms.c.cid)
            .where(self.t_synonyms.c.synonym.ilike(name))
            .limit(1)
            .scalar_subquery()
        )

        query = self._get_base_query(include_synonyms).where(
            self.t_identifiers.c.cid == subquery
        )

        rows = self.conn.execute(query).mappings().fetchall()
        return self._aggregate_rows(rows)

    def get_by_inchi_key(self, inchi_key: str, include_synonyms: bool = False):
        query = self._get_base_query(include_synonyms).where(
            self.t_identifiers.c.inchi_key == inchi_key
        )
        rows = self.conn.execute(query).mappings().fetchall()
        return self._aggregate_rows(rows)

    def get_by_cid(self, cid: int, include_synonyms: bool = False):
        query = self._get_base_query(include_synonyms).where(
            self.t_identifiers.c.cid == cid
        )
        rows = self.conn.execute(query).mappings().fetchall()
        return self._aggregate_rows(rows)

    def close(self):
        self.conn.close()


pubchem_retriever = PubChemRetriever()

# if __name__ == "__main__":
#     result = pubchem_retriever.get_by_name("aspirin")
#     print(result)
