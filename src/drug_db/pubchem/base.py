from sqlalchemy import Column, Integer, String
from sqlalchemy.orm import declarative_base

PubChemBase = declarative_base()


class Identifier(PubChemBase):
    __tablename__ = "identifiers"

    cid = Column(Integer, primary_key=True, index=True)
    inchi = Column(String)
    inchi_key = Column(String, index=True)


class Structure(PubChemBase):
    __tablename__ = "structures"

    cid = Column(Integer, primary_key=True)
    smiles = Column(String)


class Synonym(PubChemBase):
    __tablename__ = "synonyms"

    id = Column(Integer, primary_key=True, autoincrement=True)
    cid = Column(Integer, index=True)
    synonym = Column(String, index=True)
