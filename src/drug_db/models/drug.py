from types import NotImplementedType
from typing import TYPE_CHECKING

from sqlalchemy import JSON
from sqlalchemy.orm import (
    Mapped,
    mapped_column,
    relationship,
)

from .base import Base

if TYPE_CHECKING:
    from .source import SourceRecord


class Drug(Base):  # normalised (core data model)
    __tablename__ = "drugs"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True, init=False)
    inchi_key: Mapped[str] = mapped_column(
        unique=True, index=True
    )  # derived from IUPAC
    inchi: Mapped[str]
    chem_formula: Mapped[str]
    mol_weight: Mapped[float]
    generic_name: Mapped[str | None] = mapped_column(default=None)  # INN
    smiles: Mapped[str | None] = mapped_column(default=None)  # canonical
    drugbank_id: Mapped[str | None] = mapped_column(default=None)
    chembl_id: Mapped[str | None] = mapped_column(default=None)
    synonyms: Mapped[list[str]] = mapped_column(
        JSON, default=list
    )  # USAN, brand names, regional names, etc

    sources: Mapped[list["SourceRecord"]] = relationship(
        back_populates="drug", cascade="all, delete-orphan", default_factory=list
    )

    def __str__(self) -> str:
        return f"""
        Drug Class Obj
        --------------- 
        Id: {self.id}
        Generic Name: {self.generic_name}
        Synonyms: {", ".join(self.synonyms)}
        InChIKey: {self.inchi_key}
        InChI: {self.inchi}
        Canonical SMILES: {self.smiles}
        Chemical Formula: {self.chem_formula}
        Molecular Weight: {round(self.mol_weight, 3)}
        DrugBank Id: {self.drugbank_id}
        ChEMBL Id: {self.chembl_id}
        """

    def __eq__(self, other: object) -> bool | NotImplementedType:
        if not isinstance(other, Drug):
            return NotImplemented
        return self.inchi_key == other.inchi_key

    def __hash__(self) -> int:
        return hash(self.inchi_key)
