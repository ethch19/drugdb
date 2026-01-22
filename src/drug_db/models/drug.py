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


class Drug(Base):
    __tablename__ = "drugs"

    generic_name: Mapped[str]  # INN/USAN
    iupac_name: Mapped[str]  # PubChem
    brand_names: Mapped[list[str] | None] = mapped_column(JSON)
    inchi: Mapped[str]
    inchi_key: Mapped[str] = mapped_column(primary_key=True)
    smiles: Mapped[str | None]  # canonical
    chem_formula: Mapped[str]
    mol_weight: Mapped[float]
    drugbank_id: Mapped[str | None]

    sources: Mapped[list["SourceRecord"]] = relationship(
        back_populates="drug", cascade="all, delete-orphan"
    )

    def __str__(self) -> str:
        return f"""
        Drug Class Obj
        --------------- 
        Generic Name: {self.generic_name}
        IUPAC Name: {self.iupac_name}
        Brand Names: {", ".join(self.brand_names) if self.brand_names else "None"}
        InChI: {self.inchi}
        SMILES: {self.smiles}
        Chemical Formula: {self.chem_formula}
        Molecular Weight: {round(self.mol_weight, 3)}
        """

    def __eq__(self, other: object) -> bool | NotImplementedType:
        if not isinstance(other, Drug):
            return NotImplemented
        return self.inchi_key == other.inchi_key

    def __hash__(self) -> int:
        return hash(self.inchi_key)
