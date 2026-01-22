from types import NotImplementedType
from typing import TYPE_CHECKING

from sqlalchemy import ForeignKey, UniqueConstraint
from sqlalchemy.orm import Mapped, mapped_column, relationship

from .base import Base

if TYPE_CHECKING:
    from .drug import Drug


class Source(Base):
    __tablename__ = "sources"

    id: Mapped[int] = mapped_column(primary_key=True)
    name: Mapped[str] = mapped_column(unique=True)
    url: Mapped[str | None]

    records: Mapped[list["SourceRecord"]] = relationship(
        back_populates="source", cascade="all, delete-orphan"
    )

    def __str__(self) -> str:
        return f"""
        Source Class Obj
        --------------- 
        id: {self.id}
        name: {self.name}
        url: {self.url}
        """

    def __eq__(self, other: object) -> bool | NotImplementedType:
        if not isinstance(other, Source):
            return NotImplemented
        return self.id == other.id

    def __hash__(self) -> int:
        return hash(self.name)


class SourceRecord(Base):
    __tablename__ = "source_records"
    __table_args__ = (
        UniqueConstraint("source_id", "drug_inchi_key", name="unique_source_drug"),
    )

    id: Mapped[int] = mapped_column(primary_key=True)
    drug_inchi_key: Mapped[str] = mapped_column(
        ForeignKey("drugs.inchi_key", ondelete="CASCADE")
    )
    source_id: Mapped[int] = mapped_column(ForeignKey("sources.id"))
    data_json: Mapped[str]

    drug: Mapped["Drug"] = relationship(back_populates="sources")
    source: Mapped["Source"] = relationship(back_populates="records")

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
