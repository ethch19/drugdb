from types import NotImplementedType
from typing import TYPE_CHECKING

from sqlalchemy import JSON, ForeignKey, UniqueConstraint
from sqlalchemy.orm import Mapped, mapped_column, relationship

from .base import Base

if TYPE_CHECKING:
    from .drug import Drug


class Source(Base):
    __tablename__ = "sources"
    __table_args__ = (
        UniqueConstraint("name", "version", name="unique_source_version"),
    )

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True, init=False)
    name: Mapped[str]
    version: Mapped[str]
    date_accessed: Mapped[str]  # ISO8601 "YYYY-MM-DD HH:MM:SS.SSS"
    url: Mapped[str | None] = mapped_column(default=None)

    records: Mapped[list["SourceRecord"]] = relationship(
        back_populates="source", cascade="all, delete-orphan", default_factory=list
    )

    def __str__(self) -> str:
        return f"""
        Source Class Obj
        --------------- 
        id: {self.id}
        name: {self.name}
        version: {self.version}
        date accessed: {self.date_accessed}
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

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True, init=False)
    drug_inchi_key: Mapped[str] = mapped_column(
        ForeignKey("drugs.inchi_key", ondelete="CASCADE"), index=True
    )
    source_id: Mapped[int] = mapped_column(ForeignKey("sources.id"), index=True)
    data_json: Mapped[dict] = mapped_column(
        JSON
    )  # for changes detection, must be fully reassigned

    drug: Mapped["Drug"] = relationship(back_populates="sources", init=False)
    source: Mapped["Source"] = relationship(back_populates="records", init=False)

    def __str__(self) -> str:
        return f"""
        SourceRecord Class Obj
        --------------- 
        Id: {self.id}
        Source Id: {self.source_id}
        InChI Key: {self.drug_inchi_key}
        """
        # Data Json: {json.dumps(self.data_json, indent=2)}

    def __eq__(self, other: object) -> bool | NotImplementedType:
        if not isinstance(other, SourceRecord):
            return NotImplemented
        return self.drug_inchi_key == other.drug_inchi_key
