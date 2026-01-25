import os
from abc import ABC, abstractmethod
from datetime import datetime, timezone
from typing import Iterator

import requests
from sqlalchemy.orm import Session

from drug_db.models.drug import Drug
from drug_db.models.schemas import ScraperPayload
from drug_db.models.source import Source, SourceRecord

from .adapters import BaseAdapter


class BaseScraper(ABC):
    def __init__(self, session: Session, source_name: str, version: str):
        self.session = session
        self.source_name = source_name
        self.version = version

        self.source_id = self._register_source()
        self.adapter = self._get_adapter()

    @abstractmethod
    def _extract(self) -> Iterator[any]:
        pass

    @abstractmethod
    def _get_adapter(self) -> BaseAdapter:
        pass

    def run(self):
        data = self._extract()
        count = 0
        for x in data:
            clean_x = self._standardise(x)
            if clean_x:
                self._load(clean_x)
                count += 1

                if count % 1000 == 0:
                    self.session.commit()
        self.session.commit()

    def _load(self, payload: ScraperPayload):
        drug_data = payload.drug.model_dump()

        init_drug = self.session.query(Drug).get(payload.drug.inchi_key)

        if init_drug:
            final_drug = init_drug
        else:
            final_drug = Drug(**drug_data)
            self.session.add(final_drug)

        source_record = SourceRecord(
            source_id=self.source_id,
            drug_inchi_key=final_drug.inchi_key,
            data_json=payload.record.data_json,
        )
        self.session.add(source_record)

    def _standardise(self, raw_data: any) -> any:
        try:
            return self.adapter.standardise(raw_data)
        except Exception as e:
            print(f"ERROR: Cannot standardise data from adapter: {e}")
            return None

    def _register_source(self) -> int:
        src = (
            self.session.query(Source)
            .filter_by(name=self.source_name, version=self.version)
            .first()
        )
        if not src:
            current_date = datetime.now(timezone.utc).date()
            new_src = Source(
                name=self.source_name,
                version=self.version,
                url=None,
                date_accessed=current_date.isoformat(),
                records=[],
            )
            self.session.add(new_src)
            self.session.commit()
            return new_src.id
        return src.id


class BaseFileScraper(BaseScraper, ABC):
    def __init__(
        self, session: Session, source_name: str, version: str, file_path: str
    ):
        super().__init__(session, source_name, version)
        self.file_path = file_path

        if not os.path.exists(file_path):
            raise FileNotFoundError(f"ERROR: source file not found: {file_path}")


class BaseApiScraper(BaseScraper, ABC):
    def __init__(self, session: Session, source_name: str, version: str, base_url: str):
        super().__init__(session, source_name, version)
        self.base_url = base_url
        self.http_sesh = requests.Session()

    def _register_source(self) -> int:
        src = (
            self.session.query(Source)
            .filter_by(name=self.source_name, version=self.version)
            .first()
        )
        if not src:
            current_date = datetime.now(timezone.utc).date()
            new_src = Source(
                name=self.source_name,
                version=self.version,
                url=self.base_url,
                date_accessed=current_date.isoformat(),
                records=[],
            )
            self.session.add(new_src)
            self.session.commit()
            return new_src.id
        return src.id


class BaseWebScraper(BaseScraper, ABC):
    def __init__(self, session: Session, source_name: str, version: str, base_url: str):
        super().__init__(session, source_name, version)
        self.base_url = base_url
        # selenium implem, not yet installed

    def _register_source(self) -> int:
        src = (
            self.session.query(Source)
            .filter_by(name=self.source_name, version=self.version)
            .first()
        )
        if not src:
            current_date = datetime.now(timezone.utc).date()
            new_src = Source(
                name=self.source_name,
                version=self.version,
                url=self.base_url,
                date_accessed=current_date.isoformat(),
                records=[],
            )
            self.session.add(new_src)
            self.session.commit()
            return new_src.id
        return src.id
