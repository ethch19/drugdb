import json
import os
import time
from abc import ABC, abstractmethod
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator

import requests
from sqlalchemy.orm import Session
from tqdm import tqdm

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
        timers = {
            "extraction": 0.0,
            "standardisation": 0.0,
            "db_load": 0.0,
        }
        count = 0
        skipped_count = 0
        data = self._extract()
        iterator = iter(tqdm(data, desc="Processing database", unit="record"))

        while True:
            try:
                t0 = time.perf_counter()
                x = next(iterator)
                timers["extraction"] += time.perf_counter() - t0

                t1 = time.perf_counter()
                clean_x = self._standardise(x)
                timers["standardisation"] += time.perf_counter() - t1

                if clean_x:
                    t2 = time.perf_counter()
                    self._load(clean_x)
                    timers["db_load"] += time.perf_counter() - t2
                    count += 1

                    if count % 1000 == 0:
                        t3 = time.perf_counter()
                        self.session.commit()
                        timers["db_load"] += time.perf_counter() - t3

                    if (count + skipped_count) > 1000:
                        break
                else:
                    skipped_count += 1

            except StopIteration:
                break

        self.session.commit()

        total_time = (
            timers["extraction"] + timers["standardisation"] + timers["db_load"]
        )
        total_processed = count + skipped_count
        print("\n" + "=" * 40)
        print("PERFORMANCE REPORT")
        print("=" * 40)
        print(f"Total Time:     {total_time:.2f}s")
        print(f"Processed:      {total_processed} records")
        print(f"  - Saved:      {count}")
        print(f"  - Skipped:    {skipped_count}")
        print(f"Speed:      {count / total_time:.1f} saved/sec")
        print("=" * 40)
        print(
            f"1. Extraction:      {timers['extraction']:.2f}s ({(timers['extraction'] / total_time) * 100:.1f}%)"
        )
        print(
            f"2. Standardisation: {timers['standardisation']:.2f}s ({(timers['standardisation'] / total_time) * 100:.1f}%)"
        )
        print(
            f"3. DB load:         {timers['db_load']:.2f}s ({(timers['db_load'] / total_time) * 100:.1f}%)"
        )

        print("=" * 40)
        if hasattr(self.adapter, "skipped") and self.adapter.skipped:
            skipped_count = len(self.adapter.skipped)
            print(f"Saving {skipped_count} records that are skipped")
            repo_root = Path(__file__).resolve().parents[3]
            skipped_json = (
                repo_root
                / "data"
                / (self.source_name.strip().lower() + "_skipped.json")
            )
            with open(skipped_json, "w", encoding="utf-8") as f:
                json.dump(self.adapter.skipped, f, indent=2)

    def _load(self, payload: ScraperPayload):
        drug_data = payload.drug.model_dump()

        init_drug = (
            self.session.query(Drug)
            .filter(Drug.inchi_key == payload.drug.inchi_key)
            .first()
        )

        if init_drug:
            final_drug = init_drug
        else:
            final_drug = Drug(**drug_data)
            self.session.add(final_drug)

        existing_record = (
            self.session.query(SourceRecord)
            .filter_by(source_id=self.source_id, drug_inchi_key=final_drug.inchi_key)
            .first()
        )

        new_data = payload.record.data_json
        if existing_record:
            print(
                f"\nDuplicate {payload.drug.inchi_key} ({drug_data.get('generic_name')})"
            )

            current_record = existing_record.data_json or {}

            for key, new_val in new_data.items():
                if key in current_record:
                    cur_val = current_record[key]
                    if isinstance(cur_val, list) and isinstance(new_val, list):
                        try:
                            combined = list(set(cur_val + new_val))
                            current_record[key] = combined
                        except TypeError:
                            current_record[key].extend(new_val)
                    elif isinstance(cur_val, dict) and isinstance(new_val, dict):
                        cur_val.update(new_val)
                        current_record[key] = cur_val
                else:
                    current_record[key] = new_val

            existing_record.data_json = dict(current_record)
        else:
            source_record = SourceRecord(
                source_id=self.source_id,
                drug_inchi_key=final_drug.inchi_key,
                data_json=new_data,
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
