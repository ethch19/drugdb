from abc import ABC, abstractmethod

from drug_db.scraper.base import ScraperPayload


class BaseAdapter(ABC):
    @abstractmethod
    def standardise(self, raw_data: any) -> ScraperPayload | None:
        pass
