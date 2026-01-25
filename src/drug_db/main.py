from drug_db.database.manager import db
from drug_db.models.base import Base


def main():
    Base.metadata.drop_all(db.engine)
    db.create_tables()


if __name__ == "__main__":
    main()
