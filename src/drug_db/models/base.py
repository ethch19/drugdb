from sqlalchemy.orm import (
    DeclarativeBase,
    MappedAsDataclass,
)


class Base(MappedAsDataclass, DeclarativeBase):  # shared registry
    pass
