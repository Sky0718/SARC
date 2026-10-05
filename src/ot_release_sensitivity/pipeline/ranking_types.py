from __future__ import annotations
from typing import Any
from typing import Sequence

class RankingSensitivityBuildError(RuntimeError):
    pass

class DuckDbDataset:
    def __init__(
        self, table: str, fields: Sequence[str], types: Sequence[str], batch_size: int
    ):
        self.table = table
        self.fields = tuple(fields)
        self.types = tuple(types)
        self.batch_size = batch_size
        self.row_count = 0
        self.persisted_count = 0
        self.pending_rows: list[tuple[Any, ...]] = []
        self.flush_count = 0
        self.flushed_batch_sizes: list[int] = []

    @property
    def pending_count(self) -> int:
        return len(self.pending_rows)
