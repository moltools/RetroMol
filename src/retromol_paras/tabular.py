"""Minimal reader for the tab-separated data files PARAS ships (vendored from parasect.core.tabular)."""

from collections import OrderedDict
from pathlib import Path


class Tabular:
    """A tab-separated file with a header row and a unique ID in the first column."""

    def __init__(self, path_in: str | Path, separator: str = "\t") -> None:
        path_in = Path(path_in)
        if not path_in.exists():
            raise FileNotFoundError(f"file not found: {path_in}")

        self.rows: "OrderedDict[str, OrderedDict[str, str]]" = OrderedDict()

        with open(path_in, "r") as fo:
            self.column_names = [n.strip() for n in fo.readline().strip().split(separator)]

            for line_idx, line in enumerate(fo):
                row = [v.strip() for v in line.strip().split(separator)]
                if len(row) != len(self.column_names):
                    raise ValueError(f"row {line_idx + 2} has a different number of columns than the header ({path_in})")

                row_id = row[0]
                if row_id in self.rows:
                    raise ValueError(f"duplicate row ID when reading {path_in}: {row_id}")

                self.rows[row_id] = OrderedDict(zip(self.column_names, row))

    def get_row_value(self, row_id: str, column_name: str) -> str:
        return self.rows[row_id][column_name]

    def get_row_values(self, row_id: str) -> list[str]:
        return list(self.rows[row_id].values())
