"""Force a headless matplotlib backend before pyplot is imported anywhere.

Importing this module (instead of matplotlib.pyplot directly) guarantees the Agg
backend, so teloviz renders on login/compute nodes with no display.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402  (must follow use())

__all__ = ["plt", "matplotlib"]
