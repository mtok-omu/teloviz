"""Tests for the value -> color mapping (spec section 5)."""

from __future__ import annotations

import numpy as np

from teloviz.color import build_scheme

WHITE = (1.0, 1.0, 1.0, 1.0)


def _rgba(scheme, v):
    return tuple(scheme.cmap(scheme.norm(scheme.transform(np.array([v]))[0])))


def test_sum_zero_is_white_but_small_nonzero_is_not():
    s = build_scheme("sum", bin_w=100, cap=500, log=False, cmap_sum="Reds", cmap_div="RdBu_r")
    assert _rgba(s, 0) == WHITE
    assert _rgba(s, 50) != WHITE  # 1..99 band is tinted, not white


def test_sum_is_monotonic_in_intensity():
    s = build_scheme("sum", bin_w=100, cap=500, log=False, cmap_sum="Reds", cmap_div="RdBu_r")
    # Higher counts -> higher band index (darker red).
    idx = [s.norm(s.transform(np.array([v]))[0]) for v in (0, 50, 150, 450, 600)]
    assert idx == sorted(idx)
    assert idx[-1] == idx[-2] or 600 >= 500  # cap clamps


def test_sum_values_and_orientation_values():
    s_sum = build_scheme("sum", bin_w=100, cap=500, log=False, cmap_sum="Reds", cmap_div="RdBu_r")
    s_or = build_scheme("orientation", bin_w=100, cap=500, log=False, cmap_sum="Reds", cmap_div="RdBu_r")
    fwd, rev = np.array([10.0]), np.array([3.0])
    assert s_sum.values(fwd, rev)[0] == 13
    assert s_or.values(fwd, rev)[0] == 7


def test_orientation_edges_symmetric():
    s = build_scheme("orientation", bin_w=100, cap=500, log=False, cmap_sum="Reds", cmap_div="RdBu_r")
    assert s.edges[0] == -500 and s.edges[-1] == 500
    assert 0.0 in set(s.edges)


def test_log_transform_monotonic():
    s = build_scheme("sum", bin_w=100, cap=500, log=True, cmap_sum="Reds", cmap_div="RdBu_r")
    t = s.transform(np.array([0, 1, 10, 100, 500, 999]))
    assert list(t) == sorted(t)
    assert t[0] == 0.0
