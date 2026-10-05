from __future__ import annotations

EPSILONS = [0.0, 1e-12, 1e-08, 1e-06, 0.0001, 0.001]
TAXONOMY = [
    "IDENTICAL",
    "TIE_EXPANSION",
    "TIE_CONTRACTION",
    "OVERLAPPING_LEADER_SET",
    "UNIQUE_LEADER_REPLACEMENT",
    "DISJOINT_REPLACEMENT",
]
MARGIN_LABELS = [
    "0",
    "(0, 1e-12]",
    "(1e-12, 1e-8]",
    "(1e-8, 1e-6]",
    "(1e-6, 1e-4]",
    "(1e-4, 1e-3]",
    "(1e-3, 1e-2]",
    "(1e-2, 5e-2]",
    "(5e-2, 1]",
]
