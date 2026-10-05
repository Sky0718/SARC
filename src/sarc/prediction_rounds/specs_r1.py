from ..evaluation.prediction_specs import (
    SEED,
    SELECTION_RULE,
    FEATURES,
    TRACKS,
    specifications,
    number,
    vector,
    build_estimator,
)

METHODS = ("DawnRank", "PersonaDrive", "PRODIGY")
MASTERS = (104729, 130363, 155921, 196613, 228017)
COHORTS = ("COAD_CCLE", "LUAD_CCLE")
TRANSITIONS = ("11_0_to_11_5", "11_5_to_12_0")
EXPECTED_RUNTIME = {
    "numpy": "2.2.6",
    "scipy": "1.15.3",
    "sklearn": "1.7.2",
    "threadpoolctl": "3.6.0",
    "xgboost": "3.2.0",
    "catboost": "1.2.10",
}
EXPECTED_PYTHON = (3, 10, 6)
