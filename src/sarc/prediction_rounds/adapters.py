from types import SimpleNamespace
from ..evaluation import localisation
from ..evaluation.prediction import require
from . import specs_r1, models_r3, models_r4, specs_r5, models_r6, models_r7, fit_r7

ASSISTED_TRACKS = (
    "old_output_margin",
    "input_plus_margin",
    "input_plus_localisation_and_old_output",
)
EXPECTED_TASKS = {
    "R1": (1152, 240),
    "R2": (576, 120),
    "R3": (576, 108),
    "R4": (432, 108),
    "R5": (192, 48),
    "R6": (336, 84),
    "R7": (192, 32),
}

def vector(row, names):
    require(row["method"] in specs_r1.METHODS, "Unknown method")
    values = []
    for (name) in (names):
        value = row["features"].get(name)
        if (name == "mean_consumed_cost_change" and row["method"] != "PRODIGY"):
            require(
                value is None
                and row["features"].get("cost_structurally_not_applicable") is True,
                "Binary structural cost differs",
            )
            value = 0.0
        if (value is None):
            return None
        values.append(specs_r1.number(value))
    return values + [
        float(row["method"] == "DawnRank"),
        float(row["method"] == "PersonaDrive"),
    ]

def adapter(round_id):
    require(
        round_id in EXPECTED_TASKS,
        "Only the seven finite recorded rounds are supported",
    )
    inherited = {
        name: getattr(specs_r1, name)
        for (name) in ((
            "SEED",
            "FEATURES",
            "METHODS",
            "MASTERS",
            "COHORTS",
            "TRANSITIONS",
            "EXPECTED_RUNTIME",
            "EXPECTED_PYTHON",
            "SELECTION_RULE",
        ))
    }
    model = {
        "R1": specs_r1,
        "R2": specs_r1,
        "R3": models_r3,
        "R4": models_r4,
        "R5": specs_r5,
        "R6": models_r6,
        "R7": fit_r7,
    }[round_id]
    tracks = {
        "R1": specs_r1.TRACKS,
        "R2": localisation.TRACKS,
        "R3": {**specs_r1.TRACKS, **localisation.TRACKS},
        "R4": {**specs_r1.TRACKS, **localisation.TRACKS},
        "R5": specs_r5.TRACKS,
        "R6": models_r6.TRACKS,
        "R7": models_r7.TRACKS,
    }[round_id]
    result = SimpleNamespace(
        **inherited,
        round_id = round_id,
        TRACKS = {name: list(values) for (name, values) in (tracks.items())},
        specifications = model.specifications,
    )
    result.vector = lambda row, track: vector(row, result.TRACKS[track])
    result.build_estimator = getattr(model, "build_estimator", None)
    return result

def validate_partition_order(rows, order):
    require(
        set(order) == set(specs_r1.COHORTS),
        "Both archived development-cohort priority lists required",
    )
    for (cohort) in (specs_r1.COHORTS):
        samples = {row["sample_id"] for (row) in (rows) if (row["cohort"] == cohort)}
        require(
            len(order[cohort]) == len(set(order[cohort])) == 36
            and set(order[cohort]) == samples,
            "Archived sample priorities differ from complete biological population",
        )
