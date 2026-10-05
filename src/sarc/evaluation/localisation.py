from collections import Counter
from . import prediction_specs as original_specs
from . import prediction as base
from .prediction_specs import SEED, SELECTION_RULE, specifications, build_estimator

METHODS = ("DawnRank", "PersonaDrive", "PRODIGY")
MASTERS = (104729, 130363, 155921, 196613, 228017)
COHORTS = ("COAD_CCLE", "LUAD_CCLE")
TRANSITIONS = ("11_0_to_11_5", "11_5_to_12_0")
STRICT = (
    "added_mutant_contact_edges",
    "removed_mutant_contact_edges",
    "mutant_candidates_touched_fraction",
    "mean_mutant_relative_edit_burden",
    "max_mutant_relative_edit_burden",
    "mean_mutant_deg_onehop_edit_burden",
)
ASSISTED = (
    "top10_touched_fraction",
    "next10_touched_fraction",
    "top10_mean_relative_edit_burden",
    "next10_mean_relative_edit_burden",
    "top10_minus_next10_relative_edit_burden",
    "boundary_margin_over_top20_range",
)
TRACKS = {
    "localisation_only": list(STRICT),
    "strict_plus_localisation": [*original_specs.FEATURES, *STRICT],
    "input_plus_localisation_and_old_output": [
        *original_specs.FEATURES,
        *STRICT,
        "old_native_rank_margin",
        *ASSISTED,
    ],
}
IDENTITIES = (
    "row_id",
    "cohort",
    "network_id",
    "method",
    "sample_id",
    "transition",
    "control_seed",
    "master_seeds",
)

def vector(row, track):
    names = (
        TRACKS["input_plus_localisation_and_old_output"]
        if (track == "input_plus_margin")
        else TRACKS[track]
    )
    values = []
    for (name) in (names):
        value = row["features"].get(name)
        if (name == "mean_consumed_cost_change" and row["method"] != "PRODIGY"):
            base.require(
                value is None
                and row["features"].get("cost_structurally_not_applicable") is True,
                "Binary structural cost differs",
            )
            value = 0
        if (value is None):
            return None
        values.append(original_specs.number(value))
    return values + [
        float(row["method"] == "DawnRank"),
        float(row["method"] == "PersonaDrive"),
    ]

def join_rows(original, features, production = True):
    base.validate_rows(original, production)
    base.require(len(original) == len(features), "Exact feature population required")
    by_id = {row["row_id"]: row for (row) in (features)}
    base.require(
        len(by_id) == len(features)
        and set(by_id) == {row["row_id"] for (row) in (original)},
        "Missing, extra or duplicate feature identity",
    )
    joined, counts = ([], Counter())
    for (row) in (original):
        extra = by_id[row["row_id"]]
        base.require(
            all((extra[key] == row[key] for (key) in (IDENTITIES))),
            "Feature/target identity or complete seed axis differs",
        )
        combined = dict(row["features"])
        base.require(
            not set(combined).intersection((*STRICT, *ASSISTED)),
            "New features overwrite original fields",
        )
        for (block_name, names) in ((
            ("strict_input", STRICT),
            ("old_output_assisted", ASSISTED),
        )):
            block = extra[block_name]
            base.require(
                set(block["features"]) == set(names),
                "Exact proposed feature keys required",
            )
            base.require(
                block["status"] in ("AVAILABLE", "UNAVAILABLE"),
                "Invalid feature availability",
            )
            if (block["status"] == "AVAILABLE"):
                base.require(
                    all((block["features"][name] is not None for (name) in (names))),
                    "Available feature block contains missingness",
                )
                for (value) in (block["features"].values()):
                    original_specs.number(value)
            else:
                base.require(
                    block["reason"]
                    and all((block["features"][name] is None for (name) in (names))),
                    "Unavailable feature block must preserve explicit missingness",
                )
            base.require(
                block_name != "strict_input" or block["status"] == "AVAILABLE",
                "Accepted strict features must be complete",
            )
            counts[block_name + ":" + block["status"]] += 1
            combined.update(block["features"])
        joined.append({**row, "features": combined})
    if (production):
        base.require(
            counts
            == Counter(
                {
                    "strict_input:AVAILABLE": 13824,
                    "old_output_assisted:AVAILABLE": 13440,
                    "old_output_assisted:UNAVAILABLE": 384,
                }
            ),
            "Accepted feature availability differs",
        )
    return (
        joined,
        {
            "rows": len(joined),
            "feature_blocks": dict(counts),
            "track_available": {
                track: sum((vector(row, track) is not None for (row) in (joined)))
                for (track) in (TRACKS)
            },
            "original_nonfeature_fields_unchanged": True,
            "original_feature_values_unchanged": True,
        },
    )
