import json
import math

def exact(x):
    return json.dumps(x, sort_keys = True, separators = (",", ":"), allow_nan = False)

def finite(x):
    return type(x) in (int, float) and math.isfinite(x)

def project(entry):
    r = entry["row"]
    family = entry["family"]
    summary = r.get("summary", {})
    if (family == "B"):
        assert (
            r["event_capture_direction"] == "MAXIMISE"
            and r["brier_direction"] == "MINIMISE"
        )
        contract = {
            "family": family,
            "comparison": r["comparison_key"],
            "information_regime": r["information_regime"],
            "counts": {
                k: r[k]
                for (k) in (
                    (
                        "n_biological_full",
                        "n_biological_with_available_actions",
                        "n_logical_actions_full",
                        "n_logical_actions_available",
                        "n_allocated_actions",
                        "n_events",
                    )
                )
            },
        }
        axes = [
            ("event_capture", "max", r["event_capture"]),
            ("brier", "min", r["brier"]),
        ]
        policy = {
            "policy_id": r["family"],
            "information_regime": r["information_regime"],
            "deployable": None,
            "original_row_binding": r["source_binding"],
        }
        applicable = True
    elif (family == "FIXED_POLICY_SELECTION"):
        assert (
            r["selection_loss"]["direction"] == "minimise"
            and r["comparison"]["utility_direction"] == "minimise"
        )
        contract = {
            "family": family,
            "comparison": r["comparison"],
            "uncertainty_convention": {
                k: v
                for ((k, v)) in (r.get("uncertainty", {}).items())
                if (k in ("type", "scope", "method", "resampling_unit", "exposure"))
            },
        }
        axes = [("fixed_policy_selection_loss", "min", r["selection_loss"]["value"])]
        policy = r["policy"]
        applicable = finite(r["selection_loss"]["value"])
    else:
        assert (family in ("A", "E", "C_CANDIDATE", "C_SELECTION", "C_REVIEW"))
        contract = {
            "family": family,
            "comparison": r["comparison"],
            "counts": {
                k: summary.get(k)
                for (k) in (
                    (
                        "n_full",
                        "n_actionable",
                        "n_prediction_available",
                        "n_assay_eligible",
                        "utility_denominator",
                    )
                )
            },
            "point_rule": summary.get("utility_point_rule"),
            "exposure_binding": summary.get("exposure_binding"),
        }
        axes = [
            (r["comparison"]["endpoint_construct"], "max", summary["utility_value"])
        ]
        policy = r["policy"]
        applicable = policy["applicability_status"] in (
            "applicable",
            "AVAILABLE_POST_COMPUTATION",
            "AVAILABLE",
        )
    identity = entry["catalog_id"]
    return {
        "catalog_id": identity,
        "family": family,
        "original_row_id": entry["original_row_id"],
        "source_binding": entry["source_binding"],
        "source_pointer": entry["source_pointer"],
        "panel_id": exact(contract),
        "contract": contract,
        "policy": policy,
        "axes": [{"name": n, "direction": d, "value": v} for ((n, d, v)) in (axes)],
        "applicable": applicable,
        "cost_status": "END_TO_END_UNAVAILABLE",
        "end_to_end_cost": None,
        "utility_cost_classification": "UNCLASSIFIED",
        "original_summary_retained_by_reference": True,
        "original_alias_metadata": policy.get("alias_of", policy.get("alias_metadata")),
    }

def compare(members):
    valid = [
        r
        for (r) in (members)
        if (r["applicable"] and all((finite(a["value"]) for (a) in (r["axes"]))))
    ]
    for (r) in (members):
        r["descriptive_point"] = {
            "status": "UNCLASSIFIED",
            "dominated_by": [],
            "equal_to": [],
            "comparable_available_rows": len(valid),
            "comparative_evidence": len(valid) > 1,
        }
        if (r not in valid):
            continue
        current = [
            a["value"] * (1 if (a["direction"] == "max") else -1) for (a) in (r["axes"])
        ]
        for (other) in (valid):
            if (other["catalog_id"] == r["catalog_id"]):
                continue
            assert ([(a["name"], a["direction"]) for (a) in (r["axes"])] == [
                (a["name"], a["direction"]) for (a) in (other["axes"])
            ])
            value = [
                a["value"] * (1 if (a["direction"] == "max") else -1)
                for (a) in (other["axes"])
            ]
            if (value == current):
                r["descriptive_point"]["equal_to"].append(other["catalog_id"])
            elif (all((a >= b for ((a, b)) in (zip(value, current)))) and any(
                (a > b for ((a, b)) in (zip(value, current)))
            )):
                r["descriptive_point"]["dominated_by"].append(other["catalog_id"])
        r["descriptive_point"]["status"] = (
            "DOMINATED"
            if (r["descriptive_point"]["dominated_by"])
            else "NONDOMINATED_DESCRIPTIVE_POINT"
        )
    return members
