from __future__ import annotations
from typing import Any, Mapping, Sequence
from .sensitivity_types import CrossLayerSensitivityBuildError
from .ranking import (
    finite_float,
    jaccard_overlap,
    kendall_tau_b,
    linear_quantile,
    normalised_rank_displacements,
    parse_boolean,
    spearman_correlation,
    tie_aware_top_ranked_set,
    tie_aware_top_set,
)

def joined_identifiers(values: set[str] | None) -> str:
    return "" if (values is None) else ";".join(sorted(values))

def split_identifiers(value: Any) -> set[str]:
    if (value in {None, ""}):
        return set()
    return {part for (part) in (str(value).split(";")) if (part)}

def proportion(values: Sequence[bool]) -> tuple[int, int, float | None]:
    numerator = sum((bool(value) for (value) in (values)))
    denominator = len(values)
    return (
        numerator,
        denominator,
        None if (denominator == 0) else numerator / denominator,
    )

def calculate_fixed_panel_metrics(
    rows: Sequence[Mapping[str, Any]], configuration: Mapping[str, Any], layer: str
) -> dict[str, Any]:
    if (not rows):
        raise CrossLayerSensitivityBuildError(
            "Cannot calculate an empty Cross-layer sensitivity panel"
        )
    targets = [str(row["canonical_target_id"]) for (row) in (rows)]
    if (len(targets) != len(set(targets))):
        raise CrossLayerSensitivityBuildError(
            "Cross-layer sensitivity fixed panel contains duplicate canonical targets"
        )
    old_scores = [finite_float(row["old_score"]) for (row) in (rows)]
    new_scores = [finite_float(row["new_score"]) for (row) in (rows)]
    lower = finite_float(configuration["scores"]["minimum"])
    upper = finite_float(configuration["scores"]["maximum"])
    if (any((not lower <= score <= upper for (score) in (old_scores + new_scores)))):
        raise CrossLayerSensitivityBuildError(
            "Cross-layer sensitivity fixed panel contains an out-of-bounds score"
        )
    if (layer == "STABLE_ENTITY_FIXED" and any(
        (score <= lower for (score) in (old_scores + new_scores))
    )):
        raise CrossLayerSensitivityBuildError(
            "Stable-entity overall support contains a non-positive score"
        )
    minimum_targets = int(configuration["ranking"]["minimum_targets"])
    minimum_distinct = int(
        configuration["ranking"]["minimum_distinct_scores_per_release"]
    )
    if (len(targets) < minimum_targets):
        raise CrossLayerSensitivityBuildError(
            "Cross-layer sensitivity fixed panel is below the frozen target minimum"
        )
    if (
        len(set(old_scores)) < minimum_distinct
        or len(set(new_scores)) < minimum_distinct
    ):
        raise CrossLayerSensitivityBuildError(
            "Cross-layer sensitivity fixed panel is score-degenerate"
        )
    old_ranks, new_ranks, displacements = normalised_rank_displacements(
        old_scores, new_scores
    )
    old_map = dict(zip(targets, old_scores))
    new_map = dict(zip(targets, new_scores))
    minimum_top_score = finite_float(
        configuration["ranking"]["minimum_positive_top_set_score"]
    )
    result: dict[str, Any] = {
        "analysis_layer": layer,
        "release_pair": rows[0]["release_pair"],
        "old_release": rows[0]["old_release"],
        "new_release": rows[0]["new_release"],
        "canonical_disease_id": rows[0]["canonical_disease_id"],
        "canonical_label": rows[0].get("canonical_label", ""),
        "canonical_therapeutic_area_ids": rows[0].get(
            "canonical_therapeutic_area_ids", ""
        ),
        "fixed_target_count": len(targets),
        "old_distinct_score_count": len(set(old_scores)),
        "new_distinct_score_count": len(set(new_scores)),
        "old_positive_score_count": sum(
            (score > minimum_top_score for (score) in (old_scores))
        ),
        "new_positive_score_count": sum(
            (score > minimum_top_score for (score) in (new_scores))
        ),
        "kendall_tau_b": kendall_tau_b(old_scores, new_scores),
        "spearman_correlation": spearman_correlation(old_scores, new_scores),
        "median_normalised_rank_displacement": linear_quantile(displacements, 0.5),
        "p90_normalised_rank_displacement": linear_quantile(displacements, 0.9),
        "maximum_normalised_rank_displacement": max(displacements),
        "stable_entity_member_fraction": sum(
            (
                parse_boolean(
                    row.get(
                        "stable_entity_eligible", row.get("stable_entity_member", False)
                    )
                )
                for (row) in (rows)
            )
        )
        / len(rows),
    }
    for (cutoff) in (configuration["ranking"]["top_set_cutoffs"]):
        old_set = tie_aware_top_set(old_map, int(cutoff), minimum_top_score)
        new_set = tie_aware_top_set(new_map, int(cutoff), minimum_top_score)
        result[f"top_{cutoff}_jaccard"] = jaccard_overlap(old_set, new_set)
        result[f"old_top_{cutoff}_effective_size"] = (
            None if (old_set is None) else len(old_set)
        )
        result[f"new_top_{cutoff}_effective_size"] = (
            None if (new_set is None) else len(new_set)
        )
        result[f"top_{cutoff}_changed"] = (
            None if (old_set is None or new_set is None) else old_set != new_set
        )
    old_top = tie_aware_top_ranked_set(old_map, minimum_top_score)
    new_top = tie_aware_top_ranked_set(new_map, minimum_top_score)
    result.update(
        {
            "old_top_ranked_target_ids": joined_identifiers(old_top),
            "new_top_ranked_target_ids": joined_identifiers(new_top),
            "old_top_ranked_set_size": None if (old_top is None) else len(old_top),
            "new_top_ranked_set_size": None if (new_top is None) else len(new_top),
            "top_ranked_jaccard": jaccard_overlap(old_top, new_top),
            "top_ranked_set_changed": None
            if (old_top is None or new_top is None)
            else old_top != new_top,
            "top_ranked_sets_disjoint": None
            if (old_top is None or new_top is None)
            else not bool(old_top & new_top),
            "support_denominator": layer,
        }
    )
    if (result["kendall_tau_b"] is None or result["spearman_correlation"] is None):
        raise CrossLayerSensitivityBuildError(
            "An eligible Cross-layer sensitivity panel has a non-estimable rank correlation"
        )
    return result

def classify_persistence(
    primary_changed: bool | None,
    stable_eligible: bool,
    stable_changed: bool | None,
    domain_eligible: bool,
    domain_changed: bool | None,
) -> str:
    primary = bool(primary_changed)
    stable = stable_eligible and bool(stable_changed)
    domain = domain_eligible and bool(domain_changed)
    if (primary and stable and domain):
        return "MIXED"
    if (primary and stable):
        return "PERSISTS_STABLE_ENTITY"
    if (primary and domain):
        return "PERSISTS_EVIDENCE_DOMAIN"
    if (primary and (not domain_eligible)):
        return "INSUFFICIENT_DOMAIN_SUPPORT"
    if (primary and stable_eligible and (not stable) and (not domain)):
        return "DISAPPEARS_AFTER_MAPPING_CONTROL"
    if (primary):
        return "OVERALL_ONLY"
    if (domain and (not stable)):
        return "DOMAIN_ONLY"
    return "MIXED"

def release_note_alignment_class(
    category: str, evidence_status: str, configuration: Mapping[str, Any]
) -> str:
    rules = configuration["release_note_alignment"]
    if (evidence_status == str(rules["no_specific_change_status"])):
        return "NOT TESTED"
    if (category in set(rules["directly_observable_categories"])):
        return "DIRECTLY OBSERVABLE IN RELEASE OUTPUT"
    if (category in set(rules["possibly_relevant_categories"])):
        return "POSSIBLY RELEVANT"
    if (category in set(rules["infrastructure_categories"])):
        return "INFRASTRUCTURE ONLY"
    return "NOT TESTED"
