from dataclasses import dataclass
from math import isfinite
from .comparators import compare_panels
from .contracts import ContractError, PanelSnapshot, require_text, unique_identifiers

@dataclass(frozen = True)
class PersonaDrivePanel(PanelSnapshot):
    def __post_init__(self):
        if (self.resource != "PersonaDrive"):
            raise ContractError("Unexpected replay resource")
        require_text(self.release, "release")
        require_text(self.panel_id, "sample")
        if (not isinstance(self.members, tuple) or self.identity_decidable is not True):
            raise ContractError("Replay requires explicit resolved input identities")
        unique_identifiers((gene for ((gene, _)) in (self.members)), "gene")
        for (_, value) in (self.members):
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or (not isfinite(value))
                or (value < 0)
            ):
                raise ContractError(
                    "Influence scores must be finite nonnegative numbers"
                )

def ordered_genes(scores):
    return sorted(scores, key = lambda gene: (-scores[gene], gene))

def fraction(values):
    observed = [value for (value) in (values) if (value is not None)]
    if (any((not isinstance(value, bool) for (value) in (observed)))):
        raise ContractError("Event summaries require boolean observations")
    if (not observed):
        return {
            "value": None,
            "numerator": None,
            "denominator": None,
            "support_denominator": 0,
            "reason": "NO_EVALUABLE_SAMPLE",
            "unevaluable_samples": len(values),
        }
    count = sum(observed)
    return {
        "value": count / len(observed),
        "numerator": count,
        "denominator": len(observed),
        "support_denominator": len(observed),
        "reason": None,
        "unevaluable_samples": len(values) - len(observed),
    }

def original_output_comparison(first, second, top_ks):
    (left, right) = (ordered_genes(first), ordered_genes(second))
    persistent = first.keys() & second.keys()
    first_persistent = [gene for (gene) in (left) if (gene in persistent)]
    second_persistent = [gene for (gene) in (right) if (gene in persistent)]
    top_lists = {}
    for (top_k) in (top_ks):
        enough = min(len(left), len(right)) >= top_k
        selected_first = left[:top_k] if (len(left) >= top_k) else None
        selected_second = right[:top_k] if (len(right) >= top_k) else None
        top_lists[str(top_k)] = {
            "baseline_order": selected_first,
            "followup_order": selected_second,
            "membership_changed": set(selected_first) != set(selected_second)
            if (enough)
            else None,
            "order_changed": selected_first != selected_second if (enough) else None,
            "reason": None if (enough) else "INSUFFICIENT_CANDIDATES",
        }
    return {
        "baseline_order": left,
        "followup_order": right,
        "baseline_output_state": "EMITTED"
        if (left)
        else "NOT_EMITTED_EMPTY_CANDIDATE_GRAPH",
        "followup_output_state": "EMITTED"
        if (right)
        else "NOT_EMITTED_EMPTY_CANDIDATE_GRAPH",
        "complete_ordered_output_changed": left != right,
        "candidate_roster_changed": set(left) != set(right),
        "strict_leader_changed": left[0] != right[0] if (left and right) else None,
        "persistent_gene_order_changed": first_persistent != second_persistent
        if (persistent)
        else None,
        "persistent_order_reason": None if (persistent) else "NO_COMMON_CANDIDATES",
        "baseline_persistent_order": first_persistent,
        "followup_persistent_order": second_persistent,
        "top_k": top_lists,
    }

def compare_cohort(baseline, followup, cohort, settings):
    expected = unique_identifiers(cohort, "cohort sample")
    if (set(baseline) != expected or set(followup) != expected):
        raise ContractError("Every original sample needs an explicit completed output")
    rules = settings["common_comparator"]
    minimums = sorted(
        {rules["primary_minimum_roster"], *rules["secondary_minimum_rosters"]}
    )
    top_ks = tuple(rules["top_k"])
    if (not top_ks or any(
        (isinstance(k, bool) or not isinstance(k, int) or k < 1 for (k) in (top_ks))
    )):
        raise ContractError("Top-k values must be positive integers")
    records = []
    for (sample) in (cohort):
        if (not isinstance(baseline[sample], dict) or not isinstance(
            followup[sample], dict
        )):
            raise ContractError("Sample output requires a gene-to-score object")
        first = PersonaDrivePanel(
            "PersonaDrive", "STRING_11.5", sample, tuple(baseline[sample].items())
        )
        second = PersonaDrivePanel(
            "PersonaDrive", "STRING_12.0", sample, tuple(followup[sample].items())
        )
        original = original_output_comparison(
            baseline[sample], followup[sample], top_ks
        )
        reports = {
            str(minimum): compare_panels(
                first,
                second,
                minimum_roster = minimum,
                top_ks = top_ks,
                rbo_p = tuple(rules["RBO_persistence"]),
                rbo_variants = tuple(rules["RBO_variants"]),
            )
            for (minimum) in (minimums)
        }
        union = sorted(baseline[sample].keys() | followup[sample].keys())
        transitions = []
        for (gene) in (union):
            left = baseline[sample].get(gene)
            right = followup[sample].get(gene)
            transitions.append(
                {
                    "gene_id": gene,
                    "baseline_score": left,
                    "followup_score": right,
                    "baseline_state": "ABSENT_CANDIDATE"
                    if (left is None)
                    else "PRESENT_ZERO"
                    if (left == 0)
                    else "PRESENT_POSITIVE",
                    "followup_state": "ABSENT_CANDIDATE"
                    if (right is None)
                    else "PRESENT_ZERO"
                    if (right == 0)
                    else "PRESENT_POSITIVE",
                    "signed_score_change": None
                    if (left is None or right is None)
                    else right - left,
                }
            )
        records.append(
            {
                "sample_id": sample,
                "original_task": original,
                "gene_states": transitions,
                "comparisons": reports,
            }
        )
    summary = {
        name: fraction([row["original_task"][name] for (row) in (records)])
        for (name) in (
            (
                "complete_ordered_output_changed",
                "candidate_roster_changed",
                "strict_leader_changed",
                "persistent_gene_order_changed",
            )
        )
    }
    for (top_k) in (top_ks):
        for (event) in (("membership_changed", "order_changed")):
            summary[f"strict_top_{top_k}_{event}"] = fraction(
                [
                    row["original_task"]["top_k"][str(top_k)][event]
                    for (row) in (records)
                ]
            )
    shared_summary = {}
    for (minimum) in (minimums):
        section = {}
        for (support) in (("native", "fixed")):
            for (endpoint) in (("leader", *[f"top_{k}" for (k) in (top_ks)])):
                values = []
                for (row) in (records):
                    report = row["comparisons"][str(minimum)][support]
                    value = (
                        report["leader"]
                        if (endpoint == "leader")
                        else report["top_k"][endpoint.removeprefix("top_")]
                    )
                    values.append(value["value"] if (value["estimable"]) else None)
                section[f"{support}_{endpoint}_changed"] = fraction(values)
        shared_summary[str(minimum)] = section
    return {
        "resource": "PersonaDrive",
        "case": settings["case"],
        "samples": len(cohort),
        "original_task_summary": summary,
        "common_comparator_summary": shared_summary,
        "score_semantics": settings["score_semantics"],
        "fixed_semantics": settings["fixed_semantics"],
        "uncertainty": settings["uncertainty"],
        "records": records,
    }
