from fractions import Fraction

from ..precision.io import require

def statistics(values):
    require(bool(values), "Empty donor vector")
    denominator = len(values)
    positive = sum((value for (value) in (values) if (value > 0)), Fraction(0))
    negative = sum((-value for (value) in (values) if (value < 0)), Fraction(0))
    total = sum(values, Fraction(0))
    absolute = sum((abs(value) for (value) in (values)), Fraction(0))
    ordered = sorted(values)
    middle = denominator // 2
    median = (
        ordered[middle] if (denominator % 2) else (ordered[middle - 1] + ordered[middle]) / 2
    )
    require(
        total == positive - negative and absolute == positive + negative,
        "Signed mass identity failed",
    )
    return {
        "positive_sum": positive,
        "negative_sum": negative,
        "absolute_sum": absolute,
        "signed_sum": total,
        "positive_mass": positive / denominator,
        "negative_mass": negative / denominator,
        "absolute_mass": absolute / denominator,
        "mean": total / denominator,
        "minimum": ordered[0],
        "median": median,
        "maximum": ordered[-1],
        "n_negative": sum(value < 0 for (value) in (values)),
        "n_zero": sum(value == 0 for (value) in (values)),
        "n_positive": sum(value > 0 for (value) in (values)),
    }

def summarise(donor_rows):
    keys = {(row["h"], row["category"]) for (row) in (donor_rows)}
    summaries, vectors = [], []
    for (h, category) in (sorted(keys)):
        selected = [
            row for (row) in (donor_rows) if (row["h"] == h and row["category"] == category)
        ]
        require(
            len(selected) == 83 and sum(row["models"] for (row) in (selected)) == 85,
            "Incomplete donor population",
        )
        for (term) in (("I_total", "I_score_reference", "I_decision_reference")):
            values = [row[term] for (row) in (selected)]
            summaries.append(
                {
                    "h": h,
                    "category": category,
                    "term": term,
                    "n_donors": 83,
                    "n_models": 85,
                    **statistics(values),
                }
            )
            for (row, value) in (zip(selected, values)):
                vectors.append(
                    {
                        "donor_id": row["donor_id"],
                        "models": row["models"],
                        "h": h,
                        "category": category,
                        "term": term,
                        "value": value,
                        "positive": max(value, Fraction(0)),
                        "negative": max(-value, Fraction(0)),
                        "absolute": abs(value),
                    }
                )
    return summaries, vectors
