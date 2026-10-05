from fractions import Fraction

TERMS = ("I_total", "I_score_reference", "I_decision_reference")
CATEGORIES = {"ALL", "BOTH_CORE", "NOT_BOTH_CORE"}

def checked_rows(rows):
    keys = [(r["donor_id"], r["h"], r["category"]) for (r) in (rows)]
    if ((len(rows) != 498) or (len(set(keys)) != 498)):
        raise ValueError("Incomplete or duplicate donor-state-category support")
    roster = {r["donor_id"]: int(r["models"]) for (r) in (rows)}
    if ((len(roster) != 83) or (sum(roster.values()) != 85)):
        raise ValueError("Invalid 85-model 83-donor population")
    if (any(count < 1 for (count) in (roster.values()))):
        raise ValueError("Every donor must retain at least one model")
    for (h) in (("0", "1")):
        for (category) in (sorted(CATEGORIES)):
            selected = [r for (r) in (rows) if ((r["h"] == h) and (r["category"] == category))]
            if ({r["donor_id"] for (r) in (selected)} != set(roster)):
                raise ValueError("Background/category roster mismatch")
            for (row) in (selected):
                if (int(row["models"]) != roster[row["donor_id"]]):
                    raise ValueError("Within-donor model count mismatch")
                values = [Fraction(row[k]) for (k) in (TERMS)]
                if (values[0] != sum(values[1:])):
                    raise ValueError("Saved signed decomposition mismatch")
    return roster

def statistics(values):
    if (not values):
        raise ValueError("Empty donor vector")
    denominator = len(values)
    positive = sum((v for (v) in (values) if (v > 0)), Fraction())
    negative = sum((-v for (v) in (values) if (v < 0)), Fraction())
    total = sum(values, Fraction())
    absolute = sum((abs(v) for (v) in (values)), Fraction())
    ordered = sorted(values)
    middle = denominator // 2
    median = (
        ordered[middle] if (denominator % 2) else (ordered[middle - 1] + ordered[middle]) / 2
    )
    if ((total != positive - negative) or (absolute != positive + negative)):
        raise ValueError("Mass identity failed")
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
        "n_negative": sum(v < 0 for (v) in (values)),
        "n_zero": sum(v == 0 for (v) in (values)),
        "n_positive": sum(v > 0 for (v) in (values)),
    }

def summarise(rows):
    checked_rows(rows)
    selected = [row for (row) in (rows) if (row["category"] == "BOTH_CORE")]
    vectors, summaries = [], []
    for (background) in (("0", "1")):
        donors = sorted(
            [row for (row) in (selected) if (row["h"] == background)],
            key = lambda row: row["donor_id"],
        )
        for (term) in (TERMS):
            values = [Fraction(row[term]) for (row) in (donors)]
            summaries.append(
                {
                    "h": background,
                    "category": "BOTH_CORE",
                    "term": term,
                    "n_donors": 83,
                    "n_models": 85,
                    **statistics(values),
                }
            )
            for (row, value) in (zip(donors, values)):
                vectors.append(
                    {
                        "donor_id": row["donor_id"],
                        "models": row["models"],
                        "h": background,
                        "category": "BOTH_CORE",
                        "term": term,
                        "value_exact": value,
                        "positive_exact": max(value, Fraction()),
                        "negative_exact": max(-value, Fraction()),
                        "absolute_exact": abs(value),
                    }
                )
    return {
        "population": {"models": 85, "donors": 83, "backgrounds": 2},
        "summary": summaries,
    }, vectors
