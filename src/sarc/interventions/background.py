import math
from collections import Counter
from fractions import Fraction

COHORTS = {"COAD_CCLE": 36, "LUAD_CCLE": 36, "COAD_TCGA": 396}
TRANSITIONS = ("11_0_to_11_5", "11_5_to_12_0")
CLASSES = ("PB", "PW", "DB", "DW", "PB_PW")
BUDGETS = (1, 5, 10, 20)
DIRECTIONS = ("insert", "retract")

def require(condition, message):
    if (not condition):
        raise ValueError(message)

def rational(value):
    if (isinstance(value, Fraction)):
        return value
    if (type(value) is int):
        return Fraction(value)
    require(
        isinstance(value, dict)
        and set(value)
        in ({"numerator", "denominator"}, {"numerator", "denominator", "value"}),
        "Exact rational required; unsupported wrapper/float/bool forbidden",
    )
    require(
        type(value["numerator"]) is int
        and type(value["denominator"]) is int
        and (value["denominator"] > 0),
        "Invalid exact rational",
    )
    result = Fraction(value["numerator"], value["denominator"])
    if ("value" in value):
        require(
            result.numerator == value["numerator"]
            and result.denominator == value["denominator"],
            "Noncanonical accepted-primary exact rational",
        )
        require(
            type(value["value"]) in (int, float)
            and math.isfinite(value["value"])
            and (value["value"] == float(result)),
            "Inconsistent accepted-primary rational cache",
        )
    return result

def nullable(value):
    return None if (value is None) else rational(value)

def valid_genes(value):
    return (
        isinstance(value, list)
        and all(isinstance(gene, str) and gene for (gene) in (value))
        and value == sorted(set(value))
    )

def endpoint(record, k):
    comparison = record["comparison"]
    require(
        comparison["method"] == "DawnRank" and comparison["master_seeds"] == [None],
        "Dawn single-state policy required",
    )
    ep = comparison["endpoints"]["common::NCG6_primary_all::" + str(k)]
    require(
        valid_genes(ep["support_genes"]) and valid_genes(ep["positive_genes"]),
        "Missing explicit endpoint gene sets",
    )
    sides = {}
    for (name) in (("baseline", "intervention")):
        item = ep[name]
        require(
            item["status"] in ("AVAILABLE", "UNAVAILABLE"),
            "Unknown endpoint availability",
        )
        hits = nullable(item["hits"])
        require(
            (hits is not None) == (item["status"] == "AVAILABLE"),
            "Availability/hit mismatch",
        )
        require(
            len(item["seed_hits"]) == 1 and nullable(item["seed_hits"][0]) == hits,
            "Dawn seed reduction differs",
        )
        require(hits is None or 0 <= hits <= k, "Hit level outside offered quota")
        sides[name] = {
            "hits": hits,
            "status": item["status"],
            "native_output_count": nullable(item["native_output_count"]),
            "eligible_output_count": nullable(item["eligible_output_count"]),
            "eligible_count": item["eligible_count"],
            "positive_count": item["positive_count"],
        }
    delta = nullable(ep["effect"]["signed_mean_delta"])
    available = all((sides[name]["hits"] is not None for (name) in (sides)))
    require(
        (delta is not None) == available,
        "Saved effect availability differs from its corners",
    )
    if (available):
        require(
            delta == sides["intervention"]["hits"] - sides["baseline"]["hits"],
            "Saved effect/corner identity mismatch",
        )
    return {
        "gene_sets": (ep["support_genes"], ep["positive_genes"]),
        "delta": delta,
        "sides": sides,
    }

def classify(a, b):
    old = None if (a is None) else "LOSS" if (a <= -1) else "PASS"
    new = None if (b is None) else "LOSS" if (b <= -1) else "PASS"
    complete = a is not None and b is not None
    return {
        "old_background": old,
        "new_background": new,
        "cell": old + "__" + new if (complete) else None,
        "false_pass": old == "PASS" and new == "LOSS" if (complete) else None,
        "reverse_compensation": old == "LOSS" and new == "PASS" if (complete) else None,
        "strong_false_pass": a >= 0 and b <= -1 if (complete) else None,
        "strong_reverse_compensation": a <= -1 and b >= 0 if (complete) else None,
    }

def paired_case(insert, retract, k):
    fields = ("cohort", "method", "sample_id", "kind", "class_id", "control_seed")
    require(
        all((insert[field] == retract[field] for (field) in (fields))),
        "Case identity mismatch",
    )
    require(
        insert["method"] == "DawnRank"
        and insert["kind"] == "real"
        and (insert["control_seed"] is None),
        "Only real Dawn rows allowed",
    )
    (first, second) = (endpoint(insert, k), endpoint(retract, k))
    require(
        first["gene_sets"] == second["gene_sets"],
        "Cross-background support/positive gene sets differ",
    )
    a = first["delta"]
    b = None if (second["delta"] is None) else -second["delta"]
    return {
        "sample_id": insert["sample_id"],
        "a": a,
        "b": b,
        "theta": b - a if (a is not None and b is not None) else None,
        "status": "AVAILABLE" if (a is not None and b is not None) else "UNAVAILABLE",
        "support_genes": first["gene_sets"][0],
        "positive_genes": first["gene_sets"][1],
        "classification": classify(a, b),
        "corners": {
            "G00": first["sides"]["baseline"],
            "G10": first["sides"]["intervention"],
            "G01": second["sides"]["intervention"],
            "G11": second["sides"]["baseline"],
        },
    }

def distribution(values, total):
    available = [value for (value) in (values) if (value is not None)]
    counts = Counter(available)
    cumulative = 0
    rows = []
    for (value) in (sorted(counts)):
        cumulative += counts[value]
        rows.append(
            {
                "value": value,
                "count": counts[value],
                "cdf_among_available": Fraction(cumulative, len(available)),
                "observed_mass_over_full_population": Fraction(counts[value], total),
            }
        )
    return {
        "population_total": total,
        "available": len(available),
        "missing": total - len(available),
        "points": rows,
    }

def summarise(rows):
    total = len(rows)
    require(
        total > 0 and len({row["sample_id"] for (row) in (rows)}) == total,
        "Unique nonempty group required",
    )
    complete = [row for (row) in (rows) if (row["status"] == "AVAILABLE")]
    table = {
        a + "__" + b: 0 for (a) in (("PASS", "LOSS")) for (b) in (("PASS", "LOSS"))
    }
    for (row) in (complete):
        table[row["classification"]["cell"]] += 1
    old_pass = table["PASS__PASS"] + table["PASS__LOSS"]
    old_loss = table["LOSS__PASS"] + table["LOSS__LOSS"]
    events = {}
    for (key) in ((
        "false_pass",
        "reverse_compensation",
        "strong_false_pass",
        "strong_reverse_compensation",
    )):
        n = sum((row["classification"][key] is True for (row) in (complete)))
        events[key] = {
            "count": n,
            "rate_among_joint_available": Fraction(n, len(complete))
            if (complete)
            else None,
            "observed_count_over_full_population": Fraction(n, total),
            "full_population_rate_if_complete": Fraction(n, total)
            if (len(complete) == total)
            else None,
        }
    events["false_pass"]["conditional_on_old_pass_joint_available"] = (
        Fraction(table["PASS__LOSS"], old_pass) if (old_pass) else None
    )
    events["reverse_compensation"]["conditional_on_old_loss_joint_available"] = (
        Fraction(table["LOSS__PASS"], old_loss) if (old_loss) else None
    )
    return {
        "population_total": total,
        "joint_available": len(complete),
        "unavailable_sample_ids": [
            r["sample_id"] for (r) in (rows) if (r["status"] != "AVAILABLE")
        ],
        "table": table,
        "events": events,
        "distributions": {
            name: distribution([row[name] for (row) in (rows)], total)
            for (name) in (("a", "b", "theta"))
        },
        "inference": "DESCRIPTIVE_EXACT_ORIGINAL_SAMPLE_WEIGHTS_NO_BOOTSTRAP_OR_PVALUE",
        "cohorts_pooled": False,
    }
