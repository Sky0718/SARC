from fractions import Fraction

import numpy as np

from .endpoint import as_fraction, require
from .inputs import CATEGORIES, POPULATIONS

FIELDS = ("R", "H", "Z", "U", "absent", "explicit_null", "identity_conflict", "call_conflict")
FOCAL = ("BOTH_CORE", "ONCOGENE")

def rational(value):
    return as_fraction(value)

def validate_mass(row):
    values = [rational(row[field]) for (field) in (FIELDS)]
    if (any((value < 0) for (value) in (values))):
        raise ValueError("Negative mass")
    if (rational(row["R"]) != rational(row["H"]) + rational(row["Z"]) + rational(row["U"])):
        raise ValueError("Mass conservation")
    if (rational(row["U"]) != sum((rational(row[field]) for (field) in (FIELDS[4:])), Fraction(0))):
        raise ValueError("Unknown conservation")

def aggregate_donors(model_rows):
    grouped = {}
    identities = set()
    for (row) in (model_rows):
        identity = tuple(
            row[field] for (field) in (["context", "policy", "category", "state", "sample_id"])
        )
        if (identity in identities):
            raise ValueError("Duplicate model row")
        identities.add(identity)
        key = tuple(
            row[field] for (field) in (["context", "policy", "category", "state", "donor_id"])
        )
        grouped.setdefault(key, []).append(row)
    rows = []
    for (key, source) in (sorted(grouped.items())):
        row = dict(zip(["context", "policy", "category", "state", "donor_id"], key))
        row["model_count"] = len(source)
        row.update(
            {
                field: sum((rational(item[field]) for (item) in (source)), Fraction(0))
                / len(source)
                for (field) in (FIELDS)
            }
        )
        validate_mass(row)
        rows.append(row)
    return rows

def ratio_constant(r0, h0, r1, h1, value):
    for (i) in (range(len(r0))):
        if ((h1[i] * r0[i] - h0[i] * r1[i] - value * r1[i] * r0[i]) != 0):
            return False
        for (j) in (range(i)):
            coefficient = h1[i] * r0[j] + h1[j] * r0[i] - h0[i] * r1[j] - h0[j] * r1[i]
            coefficient -= value * (r1[i] * r0[j] + r1[j] * r0[i])
            if (coefficient != 0):
                return False
    return True

def endpoint(r0, h0, r1, h1, kind):
    arrays = [[rational(value) for (value) in (array)] for (array) in ([r0, h0, r1, h1])]
    r0, h0, r1, h1 = arrays
    count = len(r0)
    if ((count < 2) or any((len(array) != count) for (array) in (arrays))):
        raise ValueError("Invalid dimension")
    if (any(
        ((r < 0) or (h < 0) or (h > r))
        for (returned, positive) in ([(r0, h0), (r1, h1)])
        for (r, h) in (zip(returned, positive))
    )):
        raise ValueError("Invalid binary masses")
    means = [sum(array, Fraction(0)) / count for (array) in (arrays)]
    m0, b0, m1, b1 = means
    if (kind == "delta_R"):
        value = m1 - m0
        influences = [new - old - value for (old, new) in (zip(r0, r1))]
        constant = all((influence == 0) for (influence) in (influences))
    elif (kind == "delta_p"):
        if ((m0 == 0) or (m1 == 0)):
            return {"status": "UNDEFINED_ORIGINAL_DENOMINATOR", "value": None}
        p0 = b0 / m0
        p1 = b1 / m1
        value = p1 - p0
        influences = [
            (hn - p1 * rn) / m1 - (ho - p0 * ro) / m0 for (ro, ho, rn, hn) in (zip(*arrays))
        ]
        constant = ratio_constant(r0, h0, r1, h1, value)
    else:
        raise ValueError("Unknown endpoint")
    if (sum(influences, Fraction(0)) != 0):
        raise ValueError("Influence centring")
    variance = sum((value * value for (value) in (influences)), Fraction(0)) / (count * (count - 1))
    if ((variance == 0) and not constant):
        return {
            "status": "UNAVAILABLE_DEGENERATE_LINEARISATION",
            "value": str(value),
            "value_float": float(value),
            "scale_squared": "0",
            "scale": 0.0,
            "influences": [str(item) for (item) in (influences)],
            "constant_certificate": False,
        }
    if (constant and (variance != 0)):
        raise ValueError("Constant with nonzero influence")
    scale = float(variance) ** 0.5
    if (not np.isfinite(scale) or ((variance > 0) and (scale == 0))):
        raise ValueError("Invalid binary64 scale")
    return {
        "status": "CONSTANT" if (constant) else "AVAILABLE",
        "value": str(value),
        "value_float": float(value),
        "scale_squared": str(variance),
        "scale": scale,
        "influences": [str(value) for (value) in (influences)],
        "constant_certificate": constant,
    }

def count_draws(draws, count):
    draws = np.asarray(draws)
    if (
        (draws.ndim != 2)
        or (draws.shape[1] != count)
        or not np.issubdtype(draws.dtype, np.integer)
    ):
        raise ValueError("Invalid donor draws")
    if (np.any(draws < 0) or np.any(draws >= count)):
        raise ValueError("Out-of-range donor index")
    counts = np.zeros((len(draws), count), dtype = np.int64)
    np.add.at(counts, (np.arange(len(draws))[:, None], draws), 1)
    return counts

def resample_values(counts, arrays, kind):
    totals = (
        counts
        @ np.asarray(
            [[float(value) for (value) in (array)] for (array) in (arrays)], dtype = np.float64
        ).T
    )
    if (not np.all(np.isfinite(totals))):
        raise ValueError("Nonfinite draw aggregate")
    undefined = np.zeros(len(counts), dtype = bool)
    if (kind == "delta_R"):
        values = (totals[:, 2] - totals[:, 0]) / counts.sum(axis = 1)
    else:
        undefined = (totals[:, 0] == 0) | (totals[:, 2] == 0)
        values = np.full(len(counts), np.nan)
        valid = ~undefined
        values[valid] = (
            totals[valid, 3] / totals[valid, 2] - totals[valid, 1] / totals[valid, 0]
        )
    return values, undefined

def calibrate(records, values):
    if (any((row["status"] == "UNDEFINED_ORIGINAL_DENOMINATOR") for (row) in (records))):
        return {"status": "UNAVAILABLE_ORIGINAL_DENOMINATOR", "q": None}
    if (any(row["status"] == "UNAVAILABLE_DEGENERATE_LINEARISATION" for (row) in (records))):
        return {"status": "UNAVAILABLE_DEGENERATE_LINEARISATION", "q": None}
    if (not np.all(np.isfinite(values))):
        return {"status": "UNAVAILABLE_UNDEFINED_DRAW", "q": None}
    varying = [index for (index, row) in (enumerate(records)) if (row["status"] != "CONSTANT")]
    if (varying):
        centres = np.asarray([records[index]["value_float"] for (index) in (varying)])
        scales = np.asarray([records[index]["scale"] for (index) in (varying)])
        maxima = np.max(np.abs(values[:, varying] - centres) / scales, axis = 1)
    else:
        maxima = np.zeros(len(values))
    quantile = float(np.quantile(maxima, 0.95, method = "linear"))
    return {
        "status": "AVAILABLE",
        "q": quantile,
        "maxima": maxima,
        "constant_coordinates": len(records) - len(varying),
    }

def prepare_levels(model_levels, context):
    aliases = {"DawnRank_strict": "DawnRank"}
    policies = (
        ("PRODIGY", "DawnRank", "degree") if (context == "CRC") else ("DawnRank", "degree")
    )
    mapping = {
        "returned": "R",
        "positive": "H",
        "zero": "Z",
        "unknown": "U",
        "absent": "absent",
        "null": "explicit_null",
        "identity_conflict": "identity_conflict",
        "call_conflict": "call_conflict",
    }
    rows, model_roster = [], {}
    for (source) in (model_levels):
        require((source["context"] == context), "Mixed category contexts")
        policy = aliases.get(source["policy"], source["policy"])
        states = {"G0": "0", "G1": "1"} if (policy == "degree") else {"G0H0": "0", "G1H1": "1"}
        if (source["state"] not in states):
            continue
        require((policy in policies), "Unexpected category policy")
        category = source["category"]
        require(
            (category in {item for (group) in (CATEGORIES.values()) for (item) in (group)}),
            "Unknown category",
        )
        mid = source["model_id"]
        identity = (source["sample_id"], source["donor_id"])
        require(
            (mid not in model_roster or model_roster[mid] == identity),
            "Model identity mismatch",
        )
        model_roster[mid] = identity
        mass = {target: rational(source[field]) for (field, target) in (mapping.items())}
        validate_mass(mass)
        require(
            (category not in FOCAL or mass["U"] == 0),
            "Focal unknown mass is outside this diagnostic",
        )
        rows.append(
            {
                "context": context,
                "policy": policy,
                "category": category,
                "state": states[source["state"]],
                "model_id": mid,
                "sample_id": source["sample_id"],
                "donor_id": source["donor_id"],
                **mass,
            }
        )
    expected_models, expected_donors = POPULATIONS[context]
    require(
        (
            len(model_roster) == expected_models
            and len({value[0] for (value) in (model_roster.values())}) == expected_models
        ),
        "Full model population required",
    )
    require(
        (len({value[1] for (value) in (model_roster.values())}) == expected_donors),
        "Full donor population required",
    )
    indexed = {
        (row["model_id"], row["policy"], row["category"], row["state"]): row for (row) in (rows)
    }
    categories = [item for (group) in (CATEGORIES.values()) for (item) in (group)]
    require(
        (len(indexed) == len(rows) == expected_models * len(policies) * len(categories) * 2),
        "Complete unique model category grid required",
    )
    for (mid) in (model_roster):
        for (policy) in (policies):
            for (state) in (("0", "1")):
                grouped = {
                    category: indexed[(mid, policy, category, state)] for (category) in (categories)
                }
                for (partition) in (("core", "role")):
                    for (field) in (FIELDS):
                        require(
                            (
                                sum(
                                    (
                                        grouped[category][field]
                                        for (category) in (CATEGORIES[partition])
                                    ),
                                    Fraction(),
                                )
                                == grouped["ALL"][field]
                            ),
                            "Category partition conservation",
                        )
    return rows, aggregate_donors(rows)

def family(model_levels, donor_ids, draws, context):
    model_rows, donors = prepare_levels(model_levels, context)
    policies = (
        ("PRODIGY", "DawnRank", "degree") if (context == "CRC") else ("DawnRank", "degree")
    )
    require(
        (donor_ids == sorted(donor_ids) and len(set(donor_ids)) == POPULATIONS[context][1]),
        "Full ordered donor axis required",
    )
    require(
        (set(donor_ids) == {row["donor_id"] for (row) in (donors)}), "Draw donor identity mismatch"
    )
    require((len(draws) == 10000), "All 10000 paired donor draws required")
    counts = count_draws(draws, len(donor_ids))
    indexed = {
        (row["policy"], row["category"], row["state"], row["donor_id"]): row for (row) in (donors)
    }
    records, values, levels, undefined_rows = [], [], [], []
    for (policy) in (policies):
        for (category) in (FOCAL):
            old = [indexed[(policy, category, "0", donor)] for (donor) in (donor_ids)]
            new = [indexed[(policy, category, "1", donor)] for (donor) in (donor_ids)]
            arrays = [
                [row[field] for (row) in (rows)]
                for (rows, field) in (((old, "R"), (old, "H"), (new, "R"), (new, "H")))
            ]
            for (state, rows) in ((("0", old), ("1", new))):
                mass = {
                    field: sum((row[field] for (row) in (rows)), Fraction()) / len(rows)
                    for (field) in (FIELDS)
                }
                levels.append(
                    {
                        "context": context,
                        "policy": policy,
                        "category": category,
                        "state": state,
                        **mass,
                        "pooled_positive_fraction": mass["H"] / mass["R"]
                        if (mass["R"] != 0)
                        else None,
                    }
                )
            for (kind) in (("delta_R", "delta_p")):
                record = endpoint(*arrays, kind)
                identifier = "__".join((policy, category, kind))
                record.update(
                    {"id": identifier, "policy": policy, "category": category, "endpoint": kind}
                )
                sampled, undefined = resample_values(counts, arrays, kind)
                for (index) in (np.flatnonzero(undefined)):
                    undefined_rows.append(
                        {
                            "context": context,
                            "endpoint_id": identifier,
                            "draw_index": int(index),
                            "donor_indices": draws[index],
                        }
                    )
                records.append(record)
                values.append(sampled)
    matrix = np.asarray(values).T
    calibrated = calibrate(records, matrix)
    maxima = calibrated.pop("maxima", None)
    for (record) in (records):
        q = calibrated["q"]
        if (q is not None):
            record["lower"] = record["value_float"] - q * record["scale"]
            record["upper"] = record["value_float"] + q * record["scale"]
            record["direction"] = (
                "NEGATIVE"
                if (record["upper"] < 0)
                else "POSITIVE"
                if (record["lower"] > 0)
                else "ZERO"
                if (record["lower"] == record["upper"] == 0)
                else "CROSSES_ZERO"
            )
        else:
            record.update({"lower": None, "upper": None, "direction": "UNAVAILABLE"})
    result = {
        "context": context,
        "donor_ids": donor_ids,
        "models": POPULATIONS[context][0],
        "draws": len(draws),
        "family_size": len(records),
        **calibrated,
        "endpoints": records,
        "undefined_draw_endpoint_count": len(undefined_rows),
    }
    samples = [
        {
            "draw_index": index,
            **{
                record["id"]: float(matrix[index, column])
                if (np.isfinite(matrix[index, column]))
                else None
                for (column, record) in (enumerate(records))
            },
            "maximum_standardised_deviation": float(maxima[index])
            if (maxima is not None)
            else None,
        }
        for (index) in (range(len(draws)))
    ]
    return result, {
        "model_levels": model_rows,
        "donor_levels": donors,
        "original_levels": levels,
        "resamples": samples,
        "undefined_draws": undefined_rows,
    }
