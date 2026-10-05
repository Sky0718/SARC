from collections import defaultdict
from fractions import Fraction

DIRECTIONS = [
    "STRICTLY_NEGATIVE",
    "STRICTLY_POSITIVE",
    "IDENTIFIED_ZERO",
    "SIGN_UNRESOLVED",
]
MASS_KEYS = [
    "pro_negative",
    "dawn_negative",
    "pro_positive",
    "dawn_positive",
    "shared_negative",
    "shared_positive",
    "pro_negative_dawn_positive",
    "pro_positive_dawn_negative",
    "union_negative",
    "union_positive",
    "union_pro_negative_dawn_positive",
    "union_pro_positive_dawn_negative",
]
PARTS = [
    "total",
    "known_positive",
    "known_negative",
    "explicit_null",
    "absent",
    "lower",
    "upper",
]
RATIOS = {
    "negative": ("shared_negative", "union_negative"),
    "positive": ("shared_positive", "union_positive"),
    "pro_negative_dawn_positive": (
        "pro_negative_dawn_positive",
        "union_pro_negative_dawn_positive",
    ),
    "pro_positive_dawn_negative": (
        "pro_positive_dawn_negative",
        "union_pro_positive_dawn_negative",
    ),
}

def fraction(value):
    if (isinstance(value, bool) or not isinstance(value, (str, int, Fraction))):
        raise ValueError("Exact integer or rational string required")
    return Fraction(value)

def identity(row):
    value = (row["sample_id"], row["gene"])
    if (any((not isinstance(x, str) or not x for (x) in (value)))):
        raise ValueError("Invalid literal identity")
    return value

def labels_from(document):
    result = {}
    for (row) in (document["calls"]):
        key = identity(row)
        state, value = (row["state"], row["value"])
        if (state == "MEASURED"):
            if (type(value) is not int or value not in (0, 1)):
                raise ValueError("Invalid measured binary label")
        elif (state in ("EXPLICIT_NULL", "ABSENT")):
            if (value is not None):
                raise ValueError("Unknown label has value")
        else:
            raise ValueError("Invalid label state")
        if (key in result):
            raise ValueError("Duplicate selected label")
        result[key] = (state, value)
    return result

def merge_labels(pro, dawn):
    if (not pro.get("source_binding") or pro["source_binding"] != dawn.get(
        "source_binding"
    )):
        raise ValueError("The same assay source record is required for both methods")
    if (pro["label_contract"] != dawn["label_contract"]):
        raise ValueError("Incompatible label contracts")
    if (pro.get("data_class") != dawn.get("data_class")):
        raise ValueError("Incompatible data classes")
    result = labels_from(pro)
    for (key, value) in (labels_from(dawn).items()):
        if (key in result and result[key] != value):
            raise ValueError("Conflicting shared label identity")
        result[key] = value
    return result

def direction(lower, upper):
    if (lower > upper):
        raise ValueError("Invalid bounds")
    if (upper < 0):
        return DIRECTIONS[0]
    if (lower > 0):
        return DIRECTIONS[1]
    if (lower == upper == 0):
        return DIRECTIONS[2]
    return DIRECTIONS[3]

def sharp(vector, labels):
    known = Fraction(0)
    low = Fraction(0)
    high = Fraction(0)
    for (key, value) in (vector.items()):
        state, call = labels[key]
        if (state == "MEASURED"):
            known += value * call
        else:
            low += min(value, 0)
            high += max(value, 0)
    return {
        "known": str(known),
        "unknown_lower": str(low),
        "unknown_upper": str(high),
        "lower": str(known + low),
        "upper": str(known + high),
        "direction": direction(known + low, known + high),
    }

def mix(vectors, weights):
    if (len(vectors) != len(weights)):
        raise ValueError("Vector/weight length mismatch")
    result = defaultdict(Fraction)
    for (vector, weight) in (zip(vectors, weights)):
        for (key, value) in (vector.items()):
            result[key] += value * weight
    return dict(result)

def difference(pro, dawn):
    return {
        key: dawn.get(key, Fraction(0)) - pro.get(key, Fraction(0))
        for (key) in (pro.keys() | dawn.keys())
    }

def coefficient_vector(record, sample, labels):
    if (record.get("available") is not True or record.get("missing")):
        raise ValueError("Unavailable coefficient contrast")
    vector = {}
    for (row) in (record["coefficients"]):
        key = identity(row)
        if (key[0] != sample or key not in labels or key in vector):
            raise ValueError("Invalid or missing coefficient identity")
        vector[key] = fraction(row["coefficient"])
    return vector

def atom_masses(pro, dawn):
    pn, dn, pp, dp = (max(-pro, 0), max(-dawn, 0), max(pro, 0), max(dawn, 0))
    sn, sp, np, pn2 = (min(pn, dn), min(pp, dp), min(pn, dp), min(pp, dn))
    return dict(
        zip(
            MASS_KEYS,
            [
                pn,
                dn,
                pp,
                dp,
                sn,
                sp,
                np,
                pn2,
                pn + dn - sn,
                pp + dp - sp,
                pn + dp - np,
                pp + dn - pn2,
            ],
        )
    )

def zero_mass_summaries():
    return {name: {part: Fraction(0) for (part) in (PARTS)} for (name) in (MASS_KEYS)}

def add_atom(summaries, masses, label, weight = Fraction(1)):
    state, value = label
    for (name, mass) in (masses.items()):
        mass *= weight
        if (mass < 0):
            raise ValueError("Negative unsigned mass")
        target = summaries[name]
        target["total"] += mass
        if (state == "MEASURED"):
            target["known_positive" if (value == 1) else "known_negative"] += mass
        else:
            target["explicit_null" if (state == "EXPLICIT_NULL") else "absent"] += mass

def final_overlap(summaries):
    output = {}
    for (name, values) in (summaries.items()):
        values = dict(values)
        values["lower"] = values["known_positive"]
        values["upper"] = (
            values["known_positive"] + values["explicit_null"] + values["absent"]
        )
        if (values["total"] != sum(
            (
                values[part]
                for (part) in (
                    ["known_positive", "known_negative", "explicit_null", "absent"]
                )
            )
        )):
            raise ValueError("Functional partition mismatch")
        output[name] = {key: str(value) for ((key, value)) in (values.items())}
    ratios = {}
    for (name, (numerator, denominator)) in (RATIOS.items()):
        n, d = (summaries[numerator]["total"], summaries[denominator]["total"])
        ratios[name] = str(n / d) if (d) else None
    return {"masses": output, "weighted_jaccard": ratios}

def combine_mass_summaries(collection, weights):
    if (len(collection) != len(weights)):
        raise ValueError("Mass/weight length mismatch")
    result = zero_mass_summaries()
    for (source, weight) in (zip(collection, weights)):
        for (name) in (MASS_KEYS):
            for (part) in (PARTS):
                result[name][part] += source[name][part] * weight
    return result

def method_models(freeze, labels, protocol, pro):
    scope = protocol["scope"]
    if (freeze["budget"] != scope["budget"] or freeze["master_seeds"] != scope["seeds"]):
        raise ValueError("Budget or complete seed coverage mismatch")
    requests = [identity(row) for (row) in (freeze["requested_labels"])]
    if (len(requests) != len(set(requests)) or set(requests) != set(labels)):
        raise ValueError("Selected-call/request coverage mismatch")
    result = {}
    samples = set()
    for (model) in (freeze["models"]):
        mid, donor, sample = (model["model_id"], model["donor_id"], model["sample_id"])
        if (
            any((not isinstance(x, str) or not x for (x) in ([mid, donor, sample])))
            or mid in result
            or sample in samples
        ):
            raise ValueError("Duplicate or invalid model mapping")
        samples.add(sample)
        cs = {}
        for (name) in (protocol["contrasts"]):
            alias = protocol.get("pro_aliases", {}).get(name, name) if (pro) else name
            cs[name] = coefficient_vector(model["contrasts"][alias], sample, labels)
        result[mid] = {"donor_id": donor, "sample_id": sample, "contrasts": cs}
    if (
        len(result) != scope["models"]
        or len({x["donor_id"] for (x) in (result.values())}) != scope["donors"]
    ):
        raise ValueError("Full population coverage mismatch")
    if (any((sample not in samples for ((sample, gene)) in (requests)))):
        raise ValueError("Orphan selected sample")
    return result

def analysis(protocol, pro_freeze, dawn_freeze, pro_labels, dawn_labels):
    labels = merge_labels(pro_labels, dawn_labels)
    pro = method_models(pro_freeze, labels_from(pro_labels), protocol, True)
    dawn = method_models(dawn_freeze, labels_from(dawn_labels), protocol, False)
    if (set(pro) != set(dawn)):
        raise ValueError("Cross-method model population mismatch")
    donors = defaultdict(list)
    sample_genes = defaultdict(list)
    for (sample, gene) in (labels):
        sample_genes[sample].append(gene)
    for (mid) in (sorted(pro)):
        if (any(
            (pro[mid][key] != dawn[mid][key] for (key) in (["sample_id", "donor_id"]))
        )):
            raise ValueError("Cross-method donor/sample mapping mismatch")
        donors[pro[mid]["donor_id"]].append(mid)
    ndonors = len(donors)
    models = []
    rows = []
    masses_by_model = {}
    gene_data = {}
    for (mid) in (sorted(pro)):
        donor, sample = (pro[mid]["donor_id"], pro[mid]["sample_id"])
        weight = Fraction(1, ndonors * len(donors[donor]))
        item = {
            "model_id": mid,
            "donor_id": donor,
            "sample_id": sample,
            "population_weight": str(weight),
            "contrasts": {},
        }
        masses_by_model[mid] = {}
        for (contrast) in (protocol["contrasts"]):
            pv, dv = (pro[mid]["contrasts"][contrast], dawn[mid]["contrasts"][contrast])
            summaries = zero_mass_summaries()
            for (gene) in (sorted(sample_genes[sample])):
                key = (sample, gene)
                pc, dc = (pv.get(key, Fraction(0)), dv.get(key, Fraction(0)))
                masses = atom_masses(pc, dc)
                add_atom(summaries, masses, labels[key])
                state, value = labels[key]
                rows.append(
                    {
                        "model_id": mid,
                        "donor_id": donor,
                        "sample_id": sample,
                        "gene": gene,
                        "contrast": contrast,
                        "population_weight": str(weight),
                        "pro_coefficient": str(pc),
                        "dawn_coefficient": str(dc),
                        "paired_difference": str(dc - pc),
                        "state": state,
                        "value": value,
                        "masses": {
                            name: str(mass) for ((name, mass)) in (masses.items())
                        },
                    }
                )
                gkey = (contrast, gene)
                if (gkey not in gene_data):
                    gene_data[gkey] = {
                        "model_count": 0,
                        "nonzero_model_count": 0,
                        "pro_coefficient": Fraction(0),
                        "dawn_coefficient": Fraction(0),
                        "masses": zero_mass_summaries(),
                    }
                agg = gene_data[gkey]
                agg["model_count"] += 1
                agg["nonzero_model_count"] += int(pc != 0 or dc != 0)
                agg["pro_coefficient"] += pc * weight
                agg["dawn_coefficient"] += dc * weight
                add_atom(agg["masses"], masses, labels[key], weight)
            masses_by_model[mid][contrast] = summaries
            item["contrasts"][contrast] = {
                "pro": sharp(pv, labels),
                "dawn": sharp(dv, labels),
                "paired_difference": sharp(difference(pv, dv), labels),
                "overlap": final_overlap(summaries),
            }
        models.append(item)
    donor_results = []
    for (donor, mids) in (sorted(donors.items())):
        item = {"donor_id": donor, "model_ids": mids, "contrasts": {}}
        weights = [Fraction(1, len(mids))] * len(mids)
        for (contrast) in (protocol["contrasts"]):
            pv = mix([pro[mid]["contrasts"][contrast] for (mid) in (mids)], weights)
            dv = mix([dawn[mid]["contrasts"][contrast] for (mid) in (mids)], weights)
            masses = combine_mass_summaries(
                [masses_by_model[mid][contrast] for (mid) in (mids)], weights
            )
            item["contrasts"][contrast] = {
                "pro": sharp(pv, labels),
                "dawn": sharp(dv, labels),
                "paired_difference": sharp(difference(pv, dv), labels),
                "overlap": final_overlap(masses),
            }
        donor_results.append(item)
    population = {}
    contingency = {}
    mids = sorted(pro)
    weights = [
        Fraction(1, ndonors * len(donors[pro[mid]["donor_id"]])) for (mid) in (mids)
    ]
    for (contrast) in (protocol["contrasts"]):
        pv = mix([pro[mid]["contrasts"][contrast] for (mid) in (mids)], weights)
        dv = mix([dawn[mid]["contrasts"][contrast] for (mid) in (mids)], weights)
        masses = combine_mass_summaries(
            [masses_by_model[mid][contrast] for (mid) in (mids)], weights
        )
        population[contrast] = {
            "pro": sharp(pv, labels),
            "dawn": sharp(dv, labels),
            "paired_difference": sharp(difference(pv, dv), labels),
            "overlap": final_overlap(masses),
        }
        grid = []
        for (pdir) in (DIRECTIONS):
            for (ddir) in (DIRECTIONS):
                ids = [
                    donor["donor_id"]
                    for (donor) in (donor_results)
                    if (
                        donor["contrasts"][contrast]["pro"]["direction"] == pdir
                        and donor["contrasts"][contrast]["dawn"]["direction"] == ddir
                    )
                ]
                grid.append(
                    {
                        "pro_direction": pdir,
                        "dawn_direction": ddir,
                        "count": len(ids),
                        "donor_ids": ids,
                    }
                )
        contingency[contrast] = grid
    gene_aggregates = []
    for (contrast) in (protocol["contrasts"]):
        for ((c, gene), data) in (sorted(gene_data.items())):
            if (c != contrast):
                continue
            gene_aggregates.append(
                {
                    "contrast": contrast,
                    "gene": gene,
                    "model_count": data["model_count"],
                    "nonzero_model_count": data["nonzero_model_count"],
                    "pro_coefficient": str(data["pro_coefficient"]),
                    "dawn_coefficient": str(data["dawn_coefficient"]),
                    "paired_difference": str(
                        data["dawn_coefficient"] - data["pro_coefficient"]
                    ),
                    "overlap": final_overlap(data["masses"]),
                }
            )
    return {
        "schema": "cross_method_pairing_readout_v1",
        "counts": {
            "models": len(models),
            "donors": len(donor_results),
            "contrasts": len(protocol["contrasts"]),
            "selected_identities": len(labels),
            "sample_gene_rows": len(rows),
            "gene_aggregates": len(gene_aggregates),
        },
        "contrast_order": protocol["contrasts"],
        "direction_order": DIRECTIONS,
        "models": models,
        "donors": donor_results,
        "population": population,
        "donor_contingency": contingency,
        "sample_gene_rows": rows,
        "gene_aggregates": gene_aggregates,
    }
