import json
import re
from collections import Counter
from fractions import Fraction

SCALAR = re.compile("(?:null|-?(?:0|[1-9]\\d*)(?:\\.\\d+)?(?:[eE][+-]?\\d+)?)")

def require(condition, message):
    if (not condition):
        raise ValueError(message)

def fraction(value):
    require(
        type(value) is dict
        and set(value) == {"numerator", "denominator"}
        and (type(value["numerator"]) is int)
        and (type(value["denominator"]) is int)
        and (value["denominator"] > 0),
        "Exact accepted allocation fraction required",
    )
    result = Fraction(value["numerator"], value["denominator"])
    require(
        (result.numerator, result.denominator)
        == (value["numerator"], value["denominator"]),
        "Noncanonical accepted fraction",
    )
    return result

def validate_models(models, production = True):
    names = [r["sample_ID"] for (r) in (models)]
    require(
        len(names) == len(set(names))
        and all((isinstance(n, str) and n for (n) in (names))),
        "Original literal sample_ID roster required",
    )
    for (row) in (models):
        require(
            isinstance(row["study_donor_id"], str)
            and row["study_donor_id"]
            and (row["split"] in ("held_out", "development")),
            "Invalid donor/split metadata",
        )
    donors = {r["study_donor_id"] for (r) in (models)}
    require(
        all(
            (
                len({r["split"] for (r) in (models) if (r["study_donor_id"] == d)}) == 1
                for (d) in (donors)
            )
        ),
        "Related models cross original donor split",
    )
    if (production):
        expected = {"held_out": (59, 58), "development": (26, 25), "all": (85, 83)}
        for (subset, dimensions) in (expected.items()):
            rows = [r for (r) in (models) if (subset == "all" or r["split"] == subset)]
            require(
                (len(rows), len({r["study_donor_id"] for (r) in (rows)})) == dimensions,
                "Original subset population differs",
            )
    return {r["sample_ID"]: r for (r) in (models)}

def exact_weights(value):
    require(
        type(value) is dict and all((isinstance(g, str) and g for (g) in (value))),
        "Literal complete allocation menu required",
    )
    weights = {g: fraction(v) for ((g, v)) in (value.items())}
    require(
        all((0 <= v <= 1 for (v) in (weights.values())))
        and len({v for (v) in (weights.values()) if (0 < v < 1)}) <= 1,
        "Original fractional-cutoff allocation shape differs",
    )
    return weights

def pair_coefficients(native, prior):
    require(
        native is None and prior is None or (native is not None and prior is not None),
        "One-sided technical availability",
    )
    if (native is None):
        return None
    require(
        set(native) == set(prior)
        and all(
            (
                type(v) is Fraction and 0 <= v <= 1
                for (v) in (list(native.values()) + list(prior.values()))
            )
        ),
        "Same exact original menu required",
    )
    capacity = min(10, len(native))
    require(
        sum(native.values(), Fraction()) == sum(prior.values(), Fraction()) == capacity,
        "Native/F capacity differs",
    )
    coefficients = {
        g: native[g] - prior[g] for (g) in (sorted(native)) if (native[g] != prior[g])
    }
    require(
        sum(coefficients.values(), Fraction()) == 0,
        "Paired equal-capacity coefficients must cancel",
    )
    return coefficients

def carrier_record(raw, release, sample):
    require(
        raw["release"] == release and raw["sample_id"] == sample,
        "Accepted MODEL carrier identity differs",
    )
    if (raw["technical_status"] != "COMPLETE"):
        require(
            raw["technical_status"] == "TECHNICAL_UNAVAILABLE",
            "Unknown technical prediction state",
        )
        return {
            "release": release,
            "sample_ID": sample,
            "technical_status": "TECHNICAL_UNAVAILABLE",
            "coefficients": None,
            "native_F_equal": None,
            "menu": None,
            "returned_mass": None,
            "empty_menu": None,
            "all_returned_menu": None,
            "informative_menu": None,
            "native_T_equal_inherited": None,
            "D_native_T_inherited": None,
        }
    native, prior = (exact_weights(raw["W_native"]), exact_weights(raw["W_F"]))
    require(
        type(raw["menu"]) is list
        and len(raw["menu"]) == len(set(raw["menu"]))
        and (set(raw["menu"]) == set(native) == set(prior)),
        "Full offered menu differs from both allocations",
    )
    coefficients = pair_coefficients(native, prior)
    capacity = min(10, len(native))
    require(
        raw["offered_count"] == len(native)
        and raw["returned_mass"] == capacity
        and (raw["native_F_equal"] is (not coefficients)),
        "Inherited allocation capacity/equality metadata differs",
    )
    require(
        raw["empty_menu"] is (len(native) == 0)
        and raw["all_returned_menu"] is (len(native) <= 10)
        and (raw["informative_menu"] is (len(native) > 10)),
        "Inherited menu/vacuity flags differ",
    )
    t_distance = None if (raw["D_native_T"] is None) else fraction(raw["D_native_T"])
    require(
        (raw["native_T_equal"] is None) == (t_distance is None)
        and (t_distance is None or raw["native_T_equal"] is (t_distance == 0)),
        "Inherited native/T limitation metadata differs",
    )
    return {
        "release": release,
        "sample_ID": sample,
        "technical_status": "COMPLETE",
        "coefficients": coefficients,
        "native_F_equal": not coefficients,
        "menu": sorted(native),
        "returned_mass": capacity,
        "empty_menu": raw["empty_menu"],
        "all_returned_menu": raw["all_returned_menu"],
        "informative_menu": raw["informative_menu"],
        "native_T_equal_inherited": raw["native_T_equal"],
        "D_native_T_inherited": t_distance,
    }

def required_cells(records):
    return sorted(
        {
            (r["sample_ID"], gene)
            for (r) in (records)
            if (r["coefficients"] is not None)
            for (gene) in (r["coefficients"])
        }
    )

def selective_binary_calls(text, wanted, samples):
    wanted = set(wanted)
    require(
        wanted
        and all((m in samples and isinstance(g, str) and g for ((m, g)) in (wanted))),
        "Only frozen nonzero-cell label requests allowed",
    )
    decoder = json.JSONDecoder()
    position = 0
    model_seen, cells, scalar_count = (set(), {}, 0)

    def whitespace():
        nonlocal position
        while (position < len(text) and text[position] in " \t\r\n"):
            position += 1

    def expect(character):
        nonlocal position
        whitespace()
        require(
            position < len(text) and text[position] == character,
            "Malformed two-level accepted binary-cache structure",
        )
        position += 1

    def string_key():
        nonlocal position
        whitespace()
        require(
            position < len(text) and text[position] == '"',
            "Literal JSON object key required",
        )
        value, end = decoder.raw_decode(text, position)
        require(isinstance(value, str) and value, "Invalid model/gene key")
        position = end
        return value

    expect("{")
    whitespace()
    while (position < len(text) and text[position] != "}"):
        model = string_key()
        require(
            model in samples and model not in model_seen,
            "Unknown/duplicate original Table5 model key",
        )
        model_seen.add(model)
        expect(":")
        expect("{")
        genes_seen = set()
        whitespace()
        while (position < len(text) and text[position] != "}"):
            gene = string_key()
            require(gene not in genes_seen, "Duplicate original Table5 gene key")
            genes_seen.add(gene)
            expect(":")
            whitespace()
            match = SCALAR.match(text, position)
            require(match is not None, "Invalid scalar token in accepted binary cache")
            token = match.group()
            position = match.end()
            scalar_count += 1
            if ((model, gene) in wanted):
                value = json.loads(token)
                require(
                    value is None or (type(value) in (int, float) and value in (0, 1)),
                    "Selected accepted Table5 call is nonbinary",
                )
                cells[model, gene] = {
                    "mask": "EXPLICIT_NULL" if (value is None) else "PRESENT",
                    "call": None if (value is None) else int(value),
                }
            whitespace()
            require(
                position < len(text) and text[position] in ",}",
                "Invalid scalar termination",
            )
            if (text[position] == ","):
                position += 1
                whitespace()
                require(
                    position < len(text) and text[position] == '"',
                    "Trailing or missing gene key",
                )
            else:
                break
        expect("}")
        whitespace()
        require(
            position < len(text) and text[position] in ",}",
            "Invalid model-object termination",
        )
        if (text[position] == ","):
            position += 1
            whitespace()
            require(
                position < len(text) and text[position] == '"',
                "Trailing or missing model key",
            )
        else:
            break
    expect("}")
    whitespace()
    require(
        position == len(text) and model_seen == set(samples),
        "Cache has trailing data or missing original models",
    )
    for (key) in (wanted - set(cells)):
        cells[key] = {"mask": "KEY_ABSENT", "call": None}
    require(set(cells) == wanted, "Selective cache read changed requested cell set")
    return (
        cells,
        {
            "whole_cache_physical_decode_and_structural_scan": True,
            "scalar_tokens_structurally_scanned": scalar_count,
            "requested_distinct_cells": len(wanted),
            "selected_present_or_null_tokens_interpreted": sum(
                (r["mask"] != "KEY_ABSENT" for (r) in (cells.values()))
            ),
            "unrequested_scalar_tokens_numerically_interpreted": 0,
        },
    )

def value_bound(coefficients, labels):
    if (coefficients is None):
        return {
            "status": "TECHNICAL_UNAVAILABLE",
            "known_component": None,
            "lower": None,
            "upper": None,
            "point_identified": False,
            "label_states": None,
        }
    require(
        all((type(c) is Fraction and c != 0 for (c) in (coefficients.values())))
        and set(coefficients) <= set(labels),
        "Paired nonzero coefficients/selected labels required",
    )
    known, low_unknown, high_unknown = (Fraction(), Fraction(), Fraction())
    for (gene, coefficient) in (coefficients.items()):
        state = labels[gene]
        require(
            state["mask"] in ("PRESENT", "EXPLICIT_NULL", "KEY_ABSENT")
            and (
                state["mask"] == "PRESENT"
                and type(state["call"]) is int
                and (state["call"] in (0, 1))
                or (state["mask"] != "PRESENT" and state["call"] is None)
            ),
            "Selected binary-call mask/value inconsistent",
        )
        if (state["call"] is None):
            low_unknown += min(coefficient, Fraction())
            high_unknown += max(coefficient, Fraction())
        else:
            known += coefficient * state["call"]
    lower, upper = (known + low_unknown, known + high_unknown)
    status = (
        "NATIVE_STRICTLY_HIGHER"
        if (lower > 0)
        else "FIXED_F_STRICTLY_HIGHER"
        if (upper < 0)
        else "EXACT_ZERO"
        if (lower == upper == 0)
        else "SIGN_UNRESOLVED"
    )
    return {
        "status": status,
        "known_component": known,
        "lower": lower,
        "upper": upper,
        "point_identified": lower == upper,
        "label_states": {g: labels[g] for (g) in (coefficients)},
    }

def model_weights(models):
    counts = Counter((r["study_donor_id"] for (r) in (models)))
    weights = {
        r["sample_ID"]: Fraction(1, len(counts) * counts[r["study_donor_id"]])
        for (r) in (models)
    }
    require(
        len(weights) == len(models) and sum(weights.values(), Fraction()) == 1,
        "Full original donor-derived model weights required",
    )
    return weights

def summarise(records, models, release, subset, cells):
    weights = model_weights(models)
    require(
        {r["sample_ID"] for (r) in (records)} == set(weights)
        and len(records) == len(weights),
        "Complete fixed summary population required",
    )
    available = [r for (r) in (records) if (r["coefficients"] is not None)]
    unavailable = [r["sample_ID"] for (r) in (records) if (r["coefficients"] is None)]
    coefficients = {}
    for (record) in (available):
        sample = record["sample_ID"]
        for (gene, value) in (record["coefficients"].items()):
            key = (sample, gene)
            coefficients[key] = (
                coefficients.get(key, Fraction()) + weights[sample] * value
            )
    coefficients = {k: v for ((k, v)) in (coefficients.items()) if (v)}
    bound = value_bound(coefficients, cells)
    return {
        "release": release,
        "subset": subset,
        "primary": subset == "held_out",
        "models": len(models),
        "donors": len({r["study_donor_id"] for (r) in (models)}),
        "model_weights": weights,
        "unavailable_models": unavailable,
        "unavailable_model_mass": sum(
            (weights[m] for (m) in (unavailable)), Fraction()
        ),
        "full_population_status": "TECHNICAL_UNAVAILABLE"
        if (unavailable)
        else "COMPLETE",
        "available_original_weight_known_contribution": bound["known_component"],
        "available_original_weight_lower": bound["lower"],
        "available_original_weight_upper": bound["upper"],
        "mean_known_contribution": None if (unavailable) else bound["known_component"],
        "mean_lower": None if (unavailable) else bound["lower"],
        "mean_upper": None if (unavailable) else bound["upper"],
        "full_population_contrast_status": "TECHNICAL_UNAVAILABLE"
        if (unavailable)
        else bound["status"],
        "exact_equal_allocation_models": sum(
            (r["native_F_equal"] is True for (r) in (records))
        ),
        "empty_menu_models": sum((r["empty_menu"] is True for (r) in (records))),
        "all_returned_menu_models": sum(
            (r["all_returned_menu"] is True for (r) in (records))
        ),
        "informative_menu_models": sum(
            (r["informative_menu"] is True for (r) in (records))
        ),
        "inherited_native_T_different_models": sum(
            (r["native_T_equal_inherited"] is False for (r) in (records))
        ),
        "inherited_native_T_undefined_models": sum(
            (
                r["technical_status"] == "COMPLETE"
                and r["native_T_equal_inherited"] is None
                for (r) in (records)
            )
        ),
        "bounds_sharp_for_this_summary_only": not unavailable,
        "bounds_are_confidence_intervals": False,
    }
