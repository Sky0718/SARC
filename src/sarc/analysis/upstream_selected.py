import math
import struct
from collections import Counter
from . import hypothesis_analysis as hypothesis

def require(condition, message):
    if (not condition):
        raise ValueError(message)

def labels(values):
    require(isinstance(values, list), "Complete character labels required")
    result = [
        bytes.fromhex(value).decode("utf-8") if (isinstance(value, str)) else None
        for (value) in (values)
    ]
    require(
        None not in result and len(set(result)) == len(result),
        "Missing/duplicated labels",
    )
    return result

def fields(node):
    require(node["t"] == "list", "Named list required")
    names = labels(node["names"])
    require(len(names) == len(node["values"]), "Named list length mismatch")
    return dict(zip(names, node["values"]))

def scalar(node):
    if (node["t"] == "character"):
        values = labels(node["values"])
    else:
        require(node["t"] in ("integer", "double"), "Identity scalar type unavailable")
        raw = bytes.fromhex(node["bits"])
        values = [
            v[0]
            for (v) in (
                struct.iter_unpack("<i" if (node["t"] == "integer") else "<d", raw)
            )
        ]
        require(
            node.get("states", "F" * len(values)) == "F" * len(values),
            "Nonfinite identity scalar",
        )
    require(len(values) == 1, "Identity must be scalar")
    return values[0]

def selected(node, identity, universe):
    root = fields(node)
    saved = fields(root["identity"])
    require(
        all((scalar(saved[key]) == value for ((key, value)) in (identity.items()))),
        "Selected checkpoint identity mismatch",
    )
    vector = root["differences"]
    require(
        vector["t"] == "double" and vector["dim"] is None,
        "Full selected-DEG double vector unavailable",
    )
    raw = bytes.fromhex(vector["bits"])
    require(len(raw) % 8 == 0, "Selected vector byte length invalid")
    size = len(raw) // 8
    names = [] if (size == 0 and vector["names"] is None) else labels(vector["names"])
    require(all((name != "" for (name) in (names))), "Empty selected gene label")
    require(
        len(names) == size and vector["states"] == "F" * size,
        "Incomplete/nonfinite selected-DEG vector",
    )
    require(
        all(
            (
                math.isfinite(x[0]) and x[0] >= 0
                for (x) in (struct.iter_unpack("<d", raw))
            )
        ),
        "Invalid absolute-fold-change carrier",
    )
    require(
        set(names) <= universe,
        "Selected DEG outside verified projected count-input axis",
    )
    return set(names)

def compare(u0, u1, d0, d1):
    common = u0 & u1
    result = {
        "common_projected_input_genes": len(common),
        "old_only_projected_input_genes": len(u0 - u1),
        "new_only_projected_input_genes": len(u1 - u0),
    }
    if (d0 is None or d1 is None):
        return dict(
            result,
            status = "UNAVAILABLE_SELECTED_VECTOR",
            both_selected = None,
            old_selected_new_not_selected = None,
            old_not_selected_new_selected = None,
            neither_selected = None,
        )
    require(
        d0 <= u0 and d1 <= u1, "Selected set does not belong to complete input universe"
    )
    return dict(
        result,
        status = "AVAILABLE",
        both_selected = len(common & d0 & d1),
        old_selected_new_not_selected = len((common & d0) - d1),
        old_not_selected_new_selected = len((common & d1) - d0),
        neither_selected = len(common - (d0 | d1)),
    )

def comparison_key(sample, pair, old_origin, new_origin):
    return (sample, *pair, old_origin, new_origin)

def analyse(models, axes, origins, seed_states):
    samples = hypothesis.validate_roster(models)
    require(
        set(axes) == set(hypothesis.RELEASES),
        "Complete three-release expression axes required",
    )
    universes = {}
    first = axes[hypothesis.RELEASES[0]]
    for (release, axis) in (axes.items()):
        genes = axis["genes"]
        require(
            isinstance(genes, list)
            and len(genes) == len(set(genes))
            and all(isinstance(gene, str) and gene for (gene) in (genes)),
            "Complete unique gene labels required",
        )
        require(
            axis["samples"] == first["samples"]
            and axis["normal_samples"] == first["normal_samples"],
            "Fixed expression sample and normal-pool order differs",
        )
        require(
            axis["samples"] == axis["normal_samples"] + samples
            and len(axis["normal_samples"]) == 34,
            "Expected 34 fixed normals followed by 85 fixed tumours",
        )
        require(
            axis["sample_origins"] == ["normal"] * 34 + ["tumor"] * 85,
            "Expression sample roles differ",
        )
        require(
            isinstance(axis["source_id"], str) and axis["source_id"],
            "Explicit expression source identity required",
        )
        universes[release] = set(genes)
    origin_map = hypothesis.source_index(origins)
    states = hypothesis.state_index(models, seed_states, origin_map)
    cache, origin_records, contrasts, references = {}, [], {}, []
    for (identity, source) in (origin_map.items()):
        try:
            require(
                source["status"] == "AVAILABLE",
                source.get("reason", "Unavailable checkpoint source"),
            )
            cache[identity] = selected(
                source["checkpoint"],
                hypothesis.expected_identity(source),
                universes[source["release"]],
            )
            origin_records.append(
                {
                    "source_id": identity,
                    "status": "AVAILABLE",
                    "selected_gene_count": len(cache[identity]),
                    "checkpoint_source_seed": source["checkpoint_source_seed"],
                }
            )
        except (KeyError, ValueError, TypeError, struct.error) as error:
            cache[identity] = None
            origin_records.append(
                {
                    "source_id": identity,
                    "status": "UNAVAILABLE",
                    "reason": type(error).__name__ + ": " + str(error),
                    "checkpoint_source_seed": source["checkpoint_source_seed"],
                }
            )
    for (model) in (models):
        sample = model["sample_ID"]
        for (old_release, new_release) in (hypothesis.PAIRS):
            for (seed) in (hypothesis.SEEDS):
                old, new = (
                    states[sample, old_release, seed],
                    states[sample, new_release, seed],
                )
                old_id = old.get("checkpoint_source_id")
                new_id = new.get("checkpoint_source_id")
                old_origin = (
                    old_id
                    if (old_id is not None)
                    else "UNAVAILABLE_REFERENCE::"
                    + "|".join((sample, old_release, str(seed)))
                )
                new_origin = (
                    new_id
                    if (new_id is not None)
                    else "UNAVAILABLE_REFERENCE::"
                    + "|".join((sample, new_release, str(seed)))
                )
                key = comparison_key(
                    sample, (old_release, new_release), old_origin, new_origin
                )
                if (key not in contrasts):
                    old_selected = (
                        cache.get(old_id) if (old["status"] == "AVAILABLE") else None
                    )
                    new_selected = (
                        cache.get(new_id) if (new["status"] == "AVAILABLE") else None
                    )
                    contrasts[key] = dict(
                        compare(
                            universes[old_release],
                            universes[new_release],
                            old_selected,
                            new_selected,
                        ),
                        contrast_id = len(contrasts) + 1,
                        sample = sample,
                        model_id = model["model_id"],
                        donor = model["study_donor_id"],
                        release_pair = [old_release, new_release],
                        old_origin = old_origin,
                        new_origin = new_origin,
                        seed_references = [],
                    )
                expected_available = (
                    old["status"] == new["status"] == "AVAILABLE"
                    and cache.get(old_id) is not None
                    and cache.get(new_id) is not None
                )
                require(
                    (contrasts[key]["status"] == "AVAILABLE") == expected_available,
                    "Source-identical seed references have conflicting availability",
                )
                contrasts[key]["seed_references"].append(seed)
                references.append(
                    {
                        "sample": sample,
                        "release_pair": [old_release, new_release],
                        "seed": seed,
                        "old_mode": old.get("mode"),
                        "new_mode": new.get("mode"),
                        "old_checkpoint_source_seed": old.get("checkpoint_source_seed"),
                        "new_checkpoint_source_seed": new.get("checkpoint_source_seed"),
                        "contrast_id": contrasts[key]["contrast_id"],
                    }
                )
    require(
        len(references) == 1275, "Complete matched-seed comparison population required"
    )
    rows = list(contrasts.values())
    return {
        "origins": origin_records,
        "contrast_counts": rows,
        "all_seed_references": references,
        "counts": {
            "models": 85,
            "donors": 83,
            "releases": list(hypothesis.RELEASES),
            "pairs": [list(pair) for (pair) in (hypothesis.PAIRS)],
            "seed_references": len(references),
            "unique_checkpoint_origins": len(origins),
            "unique_source_pair_contrasts": len(rows),
            "origin_status_counts": dict(
                Counter(row["status"] for (row) in (origin_records))
            ),
            "contrast_status_counts": dict(Counter(row["status"] for (row) in (rows))),
        },
        "functional_labels_used": False,
        "interpretation": "Final upstream selected-gene status on the shared projected expression domain; no signed expression-effect or cause-specific attribution",
    }
