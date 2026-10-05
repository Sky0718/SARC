import struct
from pathlib import Path

from .io import (
    load_json,
    read_matrix,
    resolve_path,
    save_json,
    unique_axis,
    write_matrix,
)

def equal_binary64(left, right):
    return all(
        struct.pack("<d", float(first)) == struct.pack("<d", float(second))
        for ((first, second)) in (zip(left, right))
    )

def restrict_evidence(record, genes):
    if (
        record["schema"] != "cdf_explicit_evidence_v1"
        or record["status"] != "AVAILABLE"
    ):
        raise ValueError("Explicit inherited evidence is not available")
    names = record["genes"]
    if (len(names) != len(set(names)) or any(not name for (name) in (names))):
        raise ValueError("Invalid evidence gene axis")
    encoded = record["binary64_le"]
    if (len(encoded) != 16 * len(names) or record["size"] != len(names)):
        raise ValueError("Evidence values and names differ in length")
    payload = bytes.fromhex(encoded)
    allowed = set(genes)
    indices = [index for ((index, gene)) in (enumerate(names)) if (gene in allowed)]
    result = dict(record)
    result["genes"] = [names[index] for (index) in (indices)]
    result["binary64_le"] = b"".join(
        payload[index * 8 : (index + 1) * 8] for (index) in (indices)
    ).hex()
    result["original_size"] = record.get("original_size", len(names))
    result["size"] = len(indices)
    result["non_null_empty"] = True
    result["restriction_preserves_original_order_and_bits"] = True
    return result

def prepare_common_domain(config_path):
    config_path = Path(config_path).resolve()
    config = load_json(config_path)
    base = config_path.parent
    output = Path(resolve_path(config["output"], base))
    counts = [read_matrix(resolve_path(path, base)) for (path) in (config["counts"])]
    mutations = [
        read_matrix(resolve_path(path, base)) for (path) in (config["mutations"])
    ]
    if (len(counts) != 2 or len(mutations) != 2):
        raise ValueError("Exactly two ordered releases are required")
    models = unique_axis(config["models"], "model")
    if (len(models) != config["expected_models"]):
        raise ValueError("Model population differs from the declared population")
    genes = [gene for (gene) in (counts[0][2]) if (gene in counts[1][2])]
    if (len(genes) != config["expected_genes"]):
        raise ValueError("Common-domain gene population differs from the declaration")
    if (counts[0][1] != counts[1][1] or any(item[1] != models for (item) in (mutations))):
        raise ValueError("Release sample order differs")
    for (gene) in (genes):
        if (gene not in mutations[0][2] or gene not in mutations[1][2]):
            raise ValueError("Common-domain mutation row is absent")
        if (not equal_binary64(counts[0][2][gene], counts[1][2][gene])):
            raise ValueError("Count values differ between releases")
        if (not equal_binary64(mutations[0][2][gene], mutations[1][2][gene])):
            raise ValueError("Mutation values differ between releases")
        if (any(float(value) not in (0.0, 1.0) for (value) in (mutations[0][2][gene]))):
            raise ValueError("Mutation indicators must be binary")
    output.mkdir(parents = True, exist_ok = False)
    write_matrix(output / "counts.csv", counts[0][0], counts[0][1], genes, counts[0][2])
    write_matrix(
        output / "mutation.csv", mutations[0][0], models, genes, mutations[0][2]
    )
    (output / "common_genes.txt").write_text("\n".join(genes) + "\n", encoding = "utf-8")
    (output / "models.txt").write_text("\n".join(models) + "\n", encoding = "utf-8")
    evidence_index = []
    for (index, item) in (enumerate(config.get("evidence", []))):
        source = load_json(resolve_path(item["path"], base))
        destination = "evidence_" + str(index + 1) + ".json"
        save_json(output / destination, restrict_evidence(source, genes))
        evidence_index.append(
            {key: value for ((key, value)) in (item.items()) if (key != "path")}
            | {"path": destination}
        )
    record = {
        "genes": len(genes),
        "models": models,
        "count_samples": counts[0][1],
        "evidence": evidence_index,
    }
    save_json(output / "prepared.json", record)
    return record
