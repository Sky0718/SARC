from pathlib import Path

from . import category, category_sensitivity, endpoint, paired, score
from .draws import colorectal_draws
from .inputs import load_categories, load_crc, load_draws, load_esca
from .inputs import read_csv, read_json, resolve, save_csv, save_json

def donor_order(path):
    document = read_json(path)
    if (isinstance(document, dict)):
        if ("models" in document):
            return sorted({row["donor_id"] for (row) in (document["models"])})
        document = document["donor_ids"] if ("donor_ids" in document) else document["donors"]
    return [row["donor_id"] if (isinstance(row, dict)) else row for (row) in (document)]

def paired_updates(base, config, output):
    models, allocations, means, degrees, calls, policies = load_crc(base, config)
    result = paired.analyse(models, allocations, means, degrees, calls)
    save_json(output / "paired_updates.json", result)
    rows = [
        {"comparison": name, "unit": unit[:-1], **row}
        for (name, reduced) in (result["update_comparator"].items())
        for (unit) in (("models", "donors"))
        for (row) in (reduced[unit])
    ]
    rows.extend(
        {"comparison": name, "unit": "population", **reduced["population"]}
        for (name, reduced) in (result["update_comparator"].items())
    )
    save_csv(output / "update_comparator.csv", rows)
    save_csv(output / "seed_opposition.csv", result["seed_rows"])
    save_csv(output / "evidence_zero.csv", result["evidence_rows"])

def endpoint_sensitivity(base, config, output):
    result = read_json(resolve(base, config["paired_updates"]))
    donors = donor_order(resolve(base, config["donor_order"]))
    draws = load_draws(resolve(base, config["draws"]), donors)
    summary, vectors, samples = endpoint.colorectal_family(result, draws, donors)
    save_json(output / "endpoint_sensitivity.json", summary)
    save_json(output / "donor_vectors.json", vectors)
    save_json(output / "resampled_endpoints.json", samples)

def category_summary(base, config, output):
    context = config["context"]
    if (context == "CRC"):
        models, allocations, means, degrees, calls, policies = load_crc(base, config)
    elif (context == "ESCA"):
        models, policies, calls = load_esca(base, config)
    else:
        raise ValueError("Unknown category context")
    categories = load_categories(resolve(base, config["categories"]), models, context)
    tables = category.derive(models, policies, calls, categories, context)
    for (name, rows) in (tables.items()):
        if (rows):
            save_csv(output / (name + ".csv"), rows)
    save_json(
        output / "category_summary.json",
        {
            "context": context,
            "models": len(models),
            "donors": len({row["donor_id"] for (row) in (models)}),
            "population_levels": tables["population_category_levels"],
            "population_contrasts": tables["population_category_contrasts"],
            "oncogene_concentration": tables["oncogene_concentration"],
        },
    )

def category_envelopes(base, config, output):
    models = read_csv(resolve(base, config["model_levels"]))
    donors = donor_order(resolve(base, config["donor_order"]))
    draws = load_draws(resolve(base, config["draws"]), donors)
    result, tables = category_sensitivity.family(models, donors, draws, config["context"])
    save_json(output / "category_sensitivity.json", result)
    for (name, rows) in (tables.items()):
        if (rows):
            save_csv(output / (name + ".csv"), rows)
        elif (name == "undefined_draws"):
            save_csv(
                output / "undefined_draws.csv",
                [],
                fields = ["context", "endpoint_id", "draw_index", "donor_indices"],
            )

def signed_accounting(base, config, output):
    rows = read_csv(resolve(base, config["donor_allocation"]))
    result, vectors = score.summarise(rows)
    save_json(output / "signed_accounting.json", result)
    save_csv(output / "signed_mass_summary.csv", result["summary"])
    save_csv(output / "signed_donor_vectors.csv", vectors)

def generate_draws(base, config, output):
    donors = donor_order(resolve(base, config["donor_order"]))
    save_json(output / "donor_order.json", {"donor_ids": donors})
    save_json(output / "draw_indices.json", colorectal_draws(donors))

OPERATIONS = {
    "paired-updates": paired_updates,
    "endpoint-sensitivity": endpoint_sensitivity,
    "category-summary": category_summary,
    "category-sensitivity": category_envelopes,
    "signed-accounting": signed_accounting,
    "colorectal-draws": generate_draws,
}

def run(config_path):
    path = Path(config_path).resolve()
    config = read_json(path)
    operation = config["operation"]
    if (operation not in OPERATIONS):
        raise ValueError("Unknown statistics operation")
    output = resolve(path.parent, config["output"])
    output.mkdir(parents = True, exist_ok = False)
    OPERATIONS[operation](path.parent, config, output)
