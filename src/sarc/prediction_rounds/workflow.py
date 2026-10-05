import json
from pathlib import Path
from ..core.io import load_json, save_json, resolve_path
from ..evaluation.prediction import require, validate_rows
from .adapters import adapter, EXPECTED_TASKS, validate_partition_order
from .fitting import Fitter, safe_diagnostics
from .nested import make_manifest, fit_round
from .preparation import extract, join
from .report import summarise

def run(config_path):
    config_path = Path(config_path).resolve()
    config = load_json(config_path)
    base = config_path.parent
    output = Path(resolve_path(config["output"], base))
    require(
        not output.exists(),
        "A new output directory is required; existing results and failures are retained",
    )
    operation, round_id = config["operation"], config["round"]
    require(
        operation in ("features", "join", "fit", "summarise")
        and round_id in EXPECTED_TASKS,
        "Undeclared finite round operation",
    )
    output.mkdir(parents = True)
    save_json(output / "configuration.json", config)
    try:
        if (operation == "features"):
            result = extract(config, base)
            save_json(output / "features.json", result)
        else:
            rows = load_json(resolve_path(config["rows"], base))
            validate_rows(rows)
            if (operation == "join"):
                features = load_json(resolve_path(config["features"], base))
                result = join(rows, features, round_id)
                save_json(output / "rows.json", result)
            elif (operation == "fit"):
                specs = adapter(round_id)
                order = (
                    load_json(resolve_path(config["partition_order"], base))
                    if (round_id in ("R6", "R7"))
                    else None
                )
                if (order is not None):
                    validate_partition_order(rows, order)
                require(
                    round_id != "R7"
                    or isinstance(config.get("data_id"), str)
                    and bool(config["data_id"]),
                    "Explicit R7 carrier identity required",
                )
                manifest = make_manifest(rows, specs)
                require(
                    (len(manifest["inner_tasks"]), len(manifest["outer_tasks"]))
                    == EXPECTED_TASKS[round_id],
                    "Complete original finite grid required",
                )
                save_json(output / "manifest.json", manifest)
                result = fit_round(
                    rows, specs, Fitter(specs, order, config.get("data_id")), output
                )
                result["round"] = round_id
                save_json(output / "result.json", result)
            else:
                result = load_json(resolve_path(config["result"], base))
                prior = {
                    label: load_json(resolve_path(path, base))
                    for (label, path) in (config.get("prior_results", {}).items())
                }
                save_json(
                    output / "summary.json", summarise(rows, round_id, result, prior)
                )
        save_json(
            output / "completion.json",
            {
                "operation": operation,
                "round": round_id,
                "status": "COMPLETE",
                "adaptive_results_are_exploratory": True,
            },
        )
    except Exception as error:
        save_json(
            output / "failure.json",
            safe_diagnostics(
                {
                    "operation": operation,
                    "round": round_id,
                    "type": type(error).__name__,
                    "message": str(error),
                }
            ),
        )
        raise
    return str(output)
