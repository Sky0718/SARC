from pathlib import Path
from tempfile import mkdtemp
import json
import re
import sys
import duckdb
from .pipeline import (
    acquisition,
    broader,
    harmonisation,
    ranking,
    robustness,
    score_revision,
    sensitivity,
    supports,
)

SOURCE_ROOT = Path(__file__).resolve().parents[2]
RESOURCE_ROOT = (
    SOURCE_ROOT
    if ((SOURCE_ROOT / "config" / "support_construction.json").is_file())
    else Path(sys.prefix) / "share" / "sarc"
)
CONFIGURATION_ROOT = RESOURCE_ROOT / "config"
DOCUMENT_ROOT = RESOURCE_ROOT / "docs"

def resources(threads, memory_limit):
    if (isinstance(threads, bool) or int(threads) < 1):
        raise ValueError("Thread count must be positive")
    value = str(memory_limit).upper()
    if (not re.fullmatch(r"[1-9][0-9]*(KB|MB|GB|TB)", value)):
        raise ValueError("Memory limit must include KB, MB, GB or TB")
    return int(threads), value

def configuration(name, config_dir, threads, memory_limit):
    path = Path(config_dir or CONFIGURATION_ROOT) / f"{name}.json"
    value = json.loads(path.read_text(encoding = "utf-8"))
    value.setdefault("resources", {}).update(threads = threads, memory_limit = memory_limit)
    return value

def work_path(root, relative):
    path = (root / relative).resolve()
    if (not path.is_relative_to(root)):
        raise ValueError(f"Path escapes the working directory: {relative}")
    return path

def input_paths(root, config):
    paths = {}
    for (name, relative) in (config["inputs"].items()):
        path = work_path(root, relative)
        if (name == "release_change_matrix" and not path.is_file()):
            path = DOCUMENT_ROOT / "release_change_matrix.csv"
        if (not path.is_file()):
            raise FileNotFoundError(path)
        paths[name] = path
    return paths

def output_paths(root, config, names = None):
    keys = list(config["outputs"]) if (names is None) else list(names)
    paths = {name: work_path(root, config["outputs"][name]) for (name) in (keys)}
    existing = [str(path) for (path) in (paths.values()) if (path.exists())]
    if (existing):
        raise FileExistsError(
            f"Use a new working directory or remove the selected existing outputs: {existing}"
        )
    for (path) in (paths.values()):
        path.parent.mkdir(parents = True, exist_ok = True)
    return paths

def temporary_directory(root, stage):
    parent = work_path(root, "data/interim")
    parent.mkdir(parents = True, exist_ok = True)
    return Path(mkdtemp(prefix = f"{stage}_", dir = parent))

def download(root, config_dir = None, threads = 8, memory_limit = "8GB", workers = 8):
    if (workers < 1):
        raise ValueError("Download concurrency must be positive")
    return acquisition.run(Path(root).resolve(), workers)

def harmonise(root, config_dir = None, threads = 8, memory_limit = "8GB", workers = 8):
    root = Path(root).resolve()
    threads, memory_limit = resources(threads, memory_limit)
    config = configuration("support_construction", config_dir, threads, memory_limit)
    paths = output_paths(
        root,
        config,
        ["disease_crosswalk", "target_crosswalk", "therapeutic_area_mapping"],
    )
    inputs = harmonisation.discover_inputs(root, config)
    temporary = temporary_directory(root, "harmonisation")
    connection = duckdb.connect()
    try:
        connection.execute(f"SET threads = {threads}")
        connection.execute(
            f"SET memory_limit = {harmonisation.quote_literal(memory_limit)}"
        )
        connection.execute(
            f"SET temp_directory = {harmonisation.quote_literal(temporary.as_posix())}"
        )
        diseases = {
            release: harmonisation.load_disease_records(
                connection, inputs[(release, "disease")]
            )
            for (release) in (config["releases"])
        }
        targets = {
            release: harmonisation.load_target_records(
                connection, inputs[(release, "target")]
            )
            for (release) in (config["releases"])
        }
        disease_rows, area_rows, _, _ = harmonisation.build_disease_crosswalk(
            diseases, config
        )
        target_rows = harmonisation.build_target_crosswalk(targets, config)
        if (not disease_rows or not target_rows):
            raise ValueError("Entity crosswalks are empty")
        harmonisation.write_csv(
            paths["disease_crosswalk"], disease_rows, list(disease_rows[0])
        )
        harmonisation.write_csv(
            paths["target_crosswalk"], target_rows, list(target_rows[0])
        )
        area_fields = [
            "release",
            "source_disease_id",
            "canonical_disease_id",
            "source_therapeutic_area_id",
            "canonical_therapeutic_area_id",
            "canonical_therapeutic_area_label",
            "therapeutic_area_stability",
            "primary_area_contrast_eligible",
            "fractional_weight",
        ]
        harmonisation.write_csv(
            paths["therapeutic_area_mapping"], area_rows, area_fields
        )
    finally:
        connection.close()
    return paths

def build_supports(root, config_dir = None, threads = 8, memory_limit = "8GB", workers = 8):
    root = Path(root).resolve()
    threads, memory_limit = resources(threads, memory_limit)
    config = configuration("support_construction", config_dir, threads, memory_limit)
    inherited = {"disease_crosswalk", "target_crosswalk", "therapeutic_area_mapping"}
    output_paths(root, config, set(config["outputs"]) - inherited)
    for (name) in (inherited):
        if (not work_path(root, config["outputs"][name]).is_file()):
            raise FileNotFoundError(work_path(root, config["outputs"][name]))
    supports.run(root, config, threads, memory_limit)

def scores(root, config_dir = None, threads = 8, memory_limit = "8GB", workers = 8):
    root = Path(root).resolve()
    threads, memory_limit = resources(threads, memory_limit)
    config = configuration("score_revision", config_dir, threads, memory_limit)
    paths = output_paths(root, config)
    inputs = input_paths(root, config)
    temporary = temporary_directory(root, "scores")
    connection = duckdb.connect(str(temporary / "analysis.duckdb"))
    try:
        score_revision.configure_duckdb_resources(connection, temporary, config)
        score_revision.create_input_views(connection, inputs, config)
        score_revision.create_analysis_tables(connection, config)
        score_revision.export_stage(connection, root / "data/derived", config)
    finally:
        connection.close()
    return paths

def ranks(root, config_dir = None, threads = 8, memory_limit = "8GB", workers = 8):
    root = Path(root).resolve()
    threads, memory_limit = resources(threads, memory_limit)
    config = configuration("ranking_sensitivity", config_dir, threads, memory_limit)
    paths = output_paths(root, config)
    inputs = input_paths(root, config)
    temporary = temporary_directory(root, "ranks")
    connection = duckdb.connect(str(temporary / "analysis.duckdb"))
    try:
        ranking.configure_duckdb_resources(connection, temporary, config)
        ranking.create_input_views(connection, inputs, config)
        datasets = ranking.build_ranking_sensitivity_datasets(connection, config)
        ranking.export_stage(connection, root / "data/derived", datasets, config)
    finally:
        connection.close()
    return paths

def cross_layer(root, config_dir = None, threads = 8, memory_limit = "8GB", workers = 8):
    root = Path(root).resolve()
    threads, memory_limit = resources(threads, memory_limit)
    config = configuration("cross_layer_sensitivity", config_dir, threads, memory_limit)
    paths = output_paths(root, config)
    inputs = input_paths(root, config)
    diagnostics = sensitivity.discover_diagnostic_inputs(root, config)
    temporary = temporary_directory(root, "cross_layer")
    connection = duckdb.connect(str(temporary / "analysis.duckdb"))
    try:
        sensitivity.configure_connection(connection, root, config, temporary.name)
        sensitivity.create_input_views(connection, inputs, config)
        datasets = sensitivity.build_cross_layer_sensitivity_datasets(
            connection, root, inputs, config, diagnostics
        )
        sensitivity.export_stage(connection, root / "data/derived", datasets, config)
    finally:
        connection.close()
    return paths

def assess_robustness(root, config_dir = None, threads = 8, memory_limit = "8GB", workers = 8):
    resources(threads, memory_limit)
    return robustness.run(
        Path(root).resolve(), threads = threads, memory_limit = memory_limit
    )

def broader_efo(root, config_dir = None, threads = 8, memory_limit = "8GB", workers = 8):
    threads, memory_limit = resources(threads, memory_limit)
    return broader.run(Path(root).resolve(), threads = threads, memory_limit = memory_limit)

def analyse(root, config_dir = None, threads = 8, memory_limit = "8GB", workers = 8):
    for (operation) in ((
        harmonise,
        build_supports,
        scores,
        ranks,
        cross_layer,
        assess_robustness,
        broader_efo,
    )):
        operation(root, config_dir, threads, memory_limit, workers)

COMMANDS = {
    "download": download,
    "harmonise": harmonise,
    "supports": build_supports,
    "scores": scores,
    "ranks": ranks,
    "sensitivity": cross_layer,
    "robustness": assess_robustness,
    "broader": broader_efo,
    "analyse": analyse,
}
