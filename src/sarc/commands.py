import argparse
from importlib import import_module
from pathlib import Path

COMMANDS = {
    "prepare": ("sarc.core.prepare", "prepare_common_domain"),
    "prepare-molecular": ("sarc.core.molecular_run", "prepare_molecular"),
    "compute": ("sarc.core.run", "run_job"),
    "export-native": ("sarc.core.export", "export_results"),
    "export-catalogue": ("sarc.core.catalogue_export", "export_catalogue"),
    "analyse": ("sarc.analysis.run", "run"),
    "evaluate": ("sarc.evaluation.run", "run"),
    "compare": ("sarc.resources.workflow", "compare_resources"),
    "compare-personadrive": ("sarc.resources.workflow", "compare_personadrive"),
    "benchmark": ("sarc.resources.workflow", "controlled_benchmark"),
    "benchmark-summary": ("sarc.resources.benchmark_summary", "summarise_benchmark"),
    "prepare-temporal": ("sarc.resources.temporal_workflow", "prepare_temporal"),
    "compare-temporal": ("sarc.commands", "compare_temporal"),
    "canonicalise-string": ("sarc.commands", "canonicalise_string"),
    "prepare-string": ("sarc.resources.string_workflow", "prepare_string"),
    "compare-string": ("sarc.resources.string_workflow", "compare_string"),
    "prepare-interventions": ("sarc.interventions.workflow", "prepare_graphs"),
    "prepare-controls": ("sarc.interventions.workflow", "prepare_matched_controls"),
    "assess-residuals": ("sarc.interventions.analysis", "assess_design"),
    "project-residuals": ("sarc.interventions.analysis", "project_residuals"),
    "background-contrasts": ("sarc.interventions.analysis", "background_contrasts"),
    "analyse-edits": ("sarc.interventions.endpoint_workflow", "run"),
    "integrate": ("sarc.integration.workflow", "run"),
    "catalogue": ("sarc.catalogue.run", "run"),
    "prediction-rounds": ("sarc.prediction_rounds.workflow", "run"),
    "select-portfolio": ("sarc.transfer.workflow", "select_portfolio"),
    "evaluate-transfer": ("sarc.transfer.workflow", "evaluate_transfer"),
    "transfer-intervals": ("sarc.transfer.workflow", "calculate_intervals"),
    "precision": ("sarc.precision.workflow", "run"),
    "verify-precision": ("sarc.precision.verification", "run"),
    "mechanism": ("sarc.mechanism.workflow", "run"),
    "oesophageal": ("sarc.oesophageal.workflow", "run"),
    "statistics": ("sarc.statistics.workflow", "run"),
    "diagnostics": ("sarc.diagnostics.workflow", "run"),
}

def compare_temporal(config_path):
    from .resources.io import read_json
    from .resources.temporal_workflow import compare_temporal as compare

    config = read_json(config_path)
    base = Path(config_path).resolve().parent
    compare(base / config["source"], base / config["output"])

def canonicalise_string(config_path):
    from .resources.io import read_json
    from .resources.string_workflow import canonicalise_links

    config = read_json(config_path)
    base = Path(config_path).resolve().parent
    canonicalise_links(base / config["source"], base / config["output"])

def main():
    parser = argparse.ArgumentParser(prog = "sarc.py")
    parser.add_argument("command", choices = tuple(COMMANDS))
    parser.add_argument("configuration", type = Path)
    arguments = parser.parse_args()
    module, function = COMMANDS[arguments.command]
    getattr(import_module(module), function)(arguments.configuration.resolve())
