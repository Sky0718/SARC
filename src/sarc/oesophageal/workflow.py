import os
import subprocess
from pathlib import Path
from .allocations import build_allocations
from .assay import select_rows, source_rows
from .contracts import freeze, new_directory, read, require, resolve, save, validate_models
from .counts import build_counts
from .functional import composition, evaluate_exact, serial, write_tables
from .mutation import build_mutations
from .scope import make_scope

def inputs(configuration, base, output):
    source_path = resolve(base, configuration["sources"])
    specification = read(source_path)
    source_base = source_path.parent
    models = specification["models"]
    validate_models(models)
    counts = build_counts(specification, models, source_base, output)
    mutation = build_mutations(
        specification["vcfs"],
        models,
        resolve(source_base, specification["gene_metadata"]),
        read(resolve(base, configuration["mutation_rule"])),
        source_base,
        output,
    )
    manifest = {
        "models": models,
        **counts,
        "tumour_counts": "tumour_htseq_counts.csv.gz",
        "normal_counts": "normal_rnaseqc_counts.csv.gz",
        "mutation_binary": mutation.name,
        "crosswalk": "gene_crosswalk.json",
    }
    save(output / "inputs.json", manifest)
    save(output / "models.json", {"models": models})
    return manifest

def prepare(configuration, base, output):
    input_path = resolve(base, configuration["inputs"])
    values = read(input_path)
    validate_models(values["models"])
    job = {
        "inputs": str(input_path),
        "output": str(output),
        **{
            name: str(resolve(input_path.parent, values[name]))
            for (name) in (("tumour_counts", "normal_counts", "mutation_binary"))
        },
        **{
            name: str(resolve(base, configuration[name]))
            for (name) in ((
                "common_genes",
                "native_genes_H0",
                "native_genes_H1",
                "edges_G0",
                "edges_G1",
                "normalise_source",
            ))
        },
    }
    job_path = output / "preparation_job.json"
    save(job_path, job)
    environment = os.environ.copy()
    if (configuration.get("r_libraries")):
        environment["R_LIBS_USER"] = os.pathsep.join(
            str(resolve(base, value)) for (value) in (configuration["r_libraries"])
        )
    for (variable) in (("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")):
        environment[variable] = "1"
    script = Path(__file__).with_name("entry.R")
    executable = configuration.get("rscript", "Rscript")
    if ("/" in executable or "\\" in executable):
        executable = str(resolve(base, executable))
    with (
        (output / "stdout.log").open("xb") as stdout,
        (output / "stderr.log").open("xb") as stderr,
    ):
        process = subprocess.run(
            [executable, "--vanilla", str(script), str(job_path)],
            stdout = stdout,
            stderr = stderr,
            env = environment,
            check = False,
        )
    require(process.returncode == 0, "Input preparation failed; inspect stderr.log")
    result = read(output / "preparation.json")
    require(result["status"] == "PREPARED", "Input preparation did not complete")
    return result

def allocations(configuration, base, output):
    models = read(resolve(base, configuration["inputs"]))["models"]
    return build_allocations(resolve(base, configuration["computation"]), models, output)

def scope(configuration, base, output):
    allocation = read(resolve(base, configuration["allocations"]))
    source_path = resolve(base, configuration["sources"])
    source = read(source_path)
    table = dict(source["table5_source"])
    table["path"] = str(resolve(source_path.parent, table["path"]))
    result = make_scope(
        allocation,
        resolve(base, configuration["crosswalk"]),
        resolve(source_path.parent, source["gene_metadata"]),
        resolve(source_path.parent, source["categories"]),
        table,
        source["excluded_sample_library_pairs"],
    )
    save(output / "scope.json", result)
    return result

def freeze_requests(configuration, base, output):
    allocation = read(resolve(base, configuration["allocations"]))
    specification = read(resolve(base, configuration["scope"]))
    frozen = freeze(allocation, specification)
    save(output / "coefficient_freeze.json", frozen)
    save(output / "requests.json", {"requests": frozen["requests"]})
    return {"model_count": len(allocation["models"]), "request_count": len(frozen["requests"])}

def evaluate(configuration, base, output):
    frozen = read(resolve(base, configuration["freeze"]))
    allocation = frozen["allocation"]
    specification = frozen["scope"]
    require(
        freeze(allocation, specification) == frozen, "Frozen coefficients or requests differ"
    )
    selected = select_rows(
        source_rows(specification["table5_source"]), frozen["requests"], specification
    )
    save(output / "selected_labels.json", selected)
    results = evaluate_exact(allocation, specification, frozen, selected)
    write_tables(output, results)
    uncertainty = composition(results["donor_contrasts"], frozen["donor_ids"], output)
    save(output / "composition.json", uncertainty)
    result = {
        "population_levels": results["population_levels"],
        "population_contrasts": results["population_contrasts"],
        "population_category_levels": results["population_category_levels"],
        "composition": uncertainty,
    }
    save(output / "results.json", serial(result))
    return result

def run(config_path):
    path = Path(config_path).resolve()
    configuration = read(path)
    operations = {
        "inputs": inputs,
        "prepare": prepare,
        "allocations": allocations,
        "scope": scope,
        "freeze": freeze_requests,
        "evaluate": evaluate,
    }
    operation = configuration["operation"]
    require(operation in operations, "Unknown oesophageal operation")
    output = new_directory(resolve(path.parent, configuration["output"]))
    return operations[operation](configuration, path.parent, output)
