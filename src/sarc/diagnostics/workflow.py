from collections import defaultdict
from fractions import Fraction

from ..analysis.run import run_operation
from . import biology, degree, export, gene_groups, influence, inputs, menus, offered, rank
from .allocation import CELLS, mean_weights, require

def normalise(pro, dawn):
    models = inputs.roster(pro)
    second = {row["model_id"]: row for (row) in (inputs.roster(dawn))}
    require(set(second) == {row["model_id"] for (row) in (models)}, "Method rosters differ")
    allocations = inputs.rank_allocations(models, [pro, dawn])
    for (model) in (models):
        other = second[model["model_id"]]
        require(
            all(model[key] == other[key] for (key) in (("sample_id", "donor_id"))),
            "Method biological identities differ",
        )
        query = model.get("query", model.get("eligible_genes"))
        require(
            query is not None
            and set(query) == set(other.get("query", other.get("eligible_genes"))),
            "Method candidate queries differ",
        )
        model["eligible_genes"] = query
        model["methods"] = {
            method: {
                cell: mean_weights(
                    [
                        allocations[(method, model["model_id"], cell, seed, 10)]
                        for (seed) in (
                            inputs.SEEDS if (method == "PRODIGY") else ("deterministic",)
                        )
                    ]
                )
                for (cell) in (CELLS)
            }
            for (method) in (("PRODIGY", "DawnRank"))
        }
    return {"models": models}

def influence_summary(models, labels, methods = ("PRODIGY", "DawnRank")):
    models = inputs.roster(models)
    calls, origins = inputs.calls(labels)
    output = {"methods": {}, "tables": defaultdict(list)}
    for (method) in (methods):
        require(all(method in row["methods"] for (row) in (models)), "Method unavailable")
        summary, tables = influence.analyse_method(method, models, calls, origins)
        grouped, curve = gene_groups.group_diagnostic(
            tables["gene_influence"],
            *(summary["joint"][field] for (field) in (("known", "lower", "upper"))),
        )
        summary["group_gene_influence"] = grouped
        output["methods"][method] = summary
        for (name, rows) in (tables.items()):
            output["tables"][name].extend(rows)
        output["tables"]["group_gene_curves"].extend({"method": method, **row} for (row) in (curve))
    return output

def rank_summary(models, allocations, labels):
    models = inputs.roster(models)
    allocations = inputs.rank_allocations(models, allocations)
    calls, origins = inputs.calls(labels)
    inputs.cover_labels(models, allocations, calls)
    return rank.run(models, allocations, calls)

def degree_comparison(models, edges, domain, native, labels, budget = 10):
    models = inputs.roster(models)
    require(set(edges) == {"G0", "G1"}, "Two graph states required")
    require(budget == 10, "Degree comparison uses ten offered slots")
    degrees = {}
    graph_records = {}
    for (graph) in (("G0", "G1")):
        degrees[graph], graph_records[graph] = degree.degree_rows(edges[graph], domain)
    calls, origins = inputs.calls(labels)
    indexed = {}
    methods = set()
    for (row) in (native):
        method = row["method"]
        methods.add(method)
        cell = row["cell"].replace("_", "")
        seed = row.get("seed", "deterministic")
        key = (method, row["model_id"], cell, seed)
        require(key not in indexed, "Duplicate native menu/state")
        indexed[key] = row
    require(methods and methods <= {"PRODIGY", "DawnRank"}, "Unknown native methods")
    result = {"graphs": graph_records, "degree_allocations": {"models": []}, "comparisons": {}}
    for (model) in (models):
        query = model.get("query", model.get("eligible_genes"))
        require(query is not None, "Complete upstream query is required")
        full = {
            graph: degree.allocate(query, values, budget)[0]
            for (graph, values) in (degrees.items())
        }
        result["degree_allocations"]["models"].append(
            {key: model[key] for (key) in (("model_id", "sample_id", "donor_id", "alpha"))}
            | {"query": query, "full_degree": full}
        )
    full_by_model = {
        row["model_id"]: row["full_degree"] for (row) in (result["degree_allocations"]["models"])
    }
    consumed = set()
    for (method) in (sorted(methods)):
        seeds = inputs.SEEDS if (method == "PRODIGY") else ("deterministic",)
        for (cell) in (CELLS):
            for (graph) in (("G0", "G1")):
                contrasts = {
                    name: {}
                    for (name) in (("method_minus_full", "method_minus_menu", "menu_minus_full"))
                }
                for (model) in (models):
                    query = set(model.get("query", model.get("eligible_genes")))
                    actual, menu = [], []
                    for (seed) in (seeds):
                        key = (method, model["model_id"], cell, seed)
                        require(key in indexed, "Incomplete native matrix")
                        consumed.add(key)
                        source = indexed[key]
                        emitted = source["menu"]
                        require(len(emitted) == len(set(emitted)), "Duplicate emitted gene")
                        eligible = sorted(set(emitted) & query)
                        allocation = inputs.weights(source["weights"])
                        require(set(allocation) <= set(eligible), "Allocation outside menu")
                        require(
                            sum(allocation.values(), Fraction()) == min(budget, len(eligible)),
                            "Incomplete native capacity",
                        )
                        actual.append(allocation)
                        menu.append(degree.allocate(eligible, degrees[graph], budget)[0])
                    first, second = mean_weights(actual), mean_weights(menu)
                    third = full_by_model[model["model_id"]][graph]
                    vectors = {
                        "method_minus_full": biology.diff(first, third),
                        "method_minus_menu": biology.diff(first, second),
                        "menu_minus_full": biology.diff(second, third),
                    }
                    require(
                        biology.add(vectors["method_minus_menu"], vectors["menu_minus_full"])
                        == vectors["method_minus_full"],
                        "Menu coefficient identity differs",
                    )
                    for (name, vector) in (vectors.items()):
                        contrasts[name][model["model_id"]] = vector
                for (name, local) in (contrasts.items()):
                    summary, coefficients = biology.summarise_contrast(models, local, calls)
                    result["comparisons"]["_".join((method, cell, graph, name))] = summary
    require(consumed == set(indexed), "Unexpected native states or seeds")
    return result

def opportunities(records):
    from .opportunity import analyse

    prepared = []
    for (item) in (records):
        current = {"model": item["model"], "cutoff": item.get("cutoff", [])}
        for (method) in (("pro", "dawn")):
            current[method] = {}
            for (cell, value) in (item[method].items()):
                current[method][cell] = dict(value)
                current[method][cell]["weights"] = inputs.weights(value["weights"])
                current[method][cell]["menu"] = set(value["menu"])
        prepared.append(current)
    for (row, model) in (zip(prepared, inputs.roster([row["model"] for (row) in (prepared)]))):
        row["model"] = model
    return analyse(prepared)

def run(config_path):
    operations = {
        "prodigy-allocations": export.prodigy_allocations,
        "native-menus": menus.native_menus,
        "menu-records": menus.records,
        "normalise": normalise,
        "offered-budget": offered.calculate,
        "influence": influence_summary,
        "rank": rank_summary,
        "degree-menu": degree_comparison,
        "opportunity": opportunities,
    }
    return run_operation(
        config_path,
        {name: export.packaged(operation) for (name, operation) in (operations.items())},
    )
