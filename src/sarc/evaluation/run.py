from ..analysis.run import run_operation
from . import (
    continuous,
    eight_model,
    localisation,
    localisation_features,
    matched_return,
    prediction,
    representation,
    representation_design,
)

def continuous_eight(selections, endpoints, master_seeds):
    return eight_model.analyze(
        eight_model.normalize(selections), endpoints, master_seeds
    )

def matched_analysis(
    rows,
    split,
    rank_rows,
    eligibility,
    bootstrap,
    rank_binding,
    eligibility_binding,
    production = True,
):
    from ..transfer.policies import consensus

    checkpoint = lambda: None
    groups, labels, count = matched_return.load_ledger(
        rows, split, rank_binding, eligibility_binding, checkpoint, production
    )
    ranks = matched_return.rank_index(rank_rows, split, production)
    result = matched_return.analyze(
        groups,
        labels,
        ranks,
        eligibility,
        split,
        bootstrap,
        consensus,
        checkpoint,
        production,
    )
    return {"ledger_rows": count, "analysis": result}

def localised_features(baseline, changed, candidates, degs):
    def edges(rows):
        result = {
            (row["gene1"], row["gene2"]): row["combined_score"] for (row) in (rows)
        }
        if (len(result) != len(rows)):
            raise ValueError("Repeated graph edge")
        return result

    prepared = localisation_features.prepare_graph(edges(baseline), edges(changed))
    if (set(candidates) != set(degs)):
        raise ValueError("Mutation and expression sample populations differ")
    return {
        sample: localisation_features.strict_block(
            prepared, set(candidates[sample]), set(degs[sample])
        )
        for (sample) in (candidates)
    }

OPERATIONS = {
    "continuous_three": continuous.analyze,
    "continuous_eight": continuous_eight,
    "matched_return": matched_analysis,
    "representation_design": representation_design.make_panel,
    "representation": representation.analyse,
    "localisation_features": localised_features,
    "localisation_join": localisation.join_rows,
    "prediction_manifest": prediction.make_manifest,
    "selected_r2": prediction.fit_selected,
    "prediction_metrics": prediction.metric,
}

def run(config_path):
    return run_operation(config_path, OPERATIONS)
