import csv
from fractions import Fraction
from pathlib import Path
import numpy as np
from .contracts import SEEDS, STATES, require, save, validate_models

def allocate(scores, positions):
    if (len(positions) <= 10):
        return {int(position): Fraction(1) for (position) in (positions)}
    order = sorted(
        (int(position) for (position) in (positions)), key = lambda position: -scores[position]
    )
    cutoff = scores[order[9]]
    above = [position for (position) in (order) if (scores[position] > cutoff)]
    ties = [position for (position) in (order) if (scores[position] == cutoff)]
    weights = {position: Fraction(1) for (position) in (above)}
    weights.update({position: Fraction(10 - len(above), len(ties)) for (position) in (ties)})
    return weights

def matrix(path, dimensions):
    values = np.fromfile(path, dtype = "<f8")
    require(
        values.size == int(np.prod(dimensions)) and np.isfinite(values).all(),
        "Finite complete binary array required",
    )
    return values.reshape(dimensions, order = "F")

def validate_certificates(path, samples):
    with Path(path).open(encoding = "utf-8-sig", newline = "") as stream:
        rows = list(csv.DictReader(stream))
    selected = [row for (row) in (rows) if (int(row["budget"]) == 10)]
    expected = {
        (state, sample) for (state) in (("G0_H0", "G1_H0", "G0_H1", "G1_H1")) for (sample) in (samples)
    }
    keys = [(row["state"], row["sample_id"]) for (row) in (selected)]
    require(
        len(keys) == len(set(keys)) and set(keys) == expected,
        "Complete 236 cutoff certificates required",
    )
    require(
        all(str(row["certified"]).lower() in ("true", "1") for (row) in (selected)),
        "Uncertified allocation cannot enter evaluation",
    )

def build_allocations(computation, models, output):
    computation = Path(computation)
    validate_models(models)
    genes = (computation / "genes.txt").read_text(encoding = "utf-8").splitlines()
    samples = [model["sample_id"] for (model) in (models)]
    require(len(genes) == len(set(genes)) == 7388, "Complete common input domain required")
    require(
        (computation / "models.txt").read_text(encoding = "utf-8").splitlines() == samples,
        "Numerical model order differs",
    )
    validate_certificates(computation / "cutoff_certificates.csv", samples)
    dimensions = (len(genes), len(models))
    mutation = matrix(computation / "mutation.bin", dimensions)
    require(np.isin(mutation, (0, 1)).all(), "Binary mutation matrix required")
    arrays = {
        graph + evidence: matrix(
            computation / (graph + "_" + evidence + "_scores.bin"), dimensions
        )
        for (graph) in (("G0", "G1"))
        for (evidence) in (("H0", "H1"))
    }
    degree = {
        "degree_" + graph: matrix(computation / (graph + "_degree.bin"), (len(genes),))
        for (graph) in (("G0", "G1"))
    }
    require(
        all(
            np.greater_equal(value, 0).all() and np.equal(value, np.floor(value)).all()
            for (value) in (degree.values())
        ),
        "Nonnegative graph degrees required",
    )
    records = []
    with (output / "full_query_ranks.csv").open("x", encoding = "utf-8", newline = "") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames = [
                "sample_id",
                "donor_id",
                "state",
                "gene",
                "rank_position",
                "score",
                "allocation",
            ],
            lineterminator = "\n",
        )
        writer.writeheader()
        for (column, model) in (enumerate(models)):
            positions = np.flatnonzero(mutation[:, column] == 1).tolist()
            states = {}
            for (state) in (STATES):
                values = arrays[state][:, column] if (state in arrays) else degree[state]
                weights = allocate(values, positions)
                require(
                    sum(weights.values(), Fraction()) == min(10, len(positions)),
                    "Incomplete ten-slot allocation",
                )
                states[state] = {
                    "status": "AVAILABLE",
                    "weights": {
                        genes[position]: str(weight)
                        for (position, weight) in (sorted(weights.items()))
                    },
                }
                order = sorted(positions, key = lambda position: -values[position])
                for (rank, position) in (enumerate(order, 1)):
                    writer.writerow(
                        {
                            "sample_id": model["sample_id"],
                            "donor_id": model["donor_id"],
                            "state": state,
                            "gene": genes[position],
                            "rank_position": rank,
                            "score": format(values[position], ".17g"),
                            "allocation": str(weights.get(position, Fraction())),
                        }
                    )
            records.append(
                {
                    **model,
                    "query": [genes[position] for (position) in (positions)],
                    "states": states,
                }
            )
    result = {"budget": 10, "seeds": SEEDS, "models": records}
    save(output / "allocations.json", result)
    return result
