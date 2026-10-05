import csv
import io
import itertools
import json
import math
import sys
from collections import defaultdict
from fractions import Fraction

MIN_NORMAL = Fraction(1, 2**1022)
MAX_FLOAT = Fraction.from_float(sys.float_info.max)
LIMIT = 2**53

def require(condition, message):
    if (not condition):
        raise ValueError(message)

def encode(value):
    if (isinstance(value, Fraction)):
        return {"numerator": value.numerator, "denominator": value.denominator}
    raise TypeError(type(value).__name__)

def validate_codec(codec, configurations, count = 85):
    mapping = codec["pipeline_to_canonical"]
    require(
        codec["status"] == "FROZEN_BIJECTIVE_IO_ADAPTER"
        and len(mapping) == len(set(mapping.values())) == count,
        "Complete exact sample codec required",
    )
    require(
        all(
            (
                isinstance(p, str)
                and p
                and ("_" not in p)
                and ("." not in p)
                and isinstance(c, str)
                and c
                for ((p, c)) in (mapping.items())
            )
        ),
        "Unsafe/ambiguous sample codec",
    )
    for (release, config) in (configurations.items()):
        order = config["pipeline_sample_order"]
        require(
            config["schema"] == "xc_personadrive_complete_cohort_input_v2"
            and config["context"] == "COLO_ORGANOID_SANGER"
            and (config["network_id"] == release),
            "Configuration identity differs",
        )
        require(
            config["pipeline_to_canonical"] == mapping
            and len(order) == len(set(order)) == count
            and (set(order) == set(mapping))
            and ([mapping[p] for (p) in (order)] == config["sample_order"]),
            "Configuration/codec roster differs",
        )
    return mapping

def pathways_from_stream(stream):
    pathways = defaultdict(set)
    text = stream.read().decode("utf-8")
    for (line) in (io.StringIO(text, newline = None)):
        fields = line.rstrip("\n").split("\t\t\t")
        require(len(fields) >= 2, "Original pathway triple-tab interface missing")
        for (gene) in (set(fields[1].split(",")[:-1])):
            pathways[gene].add(fields[0])
    return dict(pathways)

def parse_pps(payload, samples, pd):
    raw_rows = list(csv.reader(io.StringIO(payload.decode("utf-8-sig"))))
    require(
        len(raw_rows) == len(samples) + 1
        and all((len(row) == len(samples) + 1 for (row) in (raw_rows))),
        "PPS rectangular85x85 interface differs",
    )
    require(
        len(set(raw_rows[0][1:])) == len(samples)
        and set(raw_rows[0][1:]) == set(samples)
        and (len({r[0] for (r) in (raw_rows[1:])}) == len(samples))
        and ({r[0] for (r) in (raw_rows[1:])} == set(samples)),
        "Missing/duplicate PPS sample axis",
    )
    frame = pd.read_csv(io.BytesIO(payload), index_col = 0, float_precision = "legacy")
    require(
        frame.shape == (len(samples), len(samples))
        and set(frame.columns) == set(samples)
        and (set(frame.index) == set(samples))
        and all((str(t) == "float64" for (t) in (frame.dtypes))),
        "Source PPS float64/axis contract differs",
    )
    coefficients = {
        i: {j: float(frame[i][j]) for (j) in (samples)} for (i) in (samples)
    }
    require(
        all(
            (
                math.isfinite(v) and 0 <= v <= 1
                for (row) in (coefficients.values())
                for (v) in (row.values())
            )
        ),
        "PPS coefficient outside finite nonnegative source range",
    )
    require(
        all(
            (
                coefficients[i][j] == coefficients[j][i]
                for (i) in (samples)
                for (j) in (samples)
            )
        ),
        "Native PPS symmetry differs",
    )
    return coefficients

def parse_scores(payload, canonical_samples):
    reader = csv.DictReader(io.StringIO(payload.decode("utf-8-sig")))
    require(
        reader.fieldnames == ["sample_id", "gene", "score", "rank_position"],
        "Native score header differs",
    )
    rows = {m: [] for (m) in (canonical_samples)}
    for (item) in (reader):
        require(
            None not in item
            and item["sample_id"] in rows
            and isinstance(item["gene"], str)
            and item["gene"],
            "Invalid native score identity",
        )
        score = float(item["score"])
        rank = int(item["rank_position"])
        require(
            math.isfinite(score) and score >= 0 and (rank >= 1),
            "Invalid native score/rank scalar",
        )
        rows[item["sample_id"]].append((item["gene"], score, rank))
    for (sample, values) in (rows.items()):
        require(
            [row[2] for (row) in (values)] == list(range(1, len(values) + 1))
            and len({row[0] for (row) in (values)}) == len(values),
            "Native score rank/gene identities differ: " + sample,
        )
        require(
            values == sorted(values, key = lambda row: (-row[1], row[0])),
            "Native descending-score/lexical-tie order differs: " + sample,
        )
    return rows

def graph_statistics(stream, sample, mapping, pathways, nx):
    graph = nx.read_gml(stream)
    require(
        not graph.is_directed() and (not graph.is_multigraph()),
        "Expected simple undirected bipartite GML",
    )
    partitions = {}
    outliers = {}
    for (node, data) in (graph.nodes(data = True)):
        require(
            isinstance(node, str)
            and node
            and (type(data.get("bipartite")) is int)
            and (data["bipartite"] in (0, 1))
            and (graph.degree(node) > 0),
            "Invalid or isolated GML node",
        )
        partitions[node] = data["bipartite"]
        if (data["bipartite"] == 0):
            require("_" not in node, "Ambiguous source mutation-gene label")
        else:
            parts = node.split("_")
            require(
                len(parts) == 2 and parts[0] in mapping and parts[1],
                "Ambiguous/unknown source-labelled outlier",
            )
            outliers[node] = tuple(parts)
    require(
        all((partitions[a] != partitions[b] for ((a, b)) in (graph.edges()))),
        "Edge does not join opposite partitions",
    )
    statistics = {}
    for (gene) in (sorted((node for (node) in (graph) if (partitions[node] == 0)))):
        neighbours = tuple(sorted((outliers[n] for (n) in (graph.neighbors(gene)))))
        require(
            all((d != gene for ((j, d)) in (neighbours))),
            "Forbidden source self-gene edge",
        )
        counts, minimum_positive, terms = (defaultdict(int), {}, 0)
        contributions = []
        for (donor, other_gene) in (neighbours):
            count = None
            if (gene in pathways and other_gene in pathways):
                count = len(pathways[gene] & pathways[other_gene])
                require(
                    math.isfinite(float(count)) and int(float(count)) == count,
                    "Shared-pathway count not exactly binary64 representable",
                )
                counts[donor] += count
                terms += 1
                if (count):
                    minimum_positive[donor] = min(
                        count, minimum_positive.get(donor, count)
                    )
            contributions.append((donor, other_gene, count))
        statistics[gene] = {
            "neighbours": neighbours,
            "contributions": contributions,
            "A": dict(sorted(counts.items())),
            "minimum_positive_count": minimum_positive,
            "terms": terms,
            "L": sum(counts.values()),
        }
    own = {gene for ((donor, gene)) in (outliers.values()) if (donor == sample)}
    return (
        statistics,
        {
            "nodes": len(graph),
            "edges": graph.number_of_edges(),
            "mutation_nodes": len(statistics),
            "own_labelled_outlier_count": len(own),
        },
    )

def recipient_invariance(reference, statistics, sample):
    for (gene, record) in (statistics.items()):
        if (gene in reference):
            require(
                reference[gene]["neighbours"] == record["neighbours"],
                "Recipient-labelled-neighbour invariance failed: "
                + sample
                + "/"
                + gene,
            )
            require(
                all(
                    (
                        reference[gene][key] == record[key]
                        for (key) in (("A", "L", "minimum_positive_count", "terms"))
                    )
                ),
                "Derived pathway/kernel semantics differ for same neighbours",
            )
            reference[gene]["recipients"].append(sample)
        else:
            reference[gene] = dict(record, recipients = [sample])

def dyadic_row(coefficients):
    require(
        all(
            (
                type(value) is float and math.isfinite(value) and (0 <= value <= 1)
                for (value) in (coefficients.values())
            )
        ),
        "Invalid parsed binary64 coefficients",
    )
    ratios = {j: value.as_integer_ratio() for ((j, value)) in (coefficients.items())}
    denominator = max((d for ((n, d)) in (ratios.values())))
    numerators = {j: n * (denominator // d) for ((j, (n, d))) in (ratios.items())}
    require(
        all((denominator % d == 0 for ((n, d)) in (ratios.values()))),
        "Binary64 common denominator not dyadic",
    )
    return {
        "denominator": denominator,
        "numerators": numerators,
        "total_numerator": sum(numerators.values()),
    }

def reconstruct(record, row, stored):
    require(
        type(stored) is float and math.isfinite(stored) and (stored >= 0),
        "Invalid stored native binary64 score",
    )
    denominator, numerators = (row["denominator"], row["numerators"])
    require(
        type(record["terms"]) is int
        and record["terms"] >= 0
        and all((type(a) is int and a >= 0 for (a) in (record["A"].values())))
        and (sum(record["A"].values()) == record["L"]),
        "Invalid integer source sufficient statistics",
    )
    require(set(record["A"]) <= set(numerators), "Unknown donor in A")
    h = sum((numerators[j] * value for ((j, value)) in (record["A"].items())))
    exact = Fraction(h, denominator)
    require(
        record["terms"] > 0 or h == 0,
        "Zero native pathway terms cannot have positive reconstructed score",
    )
    for (donor, count) in (record["minimum_positive_count"].items()):
        product = Fraction(numerators[donor] * count, denominator)
        require(
            product == 0 or product >= MIN_NORMAL, "Nonzero native product is subnormal"
        )
    operations = 2 * record["terms"] + 2
    require(operations < LIMIT, "Gamma operation count invalid")
    gamma = Fraction(operations, LIMIT - operations)
    require(
        (1 + gamma) * exact <= MAX_FLOAT, "Conservative no-overflow precondition failed"
    )
    saved = Fraction.from_float(stored)
    error, bound = (abs(saved - exact), gamma * exact)
    require(
        saved == 0 if (exact == 0 or record["terms"] == 0) else error <= bound,
        "All-emitted source score reconstruction envelope failed",
    )
    total = row["total_numerator"]
    require(
        total > 0 or (exact == 0 and saved == 0),
        "Zero-kernel native score must be exactly zero",
    )
    return {
        "H": h,
        "R": exact,
        "stored_binary64_exact": saved,
        "absolute_error": error,
        "gamma": gamma,
        "error_bound": bound,
        "pathway_eligible_terms": record["terms"],
        "p": operations,
        "T": Fraction(h, total) if (total) else None,
        "reconstruction_passed": True,
    }

def allocate(scores, capacity = 10):
    require(
        type(capacity) is int
        and capacity >= 0
        and all((isinstance(g, str) and g for (g) in (scores))),
        "Invalid exact allocation input",
    )
    require(
        all(
            (
                type(s) in (int, Fraction)
                and s >= 0
                or (type(s) is float and math.isfinite(s) and (s >= 0))
                for (s) in (scores.values())
            )
        ),
        "Invalid allocation scalar",
    )
    ordered = sorted(scores, key = lambda g: (-scores[g], g))
    weights, used = ({}, 0)
    for (value, block) in (itertools.groupby(ordered, key = scores.__getitem__)):
        genes = list(block)
        slots = min(len(genes), max(0, capacity - used))
        weights.update({g: Fraction(slots, len(genes)) for (g) in (genes)})
        used += len(genes)
    require(
        set(weights) == set(scores)
        and sum(weights.values(), Fraction()) == min(capacity, len(scores)),
        "Allocation lost exact menu or capacity",
    )
    return weights

def displacement(first, second):
    require(
        set(first) == set(second)
        and sum(first.values(), Fraction()) == sum(second.values(), Fraction()),
        "Allocation menus/capacities differ",
    )
    return sum((abs(first[g] - second[g]) for (g) in (first)), Fraction()) / 2

def distribution(values):
    require(bool(values), "Undefined empty summary cannot be silently averaged")
    values = sorted(values)
    n = len(values)
    return {
        "count": n,
        "minimum": values[0],
        "median": (values[(n - 1) // 2] + values[n // 2]) / 2,
        "maximum": values[-1],
        "uniform_model_mean": sum(values, Fraction()) / n,
    }

def model_record(sample, release, scores, statistics, row, eligibility, gene_stream):
    require(
        set(statistics) == {g for ((g, s, rank)) in (scores)},
        "Full emitted graph/score gene set differs: " + sample,
    )
    menu = set(statistics) & set(eligibility)
    native, prior, mathematical, gene_records = ({}, {}, {}, [])
    for (gene, stored, rank) in (scores):
        record = statistics[gene]
        try:
            algebra = reconstruct(record, row, stored)
        except Exception as error:
            raise ValueError(
                release + "/" + sample + "/" + gene + ": " + str(error)
            ) from error
        eligible = gene in menu
        if (eligible):
            native[gene], prior[gene], mathematical[gene] = (
                stored,
                record["L"],
                algebra["H"],
            )
        gene_records.append(
            {
                "release": release,
                "sample_id": sample,
                "gene": gene,
                "native_rank_position": rank,
                "native_score": stored,
                "eligible": eligible,
                "A_by_source_model": record["A"],
                "L": record["L"],
                "F": Fraction(record["L"], 85) if (eligible) else None,
                **algebra,
            }
        )
    w_native, w_f = (allocate(native), allocate(prior))
    w_t = allocate(mathematical) if (row["total_numerator"]) else None
    for (record) in (gene_records):
        g = record["gene"]
        record.update(
            {
                "W_native": w_native.get(g),
                "W_F": w_f.get(g),
                "W_T": None if (w_t is None) else w_t.get(g),
            }
        )
        gene_stream.write(json.dumps(record, default = encode, allow_nan = False) + "\n")
    return {
        "release": release,
        "sample_id": sample,
        "technical_status": "COMPLETE",
        "emitted_gene_count": len(scores),
        "emitted_reconstruction_pass_count": len(scores),
        "menu": sorted(menu),
        "offered_count": len(menu),
        "returned_mass": min(10, len(menu)),
        "empty_menu": len(menu) == 0,
        "all_returned_menu": len(menu) <= 10,
        "informative_menu": len(menu) > 10,
        "zero_kernel": row["total_numerator"] == 0,
        "kernel_sum": Fraction(row["total_numerator"], row["denominator"]),
        "W_native": w_native,
        "W_F": w_f,
        "W_T": w_t,
        "native_F_equal": w_native == w_f,
        "D_native_F": displacement(w_native, w_f),
        "native_T_equal": None if (w_t is None) else w_native == w_t,
        "D_native_T": None if (w_t is None) else displacement(w_native, w_t),
    }

def release_summary(release, records, expected_count = 85):
    require(
        len(records) == len({r["sample_id"] for (r) in (records)}) == expected_count
        and all((r["technical_status"] == "COMPLETE" for (r) in (records))),
        "No complete full-population mean after technical loss",
    )
    defined = [r for (r) in (records) if (r["D_native_T"] is not None)]
    return {
        "release": release,
        "model_count": expected_count,
        "technical_unavailable_count": 0,
        "native_F_equal_count": sum((r["native_F_equal"] for (r) in (records))),
        "native_F_equal_denominator": expected_count,
        "native_F_displacement": distribution([r["D_native_F"] for (r) in (records)]),
        "native_T_equal_count": sum((r["native_T_equal"] for (r) in (defined))),
        "native_T_defined_denominator": len(defined),
        "native_T_undefined_count": expected_count - len(defined),
        "native_T_displacement": distribution([r["D_native_T"] for (r) in (defined)])
        if (defined)
        else None,
        "empty_menu_count": sum((r["empty_menu"] for (r) in (records))),
        "all_returned_menu_count": sum((r["all_returned_menu"] for (r) in (records))),
        "informative_menu_count": sum((r["informative_menu"] for (r) in (records))),
        "zero_kernel_count": sum((r["zero_kernel"] for (r) in (records))),
        "all_emitted_reconstruction_pass_count": sum(
            (r["emitted_reconstruction_pass_count"] for (r) in (records))
        ),
    }
