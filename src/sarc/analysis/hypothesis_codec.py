import math
import struct
import sys
from fractions import Fraction

U = Fraction(1, 2**53)
MIN_NORMAL = Fraction.from_float(sys.float_info.min)
MAX_FINITE = Fraction.from_float(sys.float_info.max)

def require(condition, message):
    if (not condition):
        raise ValueError(message)

def text_vector(values):
    return (
        None
        if (values is None)
        else [
            None if (x is None) else bytes.fromhex(x).decode("utf-8")
            for (x) in (values)
        ]
    )

class Numeric:
    def __init__(self, value):
        self.kind = value["t"]
        self.names = text_vector(value["names"])
        self.dim = value["dim"]
        self.dimnames = (
            None
            if (value["dimnames"] is None)
            else [text_vector(x) for (x) in (value["dimnames"])]
        )
        self.bits = bytes.fromhex(value["bits"])
        size = 8 if (self.kind == "double") else 4
        require(len(self.bits) % size == 0, "Numeric byte length invalid")
        self.values = [
            x[0]
            for (x) in (struct.iter_unpack("<d" if (size == 8) else "<i", self.bits))
        ]
        self.states = value.get(
            "states",
            "".join(("A" if (x == -2147483648) else "F" for (x) in (self.values))),
        )
        require(len(self.states) == len(self.values), "Numeric state count invalid")
        if (self.kind == "double"):
            require(
                all(
                    (
                        s == "F"
                        and math.isfinite(x)
                        or (s in "AN" and math.isnan(x))
                        or (s == "P" and x == math.inf)
                        or (s == "M" and x == -math.inf)
                        for ((x, s)) in (zip(self.values, self.states))
                    )
                ),
                "Numeric state/IEEE mismatch",
            )
        require(
            self.names is None or len(self.names) == len(self.values),
            "Numeric names mismatch",
        )
        require(
            self.dim is None or math.prod(self.dim) == len(self.values),
            "Numeric dimensions mismatch",
        )

    def token(self, index):
        size = 8 if (self.kind == "double") else 4
        return {
            "state": self.states[index],
            "bits": self.bits[index * size : (index + 1) * size].hex(),
        }

def decode(value):
    kind = value["t"]
    if (kind == "null"):
        return None
    if (kind in ("double", "integer", "logical")):
        return Numeric(value)
    if (kind == "character"):
        return text_vector(value["values"])
    require(kind == "list", "Unexpected serialized type")
    values = [decode(x) for (x) in (value["values"])]
    names = text_vector(value["names"])
    if (names is None):
        return values
    require(
        len(names) == len(values)
        and len(set(names)) == len(names)
        and (None not in names),
        "Named list identity invalid",
    )
    return dict(zip(names, values))

def scalar(value):
    items = value.values if (isinstance(value, Numeric)) else value
    require(isinstance(items, list) and len(items) == 1, "Scalar field expected")
    return items[0]

def listify(value):
    return [] if (value is None) else value if (isinstance(value, list)) else [value]

def unique_axis(values, title):
    require(
        isinstance(values, list)
        and all((isinstance(x, str) and x for (x) in (values)))
        and (len(set(values)) == len(values)),
        title + " axis invalid",
    )
    return values

def sum_envelope(values):
    require(
        all(
            (
                type(x) in (float, int) and math.isfinite(x) and (x >= 0)
                for (x) in (values)
            )
        ),
        "Nonnegative finite summands required",
    )
    exact = [Fraction(x) for (x) in (values)]
    require(
        all((x == 0 or x >= MIN_NORMAL for (x) in (exact))),
        "Positive subnormal outside frozen sum contract",
    )
    total = sum(exact, Fraction())
    factor = (len(values) + 1) * U
    require(factor < 1, "Sum gamma undefined")
    envelope = factor / (1 - factor) * total
    require(
        total + envelope <= MAX_FINITE, "Possible aggregate overflow outside contract"
    )
    return (total, envelope)

def mask_operator(rows, genes, columns, vector = False):
    require(
        len(columns) == len(genes) and all((len(c) == len(rows) for (c) in (columns))),
        "Rectangle dimensions differ",
    )
    require(
        all((all((math.isfinite(v) and v >= 0 for (v) in (c))) for (c) in (columns))),
        "Mask input requires finite clipped values",
    )
    if (not rows or not genes):
        return (set(), "EMPTY_ANALYTICAL_RECTANGLE", {})
    if (vector):
        require(len(rows) == 1, "Native vector requires one pathway")
        return (
            set(rows),
            "NATIVE_VECTOR_OR_SCALAR",
            {g: sum_envelope(c) for ((g, c)) in (zip(genes, columns))},
        )
    if (len(genes) < 4):
        selected = list(range(len(rows)))
        branch = "MATRIX_FEWER_THAN_FOUR"
    else:
        selected = [
            i
            for (i) in (range(len(rows)))
            if (2 * sum((c[i] > 0 for (c) in (columns))) < len(genes))
        ]
        branch = (
            "FILTER_ZERO_ROWS"
            if (not selected)
            else "FILTER_ONE_ROW"
            if (len(selected) == 1)
            else "FILTER_MULTIPLE_ROWS"
        )
        if (not any((c[i] > 0 for (c) in (columns) for (i) in (selected)))):
            require(
                len(rows) > 1,
                "UNAVAILABLE_SOURCE_MATRIX_DIMENSION_DROP_IN_ALL_ZERO_FALLBACK",
            )
            selected = list(range(len(rows)))
            branch += "_ALL_ZERO_FALLBACK"
        elif (len(selected) >= 2):
            selected = [
                i for (i) in (selected) if (any((c[i] > 0 for (c) in (columns))))
            ]
            branch += "_ZERO_ONLY_ROWS_PRUNED"
    totals = {
        g: sum_envelope([columns[j][i] for (i) in (selected)])
        for ((j, g)) in (enumerate(genes))
    }
    return ({rows[i] for (i) in (selected)}, branch, totals)

def pathway_dictionary(tree):
    require(isinstance(tree, dict), "Named fixed pathway dictionary required")
    return {
        p: set(row["members"])
        if (scalar(row["status"]) == "EXACT_SOURCE_MEMBERSHIP")
        else None
        for ((p, row)) in (tree.items())
    }

def checkpoint(tree, identity):
    record = tree["identity"]
    for (key, expected) in (identity.items()):
        require(scalar(record[key]) == expected, "Checkpoint identity differs: " + key)
    genes = unique_axis(listify(tree["axes"]["mutation_genes"]), "Native mutation")
    rows = unique_axis(listify(tree["axes"]["evaluated_pathways"]), "Evaluated pathway")
    status = scalar(tree["restored"]["status"])
    array = tree["restored"]["values"]
    if (status == "SUCCESSFUL_EMPTY"):
        require(
            array is None and tree["raw_influence"] is None,
            "Successful empty checkpoint not null",
        )
        return {
            "status": "SUCCESSFUL_EMPTY",
            "genes": genes,
            "rows": rows,
            "columns": {},
            "mask": set(),
            "branch": "SUCCESSFUL_EMPTY",
            "totals": {},
            "differences": tree["differences"],
            "axis_status": status,
        }
    require(
        isinstance(array, Numeric) and array.kind == "double",
        "Double native influence required",
    )
    require(
        status
        in (
            "AXES_VERIFIED",
            "AXES_RESTORED_EXACT_NUMERICAL_VALUES",
            "VALID_ORIGINAL_VECTOR_BRANCH_PRESERVED",
            "SCALAR_GENE_IDENTITY_VERIFIED",
        ),
        "Unsupported authoritative axis status",
    )
    if (array.dim is None):
        require(
            len(rows) == 1
            and array.names == genes
            and (len(array.values) == len(genes)),
            "Native vector identity mismatch",
        )
    else:
        require(
            array.dim == [len(rows), len(genes)] and array.dimnames == [rows, genes],
            "Native matrix authoritative axes differ",
        )
    raw = tree["raw_influence"]
    returned = tree["axes"]["returned"]
    require(
        isinstance(raw, Numeric)
        and isinstance(returned, Numeric)
        and (raw.bits == returned.bits == array.bits)
        and (raw.states == returned.states == array.states),
        "Raw/returned/restored values differ",
    )
    require(
        raw.names == returned.names
        and raw.dim == returned.dim
        and (raw.dimnames == returned.dimnames),
        "Raw and returned authoritative structure differs",
    )
    if (status == "AXES_VERIFIED"):
        require(
            raw.dim == array.dim and raw.dimnames == array.dimnames,
            "Verified matrix must not change axes",
        )
    elif (status == "AXES_RESTORED_EXACT_NUMERICAL_VALUES"):
        require(
            raw.dim is None
            and len(genes) == 1
            and (len(rows) > 1)
            and (raw.names == rows),
            "Restored singleton column lacks authoritative raw names",
        )
    elif (status == "VALID_ORIGINAL_VECTOR_BRANCH_PRESERVED"):
        require(
            raw.dim is None
            and len(rows) == 1
            and (len(genes) > 1)
            and (raw.names == genes),
            "Preserved vector raw names differ",
        )
    else:
        require(
            raw.dim is None
            and len(rows) == len(genes) == 1
            and (raw.names in (None, rows, genes)),
            "Scalar raw identity differs",
        )
    columns = {}
    clipped = []
    for (j, gene) in (enumerate(genes)):
        start = j * len(rows)
        values = array.values[start : start + len(rows)]
        states = array.states[start : start + len(rows)]
        tokens = [array.token(i) for (i) in (range(start, start + len(rows)))]
        x = [
            0.0 if (state in "ANM" or value < 0) else value
            for ((value, state)) in (zip(values, states))
        ]
        columns[gene] = {"values": values, "states": states, "tokens": tokens, "X": x}
        clipped.append(x)
    mask, branch, totals = mask_operator(rows, genes, clipped, array.dim is None)
    return {
        "status": "AVAILABLE",
        "genes": genes,
        "rows": rows,
        "columns": columns,
        "mask": mask,
        "branch": branch,
        "totals": totals,
        "differences": tree["differences"],
        "axis_status": status,
    }

def aggregate_check(cp, tree):
    actual = tree["aggregate"]
    expected = {g: item for ((g, item)) in (cp["totals"].items()) if (item[0] > 0)}
    if (actual is None):
        require(
            cp["status"] == "SUCCESSFUL_EMPTY" and (not expected),
            "Unexpected null aggregate",
        )
        return []
    require(
        isinstance(actual, Numeric)
        and actual.kind == "double"
        and all((s == "F" for (s) in (actual.states))),
        "Aggregate must be finite double",
    )
    names = actual.names or []
    require(
        set(names) == set(expected)
        and len(names) == len(expected) == len(actual.values),
        "Complete positive aggregate gene set differs",
    )
    require(
        all((x > 0 for (x) in (actual.values)))
        and all((a >= b for ((a, b)) in (zip(actual.values, actual.values[1:])))),
        "Saved aggregate ordering/positivity invalid",
    )
    checks = []
    for (gene, value, i) in (zip(names, actual.values, range(len(names)))):
        exact, envelope = expected[gene]
        residual = Fraction(value) - exact
        require(
            abs(residual) <= envelope,
            "Native aggregate outside prospective sum envelope: " + gene,
        )
        checks.append(
            {
                "gene": gene,
                "saved": actual.token(i),
                "exact_sum": exact,
                "residual": residual,
                "envelope": envelope,
            }
        )
    return checks

def prizes(cp, pathway, members):
    if (pathway in cp.setdefault("_prize_cache", {})):
        return cp["_prize_cache"][pathway]
    differences = cp["differences"]
    if (
        members is None
        or not isinstance(differences, Numeric)
        or differences.kind != "double"
        or (differences.names is None)
    ):
        return {"status": "UNAVAILABLE_PRIZE_INTERFACE"}
    if (
        len(set(differences.names)) != len(differences.names)
        or None in differences.names
    ):
        return {"status": "UNAVAILABLE_PRIZE_NAMES"}
    indices = [i for ((i, gene)) in (enumerate(differences.names)) if (gene in members)]
    if (any(
        (differences.states[i] != "F" or differences.values[i] < 0 for (i) in (indices))
    )):
        return {"status": "UNAVAILABLE_NONFINITE_OR_NEGATIVE_PRIZES"}
    try:
        total, envelope = sum_envelope([differences.values[i] for (i) in (indices)])
    except ValueError as error:
        return {"status": "UNAVAILABLE_PRIZE_SUM_CONTRACT", "reason": str(error)}
    value = {
        "status": "RECOVERED_STORED_VECTOR_NO_ENRICHMENT",
        "pathway": pathway,
        "ordered": [[differences.names[i], differences.token(i)] for (i) in (indices)],
        "exact_sum": total,
        "native_sum_envelope": envelope,
    }
    cp["_prize_cache"][pathway] = value
    return value
