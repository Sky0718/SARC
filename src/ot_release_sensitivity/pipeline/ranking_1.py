from __future__ import annotations
import math
from collections import Counter
from decimal import Decimal
from pathlib import Path
from typing import Any
from typing import Iterable
from typing import Iterator
from typing import Mapping
from typing import Sequence
from .ranking_types import DuckDbDataset, RankingSensitivityBuildError

def parse_boolean(value: Any) -> bool:
    if (isinstance(value, bool)):
        return value
    if (isinstance(value, int) and value in {0, 1}):
        return bool(value)
    text = str(value).strip().lower()
    if (text in {"true", "1", "yes"}):
        return True
    if (text in {"false", "0", "no", "", "none", "null"}):
        return False
    raise ValueError(f"Cannot parse Boolean value: {value}")

def finite_float(value: Any) -> float:
    result = float(value)
    if (not math.isfinite(result)):
        raise ValueError(f"Non-finite numeric value: {value}")
    return result

def finite_decimal(value: Any) -> Decimal:
    result = Decimal(str(value))
    if (not result.is_finite()):
        raise ValueError(f"Non-finite decimal value: {value}")
    return result

def linear_quantile(values: Iterable[float], probability: float) -> float | None:
    ordered = sorted((finite_float(value) for (value) in (values)))
    if (not ordered):
        return None
    if (probability < 0.0 or probability > 1.0):
        raise ValueError("Probability must lie in [0, 1]")
    if (len(ordered) == 1):
        return ordered[0]
    position = probability * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if (lower == upper):
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction

def weighted_quantile(
    values: Sequence[float], weights: Sequence[float], probability: float
) -> float | None:
    if (len(values) != len(weights)):
        raise ValueError("Values and weights must have equal length")
    if (probability < 0.0 or probability > 1.0):
        raise ValueError("Probability must lie in [0, 1]")
    pairs = sorted(
        (
            (finite_float(value), finite_float(weight))
            for ((value, weight)) in (zip(values, weights))
            if (finite_float(weight) > 0.0)
        )
    )
    if (not pairs):
        return None
    total = sum((weight for ((_, weight)) in (pairs)))
    threshold = probability * total
    cumulative = 0.0
    for (value, weight) in (pairs):
        cumulative += weight
        if (cumulative + 1e-15 >= threshold):
            return value
    return pairs[-1][0]

def weighted_fraction(values: Sequence[bool], weights: Sequence[float]) -> float | None:
    if (len(values) != len(weights)):
        raise ValueError("Values and weights must have equal length")
    total = sum((finite_float(weight) for (weight) in (weights)))
    if (total <= 0.0):
        return None
    return (
        sum(
            (
                finite_float(weight)
                for ((value, weight)) in (zip(values, weights))
                if (value)
            )
        )
        / total
    )

def average_descending_ranks(values: Sequence[float]) -> list[float]:
    checked = [finite_float(value) for (value) in (values)]
    order = sorted(range(len(checked)), key = lambda index: (-checked[index], index))
    ranks = [0.0] * len(checked)
    start = 0
    while (start < len(order)):
        end = start + 1
        score = checked[order[start]]
        while (end < len(order) and checked[order[end]] == score):
            end += 1
        average = (start + 1 + end) / 2.0
        for (position) in (range(start, end)):
            ranks[order[position]] = average
        start = end
    return ranks

def pearson_correlation(left: Sequence[float], right: Sequence[float]) -> float | None:
    if (len(left) != len(right)):
        raise ValueError("Vectors must have equal length")
    if (len(left) < 2):
        return None
    checked_left = [finite_float(value) for (value) in (left)]
    checked_right = [finite_float(value) for (value) in (right)]
    left_mean = sum(checked_left) / len(checked_left)
    right_mean = sum(checked_right) / len(checked_right)
    numerator = sum(
        (
            (left_value - left_mean) * (right_value - right_mean)
            for ((left_value, right_value)) in (zip(checked_left, checked_right))
        )
    )
    left_scale = sum(((value - left_mean) ** 2 for (value) in (checked_left)))
    right_scale = sum(((value - right_mean) ** 2 for (value) in (checked_right)))
    denominator = math.sqrt(left_scale * right_scale)
    if (denominator == 0.0):
        return None
    return numerator / denominator

def spearman_correlation(left: Sequence[float], right: Sequence[float]) -> float | None:
    return pearson_correlation(
        average_descending_ranks(left), average_descending_ranks(right)
    )

def kendall_tau_b(left: Sequence[float], right: Sequence[float]) -> float | None:
    if (len(left) != len(right)):
        raise ValueError("Vectors must have equal length")
    if (len(left) < 2):
        return None
    pairs = [(finite_float(x), finite_float(y)) for ((x, y)) in (zip(left, right))]
    pair_count = len(pairs)
    total_pairs = pair_count * (pair_count - 1) // 2
    left_counts = Counter((x for ((x, _)) in (pairs)))
    right_counts = Counter((y for ((_, y)) in (pairs)))
    left_ties = sum((count * (count - 1) // 2 for (count) in (left_counts.values())))
    right_ties = sum((count * (count - 1) // 2 for (count) in (right_counts.values())))
    denominator = math.sqrt((total_pairs - left_ties) * (total_pairs - right_ties))
    if (denominator == 0.0):
        return None
    y_values = sorted({y for ((_, y)) in (pairs)})
    y_index = {value: index + 1 for ((index, value)) in (enumerate(y_values))}
    tree = [0] * (len(y_values) + 1)

    def add(index: int) -> None:
        while (index < len(tree)):
            tree[index] += 1
            index += index & -index

    def prefix(index: int) -> int:
        total = 0
        while (index > 0):
            total += tree[index]
            index -= index & -index
        return total

    ordered = sorted(pairs)
    concordant = 0
    discordant = 0
    processed = 0
    start = 0
    while (start < pair_count):
        end = start + 1
        while (end < pair_count and ordered[end][0] == ordered[start][0]):
            end += 1
        for (_, y) in (ordered[start:end]):
            index = y_index[y]
            concordant += prefix(index - 1)
            discordant += processed - prefix(index)
        for (_, y) in (ordered[start:end]):
            add(y_index[y])
            processed += 1
        start = end
    return (concordant - discordant) / denominator

def tie_aware_top_set(
    scores: Mapping[str, float], cutoff: int, minimum_score: float = 0.0
) -> set[str] | None:
    if (cutoff < 1):
        raise ValueError("Cutoff must be positive")
    eligible = {
        str(identifier): finite_float(score)
        for ((identifier, score)) in (scores.items())
        if (finite_float(score) > minimum_score)
    }
    if (len(eligible) < cutoff):
        return None
    threshold = sorted(eligible.values(), reverse = True)[cutoff - 1]
    return {
        identifier for ((identifier, score)) in (eligible.items()) if (score >= threshold)
    }

def tie_aware_top_ranked_set(
    scores: Mapping[str, float], minimum_score: float = 0.0
) -> set[str] | None:
    eligible = {
        str(identifier): finite_float(score)
        for ((identifier, score)) in (scores.items())
        if (finite_float(score) > minimum_score)
    }
    if (not eligible):
        return None
    maximum = max(eligible.values())
    return {
        identifier for ((identifier, score)) in (eligible.items()) if (score == maximum)
    }

def jaccard_overlap(left: set[str] | None, right: set[str] | None) -> float | None:
    if (left is None or right is None):
        return None
    union = left | right
    if (not union):
        return None
    return len(left & right) / len(union)

def classify_turnover_reason(
    old_set: set[str] | None,
    new_set: set[str] | None,
    target_context: Mapping[str, Mapping[str, Any]],
) -> str:
    if (old_set is None or new_set is None):
        return "UNRESOLVED"
    changed = old_set ^ new_set
    if (not changed):
        return "UNCHANGED"
    reasons: set[str] = set()
    for (target) in (changed):
        context = target_context.get(target, {})
        mapping = parse_boolean(context.get("mapping_affected", False))
        record_category = str(context.get("record_category", "")).upper()
        transition_class = str(context.get("transition_class", "")).upper()
        if (
            mapping
            or record_category == "MAPPING_AFFECTED"
            or "MAPPING" in transition_class
        ):
            reasons.add("MAPPING_AFFECTED")
        elif (parse_boolean(
            context.get("old_present", target in old_set)
        ) and parse_boolean(context.get("new_present", target in new_set))):
            reasons.add("PERSISTENT_PAIR_RERANKING")
        elif (not parse_boolean(
            context.get("old_present", target in old_set)
        ) and parse_boolean(context.get("new_present", target in new_set))):
            reasons.add("ENTRANT")
        elif (parse_boolean(context.get("old_present", target in old_set)) and (
            not parse_boolean(context.get("new_present", target in new_set))
        )):
            reasons.add("EXIT")
        else:
            reasons.add("UNRESOLVED")
    if (len(reasons) == 1):
        return next(iter(reasons))
    return "MIXED"

def normalised_rank_displacements(
    old_scores: Sequence[float], new_scores: Sequence[float]
) -> tuple[list[float], list[float], list[float]]:
    if (len(old_scores) != len(new_scores)):
        raise ValueError("Vectors must have equal length")
    if (len(old_scores) < 2):
        raise ValueError("At least two targets are required")
    old_ranks = average_descending_ranks(old_scores)
    new_ranks = average_descending_ranks(new_scores)
    denominator = len(old_scores) - 1
    displacements = [
        abs(new - old) / denominator for ((old, new)) in (zip(old_ranks, new_ranks))
    ]
    return (old_ranks, new_ranks, displacements)

def evaluate_signal_assessment(
    changed_top_ranked_fraction: float | None,
    native_top_10_low_overlap_fraction: float | None,
    fixed_high_displacement_fraction: float | None,
    therapeutic_area_heterogeneity: bool,
    configuration: Mapping[str, Any],
) -> dict[str, bool]:
    thresholds = configuration["signal_assessment"]
    results = {
        "changed_tie_aware_top_ranked_set_fraction": changed_top_ranked_fraction
        is not None
        and changed_top_ranked_fraction
        >= finite_float(thresholds["changed_tie_aware_top_ranked_set_fraction"]),
        "native_top_10_low_overlap_fraction": native_top_10_low_overlap_fraction
        is not None
        and native_top_10_low_overlap_fraction
        >= finite_float(thresholds["native_top_10_low_overlap_fraction"]),
        "fixed_panel_high_displacement_fraction": fixed_high_displacement_fraction
        is not None
        and fixed_high_displacement_fraction
        >= finite_float(thresholds["fixed_panel_high_displacement_fraction"]),
        "therapeutic_area_heterogeneity": bool(therapeutic_area_heterogeneity),
    }
    results["at_least_one_prespecified_signal"] = any(results.values())
    return results

def evaluate_therapeutic_area_heterogeneity(
    rows: Sequence[Mapping[str, Any]], configuration: Mapping[str, Any]
) -> dict[str, Any]:
    area_config = configuration["therapeutic_areas"]
    minimum_panels = int(area_config["well_populated_minimum_panels"])
    fixed_well_populated = [
        row
        for (row) in (rows)
        if (int(row["eligible_disease_panel_count"]) >= minimum_panels)
    ]
    native_well_populated = [
        row
        for (row) in (rows)
        if (
            int(row["native_top_10_estimable_panel_count"]) >= minimum_panels
            and row.get("native_top_10_low_overlap_fraction") not in {None, ""}
        )
    ]
    fixed_values = [
        finite_float(row["fixed_top_ranked_set_changed_fraction"])
        for (row) in (fixed_well_populated)
    ]
    native_values = [
        finite_float(row["native_top_10_low_overlap_fraction"])
        for (row) in (native_well_populated)
    ]
    fixed_range = max(fixed_values) - min(fixed_values) if (fixed_values) else 0.0
    native_range = max(native_values) - min(native_values) if (native_values) else 0.0
    minimum_changed = int(
        area_config["heterogeneity_minimum_changed_panels_in_two_areas"]
    )
    fixed_qualifying = sum(
        (
            int(row["fixed_top_ranked_set_changed_count"]) >= minimum_changed
            for (row) in (fixed_well_populated)
        )
    )
    native_qualifying = sum(
        (
            int(row["native_top_10_low_overlap_count"])
            >= int(area_config["heterogeneity_minimum_changed_panels_in_two_areas"])
            for (row) in (native_well_populated)
        )
    )
    fixed_area_count_passed = len(fixed_well_populated) >= int(
        area_config["heterogeneity_minimum_well_populated_areas"]
    )
    native_area_count_passed = len(native_well_populated) >= int(
        area_config["heterogeneity_minimum_well_populated_areas"]
    )
    range_threshold = finite_float(area_config["heterogeneity_fraction_range"])
    fixed_passed = (
        fixed_area_count_passed
        and fixed_range >= range_threshold
        and (fixed_qualifying >= 2)
    )
    native_passed = (
        native_area_count_passed
        and native_range >= range_threshold
        and (native_qualifying >= 2)
    )
    passed = fixed_passed or native_passed
    if (fixed_passed and native_passed):
        passing_metric = "BOTH"
    elif (fixed_passed):
        passing_metric = "FIXED_TOP_RANKED_SET_CHANGED_FRACTION"
    elif (native_passed):
        passing_metric = "NATIVE_TOP_10_LOW_OVERLAP_FRACTION"
    else:
        passing_metric = "NONE"
    if (fixed_passed and (not native_passed)):
        well_populated_area_count = len(fixed_well_populated)
        qualifying_changed_area_count = fixed_qualifying
    elif (native_passed and (not fixed_passed)):
        well_populated_area_count = len(native_well_populated)
        qualifying_changed_area_count = native_qualifying
    else:
        well_populated_area_count = max(
            len(fixed_well_populated), len(native_well_populated)
        )
        qualifying_changed_area_count = max(fixed_qualifying, native_qualifying)
    return {
        "well_populated_area_count": well_populated_area_count,
        "fixed_well_populated_area_count": len(fixed_well_populated),
        "native_well_populated_area_count": len(native_well_populated),
        "fixed_changed_fraction_range": fixed_range,
        "native_low_overlap_fraction_range": native_range,
        "fixed_qualifying_changed_area_count": fixed_qualifying,
        "native_qualifying_changed_area_count": native_qualifying,
        "qualifying_changed_area_count": qualifying_changed_area_count,
        "fixed_metric_passed": fixed_passed,
        "native_metric_passed": native_passed,
        "passing_metric": passing_metric,
        "passed": passed,
    }

def configure_duckdb_resources(
    connection: Any, temporary_root: Path, configuration: Mapping[str, Any]
) -> Path:
    resources = configuration["resources"]
    temp_directory = (
        temporary_root / str(resources["duckdb_temp_subdirectory"])
    ).resolve()
    if (temp_directory.parent != temporary_root.resolve()):
        raise RankingSensitivityBuildError(
            "DuckDB temporary directory escaped the Ranking sensitivity build directory"
        )
    temp_directory.mkdir(parents = True, exist_ok = True)
    connection.execute(f"SET threads = {int(resources['threads'])}")
    connection.execute(f"SET memory_limit = {sql_literal(resources['memory_limit'])}")
    connection.execute(f"SET temp_directory = {sql_literal(temp_directory.as_posix())}")
    connection.execute(
        f"SET max_temp_directory_size = {sql_literal(resources['maximum_temp_directory_size'])}"
    )
    connection.execute("SET preserve_insertion_order = false")
    return temp_directory

def require_pyarrow() -> Any:
    try:
        import pyarrow
    except ImportError as exc:
        raise RankingSensitivityBuildError(
            "pyarrow is required to execute Ranking sensitivity"
        ) from exc
    return pyarrow

def sql_literal(value: Any) -> str:
    return "'" + str(value).replace("'", "''") + "'"

def sql_identifier(value: str) -> str:
    return '"' + str(value).replace('"', '""') + '"'

def relation_columns(connection: Any, relation: str) -> set[str]:
    return {
        row[1]
        for (row) in (
            connection.execute(f"PRAGMA table_info({sql_literal(relation)})").fetchall()
        )
    }

def create_input_views(
    connection: Any, paths: Mapping[str, Path], configuration: Mapping[str, Any]
) -> None:
    for (name, path) in (paths.items()):
        function = (
            "read_parquet" if (path.suffix.lower() == ".parquet") else "read_csv_auto"
        )
        options = (
            "" if (function == "read_parquet") else ", header=true, sample_size=-1"
        )
        connection.execute(
            f"CREATE VIEW {sql_identifier(name)} AS SELECT * FROM {function}({sql_literal(path.as_posix())}{options})"
        )
        observed = relation_columns(connection, name)
        expected = set(configuration["input_contract"][name])
        missing = sorted(expected - observed)
        if (missing):
            raise RankingSensitivityBuildError(
                f"Input {name} lacks required fields: {missing}"
            )

def query_dicts(connection: Any, query: str) -> list[dict[str, Any]]:
    cursor = connection.execute(query)
    fields = [item[0] for (item) in (cursor.description)]
    return [dict(zip(fields, row)) for (row) in (cursor.fetchall())]

def query_grouped_dicts(
    connection: Any, query: str, key_fields: Sequence[str], batch_size: int = 10000
) -> Iterator[tuple[tuple[Any, ...], list[dict[str, Any]]]]:
    cursor = connection.cursor()
    try:
        cursor.execute(query)
        fields = [item[0] for (item) in (cursor.description)]
        current_key: tuple[Any, ...] | None = None
        current_rows: list[dict[str, Any]] = []
        while (True):
            batch = cursor.fetchmany(batch_size)
            if (not batch):
                break
            for (raw_row) in (batch):
                row = dict(zip(fields, raw_row))
                key = tuple((row[field] for (field) in (key_fields)))
                if (current_key is not None and key != current_key):
                    yield (current_key, current_rows)
                    current_rows = []
                current_key = key
                current_rows.append(row)
        if (current_key is not None):
            yield (current_key, current_rows)
    finally:
        cursor.close()

def create_duckdb_dataset(
    connection: Any,
    table: str,
    fields: Sequence[str],
    types: Sequence[str],
    batch_size: int,
) -> DuckDbDataset:
    if (len(fields) != len(types)):
        raise RankingSensitivityBuildError(
            "DuckDB dataset field and type declarations differ in length"
        )
    if (len(fields) != len(set(fields))):
        raise RankingSensitivityBuildError(
            "DuckDB dataset field declarations must be unique"
        )
    if (batch_size < 1):
        raise RankingSensitivityBuildError("DuckDB dataset batch size must be positive")
    supported_types = {"VARCHAR", "BIGINT", "DOUBLE", "BOOLEAN"}
    if (any((str(data_type).upper() not in supported_types for (data_type) in (types)))):
        raise RankingSensitivityBuildError(
            "DuckDB dataset contains an unsupported Arrow insertion type"
        )
    connection.execute(f"DROP TABLE IF EXISTS {sql_identifier(table)}")
    declaration = ", ".join(
        (
            f"{sql_identifier(field)} {data_type}"
            for ((field, data_type)) in (zip(fields, types))
        )
    )
    connection.execute(f"CREATE TABLE {sql_identifier(table)} ({declaration})")
    return DuckDbDataset(table, fields, types, batch_size)

def arrow_type_for_duckdb(pyarrow: Any, data_type: str) -> Any:
    mapping = {
        "VARCHAR": pyarrow.string(),
        "BIGINT": pyarrow.int64(),
        "DOUBLE": pyarrow.float64(),
        "BOOLEAN": pyarrow.bool_(),
    }
    try:
        return mapping[str(data_type).upper()]
    except KeyError as exc:
        raise RankingSensitivityBuildError(
            f"Unsupported DuckDB type for Arrow insertion: {data_type}"
        ) from exc

def persist_duckdb_dataset_batch(
    connection: Any, dataset: DuckDbDataset, rows: Sequence[tuple[Any, ...]]
) -> None:
    if (not rows):
        return
    pyarrow = require_pyarrow()
    arrow_types = [
        arrow_type_for_duckdb(pyarrow, data_type) for (data_type) in (dataset.types)
    ]
    schema = pyarrow.schema(
        [
            pyarrow.field(field, arrow_type, nullable = True)
            for ((field, arrow_type)) in (zip(dataset.fields, arrow_types))
        ]
    )
    arrays = [
        pyarrow.array([row[index] for (row) in (rows)], type = arrow_type)
        for ((index, arrow_type)) in (enumerate(arrow_types))
    ]
    source = pyarrow.Table.from_arrays(arrays, schema = schema)
    source_name = f"ranking_sensitivity_arrow_insert_{dataset.table}"
    fields_sql = ", ".join((sql_identifier(field) for (field) in (dataset.fields)))
    values_sql = ", ".join(
        (
            f"CAST({sql_identifier(field)} AS {data_type})"
            for ((field, data_type)) in (zip(dataset.fields, dataset.types))
        )
    )
    connection.register(source_name, source)
    try:
        connection.execute(
            f"INSERT INTO {sql_identifier(dataset.table)} ({fields_sql}) SELECT {values_sql} FROM {sql_identifier(source_name)}"
        )
    finally:
        connection.unregister(source_name)
    dataset.persisted_count += len(rows)
    dataset.flush_count += 1
    dataset.flushed_batch_sizes.append(len(rows))

def flush_duckdb_dataset(connection: Any, dataset: DuckDbDataset) -> None:
    while (dataset.pending_rows):
        batch = dataset.pending_rows[: dataset.batch_size]
        persist_duckdb_dataset_batch(connection, dataset, batch)
        del dataset.pending_rows[: len(batch)]
    if (dataset.pending_count != 0 or dataset.persisted_count != dataset.row_count):
        raise RankingSensitivityBuildError(
            "DuckDB dataset flush did not reconcile logical and persisted row counts"
        )

def flush_duckdb_datasets(connection: Any, datasets: Iterable[Any]) -> None:
    for (dataset) in (datasets):
        if (isinstance(dataset, DuckDbDataset)):
            flush_duckdb_dataset(connection, dataset)
            if (dataset.pending_count != 0):
                raise RankingSensitivityBuildError(
                    "DuckDB dataset retained pending rows after mandatory flush"
                )

def append_duckdb_dataset_rows(
    connection: Any,
    dataset: DuckDbDataset,
    fields: Sequence[str],
    rows: Sequence[Mapping[str, Any]],
) -> None:
    if (tuple(fields) != dataset.fields):
        raise RankingSensitivityBuildError(
            "DuckDB dataset append fields do not match its explicit schema"
        )
    dataset.pending_rows.extend(
        (tuple((row.get(field) for (field) in (dataset.fields))) for (row) in (rows))
    )
    dataset.row_count += len(rows)
    while (dataset.pending_count >= dataset.batch_size):
        batch = dataset.pending_rows[: dataset.batch_size]
        persist_duckdb_dataset_batch(connection, dataset, batch)
        del dataset.pending_rows[: dataset.batch_size]

def joined_identifiers(values: set[str] | None) -> str:
    return "" if (values is None) else ";".join(sorted(values))

def membership_state(target: str, old_set: set[str], new_set: set[str]) -> str:
    if (target in old_set and target in new_set):
        return "PERSISTENT"
    if (target in old_set):
        return "OLD_ONLY"
    return "NEW_ONLY"
