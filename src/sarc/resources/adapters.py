import re
from dataclasses import dataclass
from typing import Iterable, Mapping
from .contracts import (
    ContractError,
    ScoreObservation,
    observe_score,
    require_count,
    require_text,
)

@dataclass(frozen = True)
class FieldMap:
    fields: tuple[tuple[str, str], ...]

    def __post_init__(self):
        if (not isinstance(self.fields, tuple)):
            raise ContractError("field declarations must be an explicit tuple")
        for (canonical, source) in (self.fields):
            require_text(canonical, "canonical field")
            require_text(source, "source field")
        if (len({key for ((key, _)) in (self.fields)}) != len(self.fields)):
            raise ContractError("canonical field is declared more than once")
        if (len({value for ((_, value)) in (self.fields)}) != len(self.fields)):
            raise ContractError(
                "one source field cannot silently serve multiple meanings"
            )

    def require(self, *names: str):
        declared = {key for ((key, _)) in (self.fields)}
        if (not set(names) <= declared):
            raise ContractError("required canonical fields are not declared")

    def read(self, row: Mapping[str, object], name: str) -> object:
        aliases = dict(self.fields)
        if (name not in aliases):
            raise ContractError(f"undeclared field: {name}")
        source = aliases[name]
        if (source not in row):
            raise ContractError(f"missing declared source field: {source}")
        return row[source]

@dataclass(frozen = True)
class OTAssociationSchema:
    family: str
    fields: FieldMap
    expected_aggregation_type: str | None = None
    allowed_aggregation_types: frozenset[str] = frozenset()
    allowed_selector_values: frozenset[str | None] = frozenset()

    def __post_init__(self):
        if (self.family not in {
            "association_overall_direct",
            "association_by_datatype_direct",
            "association_by_datasource_direct",
        }):
            raise ContractError("unknown Open Targets association family")
        self.fields.require("disease_id", "target_id", "score", "evidence_count")
        if (not isinstance(self.allowed_aggregation_types, frozenset) or not isinstance(
            self.allowed_selector_values, frozenset
        )):
            raise ContractError("selector domains must be explicit frozen sets")
        for (value) in (self.allowed_aggregation_types):
            require_text(value, "aggregation type")
        for (value) in (self.allowed_selector_values):
            if (value is not None):
                require_text(value, "selector value")
        if (self.expected_aggregation_type is not None):
            self.fields.require("aggregation_type", "selector")
            if (self.expected_aggregation_type not in self.allowed_aggregation_types):
                raise ContractError(
                    "expected aggregation type is not in the declared domain"
                )
        elif (self.allowed_aggregation_types):
            raise ContractError(
                "aggregation domain requires an expected aggregation type"
            )
        if (self.allowed_selector_values):
            self.fields.require("selector")
        if (self.expected_aggregation_type is not None and (
            not self.allowed_selector_values
        )):
            raise ContractError(
                "an aggregated source requires a closed selector domain"
            )
        if (self.family != "association_overall_direct" and (
            not self.allowed_selector_values
        )):
            raise ContractError(
                "datatype and datasource families require a declared selector"
            )

@dataclass(frozen = True)
class AssociationRecord:
    resource: str
    release: str
    family: str
    disease_id: str
    target_id: str
    selector: str | None
    observation: ScoreObservation
    evidence_count: int | None
    eligible: bool
    ineligible_reason: str | None

def adapt_ot_association(
    row: Mapping[str, object],
    release: str,
    schema: OTAssociationSchema,
    *,
    skip_other_aggregations: bool = False,
) -> AssociationRecord | None:
    require_text(release, "release")
    if (not isinstance(skip_other_aggregations, bool)):
        raise ContractError("skip_other_aggregations must be bool")
    fields = schema.fields
    if (schema.expected_aggregation_type is not None):
        aggregation_type = require_text(
            fields.read(row, "aggregation_type"), "aggregation type"
        )
        if (aggregation_type not in schema.allowed_aggregation_types):
            raise ContractError("unknown aggregation type")
        if (aggregation_type != schema.expected_aggregation_type):
            if (skip_other_aggregations):
                return None
            raise ContractError("row belongs to another declared aggregation type")
    selector = (
        fields.read(row, "selector") if (schema.allowed_selector_values) else None
    )
    if (selector is not None):
        require_text(selector, "selector")
    if (
        schema.allowed_selector_values
        and selector not in schema.allowed_selector_values
    ):
        raise ContractError("selector is outside its declared domain")
    disease = require_text(fields.read(row, "disease_id"), "disease identifier")
    target = require_text(fields.read(row, "target_id"), "target identifier")
    if (re.fullmatch("ENSG[0-9]{11}", target) is None):
        raise ContractError("target identifier is not a full Ensembl gene identifier")
    raw_count = fields.read(row, "evidence_count")
    evidence_count = (
        None if (raw_count is None) else require_count(raw_count, "evidence count")
    )
    observation = observe_score(fields.read(row, "score"))
    eligible = (
        observation.state == "PRESENT_POSITIVE"
        and evidence_count is not None
        and (evidence_count >= 1)
    )
    reason = None
    if (not eligible):
        reason = (
            "MISSING_EVIDENCE_COUNT"
            if (evidence_count is None)
            else "NO_EVIDENCE"
            if (observation.state == "PRESENT_POSITIVE")
            else observation.state
        )
    return AssociationRecord(
        "OpenTargets",
        release,
        schema.family,
        disease,
        target,
        selector,
        observation,
        evidence_count,
        eligible,
        reason,
    )

@dataclass(frozen = True)
class EntityRecord:
    resource: str
    release: str
    family: str
    entity_id: str
    label: str | None

def require_label(value: object) -> str | None:
    if (value is not None and (not isinstance(value, str))):
        raise ContractError("display label must be a string or null")
    return value

def adapt_ot_entity(
    row: Mapping[str, object], release: str, family: str, fields: FieldMap
) -> EntityRecord:
    require_text(release, "release")
    if (family not in {"disease", "target"}):
        raise ContractError("unknown Open Targets entity family")
    fields.require("entity_id", "label")
    entity_id = require_text(fields.read(row, "entity_id"), "entity_id")
    label = require_label(fields.read(row, "label"))
    if (family == "target" and re.fullmatch("ENSG[0-9]{11}", entity_id) is None):
        raise ContractError("target entity requires an Ensembl gene identifier")
    return EntityRecord("OpenTargets", release, family, entity_id, label)

def string_identifier(value: object) -> str:
    identifier = require_text(value, "STRING identifier")
    if (re.fullmatch("9606\\.[A-Za-z0-9_]+", identifier) is None):
        raise ContractError("STRING identifier must retain the human taxonomic prefix")
    return identifier

def string_score(value: object) -> object:
    if (isinstance(value, str)):
        if (re.fullmatch("[0-9]+", value) is None):
            raise ContractError("STRING score text must contain decimal integer digits")
        return int(value)
    return value

@dataclass(frozen = True)
class StringEdgeRecord:
    release: str
    protein1: str
    protein2: str
    observation: ScoreObservation

def adapt_string_edge(
    row: Mapping[str, object], release: str, fields: FieldMap, *, threshold: int = 150
) -> StringEdgeRecord:
    require_text(release, "release")
    fields.require("protein1", "protein2", "score")
    require_count(threshold, "threshold")
    observation = observe_score(
        string_score(fields.read(row, "score")),
        upper = 1000,
        threshold = threshold,
        integral = True,
    )
    return StringEdgeRecord(
        release,
        string_identifier(fields.read(row, "protein1")),
        string_identifier(fields.read(row, "protein2")),
        observation,
    )

def adapt_string_entity(
    row: Mapping[str, object], release: str, fields: FieldMap
) -> EntityRecord:
    require_text(release, "release")
    fields.require("entity_id", "label")
    label = require_label(fields.read(row, "label"))
    return EntityRecord(
        "STRING",
        release,
        "protein",
        string_identifier(fields.read(row, "entity_id")),
        label,
    )

@dataclass(frozen = True)
class CanonicalEdges:
    release: str
    edges: tuple[StringEdgeRecord, ...]
    input_rows: int
    reciprocal_rows: int
    self_loop_rows: int

def canonicalise_string_edges(
    rows: Iterable[StringEdgeRecord], release: str
) -> CanonicalEdges:
    require_text(release, "release")
    seen_directions = set()
    canonical = {}
    input_rows = 0
    reciprocal_rows = 0
    self_loops = 0
    for (row) in (rows):
        input_rows += 1
        if (not isinstance(row, StringEdgeRecord) or row.release != release):
            raise ContractError(
                "canonicalisation requires records from one declared release"
            )
        (source, target) = (
            string_identifier(row.protein1),
            string_identifier(row.protein2),
        )
        direction = (source, target)
        if (direction in seen_directions):
            raise ContractError("repeated same-direction STRING record")
        seen_directions.add(direction)
        if (source == target):
            self_loops += 1
            continue
        pair = tuple(sorted(direction))
        if (pair in canonical):
            if (canonical[pair].observation != row.observation):
                raise ContractError("conflicting reciprocal STRING records")
            reciprocal_rows += 1
        else:
            canonical[pair] = StringEdgeRecord(
                release, pair[0], pair[1], row.observation
            )
    return CanonicalEdges(
        release,
        tuple((canonical[pair] for (pair) in (sorted(canonical)))),
        input_rows,
        reciprocal_rows,
        self_loops,
    )
