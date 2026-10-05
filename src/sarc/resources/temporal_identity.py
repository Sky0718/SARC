import json
import re
from collections import defaultdict
from .contracts import MappingProposal, align_identities

DISEASE_FIELDS = (
    "id",
    "name",
    "dbXRefs",
    "parents",
    "ancestors",
    "therapeuticAreas",
    "ontology",
    "obsoleteTerms",
    "obsoleteXRefs",
)
TARGET_FIELDS = (
    "id",
    "approvedSymbol",
    "approvedName",
    "biotype",
    "obsoleteSymbols",
    "obsoleteNames",
)

def reference_identifiers(value):
    if (not isinstance(value, str) or not value):
        return set()
    raw = value.strip()
    tail = raw.rsplit("/", 1)[-1]
    values = {raw, tail}
    for (item) in (tuple(values)):
        if (re.fullmatch("[A-Za-z][A-Za-z0-9]*:[A-Za-z0-9_]+", item)):
            (prefix, suffix) = item.split(":", 1)
            values.add(
                f"{('Orphanet' if (prefix.upper() == 'ORPHA') else prefix)}_{suffix}"
            )
    return values

def comparable_value(value):
    if (isinstance(value, list)):
        items = [comparable_value(item) for (item) in (value)]
        return sorted(
            items, key = lambda item: json.dumps(item, sort_keys = True, ensure_ascii = False)
        )
    if (isinstance(value, dict)):
        return {key: comparable_value(item) for ((key, item)) in (sorted(value.items()))}
    return value

def value_digest(value):
    return json.dumps(
        comparable_value(value),
        sort_keys = True,
        separators = (",", ":"),
        ensure_ascii = False,
    )

def disease_alignment(baseline, followup):
    left = set(baseline)
    right = set(followup)
    exact = left & right
    candidates = defaultdict(lambda: defaultdict(set))
    evidence = defaultdict(set)
    for (source, row) in (baseline.items()):
        if (source in exact):
            continue
        for (value) in (row.get("dbXRefs") or []):
            for (target) in (reference_identifiers(value) & right):
                candidates[source]["EXPLICIT_ONE_TO_ONE_CROSS_REFERENCE"].add(target)
                evidence[source, target, "EXPLICIT_ONE_TO_ONE_CROSS_REFERENCE"].add(
                    f"baseline.dbXRefs:{value}"
                )
    for (target, row) in (followup.items()):
        for (field) in (("obsoleteTerms", "obsoleteXRefs", "dbXRefs")):
            kind = (
                "EXPLICIT_ONE_TO_ONE_CROSS_REFERENCE"
                if (field == "dbXRefs")
                else "EXPLICIT_ONE_TO_ONE_REPLACEMENT"
            )
            for (value) in (row.get(field) or []):
                for (source) in ((reference_identifiers(value) & left) - exact):
                    candidates[source][kind].add(target)
                    evidence[source, target, kind].add(f"followup.{field}:{value}")
    selected = {}
    collisions = set()
    for (source, by_kind) in (candidates.items()):
        kind = (
            "EXPLICIT_ONE_TO_ONE_REPLACEMENT"
            if (by_kind.get("EXPLICIT_ONE_TO_ONE_REPLACEMENT"))
            else "EXPLICIT_ONE_TO_ONE_CROSS_REFERENCE"
        )
        selected[source] = (kind, by_kind[kind])
        if (by_kind[kind] & exact):
            collisions.add(source)
    proposals = [
        MappingProposal(source, target, kind)
        for ((source, (kind, targets))) in (selected.items())
        for (target) in (sorted(targets))
        if (target not in exact)
    ]
    alignment = align_identities(sorted(left), sorted(right), proposals)
    matches = {
        source: (target, kind)
        for ((source, target, kind)) in (alignment.matches)
        if (source not in collisions)
    }
    rows = []
    for (source) in (sorted(left)):
        (target, kind) = matches.get(
            source,
            (
                None,
                "AMBIGUOUS"
                if (source in alignment.ambiguous_baseline or source in collisions)
                else "UNRESOLVED",
            ),
        )
        (selected_kind, targets) = selected.get(source, (None, set()))
        support = sorted(
            {
                term
                for (candidate) in (targets)
                for (term) in (evidence[source, candidate, selected_kind])
            }
        )
        changes = []
        if (target is not None):
            changes = [
                field
                for (field) in (DISEASE_FIELDS)
                if (
                    field != "id"
                    and comparable_value(baseline[source].get(field))
                    != comparable_value(followup[target].get(field))
                )
            ]
        rows.append(
            {
                "baseline_disease_id": source,
                "followup_disease_id": target,
                "mapping_state": kind,
                "mapping_reason": "COLLISION_WITH_RESERVED_EXACT_MATCH"
                if (source in collisions)
                else "NON_BIJECTIVE_OFFICIAL_RELATION"
                if (source in alignment.ambiguous_baseline)
                else "NO_SUPPORTED_CORRESPONDENCE"
                if (target is None)
                else None,
                "candidate_ids": sorted(targets),
                "mapping_evidence": support,
                "changed_metadata_fields": changes if (target is not None) else None,
            }
        )
    ledger = [
        {
            "baseline_disease_id": source,
            "followup_candidate_id": target,
            "mapping_kind": kind,
            "evidence": sorted(terms),
            "selected_precedence": selected.get(source, (None, set()))[0] == kind,
            "reserved_exact_collision": target in exact,
        }
        for (((source, target, kind), terms)) in (sorted(evidence.items()))
    ]
    if (len(
        {
            row["followup_disease_id"]
            for (row) in (rows)
            if (row["followup_disease_id"] is not None)
        }
    ) != sum((row["followup_disease_id"] is not None for (row) in (rows)))):
        raise ValueError("Disease alignment is not globally one-to-one")
    return (rows, ledger)

def target_alignment(baseline, followup):
    rows = []
    for (identifier) in (sorted(set(baseline) | set(followup))):
        left = baseline.get(identifier)
        right = followup.get(identifier)
        changes = (
            [
                field
                for (field) in (TARGET_FIELDS)
                if (
                    field != "id"
                    and comparable_value(left.get(field))
                    != comparable_value(right.get(field))
                )
            ]
            if (left is not None and right is not None)
            else None
        )
        rows.append(
            {
                "target_id": identifier,
                "baseline_present": left is not None,
                "followup_present": right is not None,
                "identity_state": "EXACT_ID"
                if (left is not None and right is not None)
                else "BASELINE_ONLY"
                if (left is not None)
                else "FOLLOWUP_ONLY",
                "changed_metadata_fields": changes,
                "alternative_genes_comparison": "NOT_COMPARABLE_FOLLOWUP_FIELD_NOT_PROVIDED",
            }
        )
    return rows
