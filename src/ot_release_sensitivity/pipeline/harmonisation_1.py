from __future__ import annotations
from collections import defaultdict
from pathlib import Path
import csv
import json
import re
from .harmonisation_types import ACCEPTED_DISEASE_STATUSES

def load_configuration(path):
    return json.loads(Path(path).read_text(encoding = "utf-8"))

def truth(value):
    if (isinstance(value, bool)):
        return value
    return str(value).strip().casefold() == "true"

def truth_text(value):
    return "true" if (bool(value)) else "false"

def semicolon(values):
    return ";".join(sorted({str(value) for (value) in (values) if (str(value))}))

def parse_semicolon(value):
    return tuple((part for (part) in (str(value or "").split(";")) if (part)))

def flatten_values(value):
    if (value is None):
        return ()
    if (isinstance(value, dict)):
        result = []
        for (key) in (sorted(value)):
            result.extend(flatten_values(value[key]))
        return tuple(result)
    if (isinstance(value, (list, tuple, set))):
        result = []
        for (item) in (value):
            result.extend(flatten_values(item))
        return tuple(result)
    return (str(value),)

def normalise_reference_forms(value):
    raw = str(value or "").strip()
    if (not raw):
        return set()
    forms = {raw}
    if (":" in raw):
        prefix, suffix = raw.split(":", 1)
        canonical_prefix = "Orphanet" if (prefix.upper() == "ORPHA") else prefix
        forms.add(f"{canonical_prefix}_{suffix}")
    tail = raw.rsplit("/", 1)[-1]
    if (tail):
        forms.add(tail)
    return forms

def normalise_label(value):
    return " ".join(re.sub("[^0-9a-z]+", " ", str(value or "").casefold()).split())

def namespace(identifier):
    return str(identifier).split("_", 1)[0]

def write_csv(path, rows, fieldnames):
    destination = Path(path)
    destination.parent.mkdir(parents = True, exist_ok = True)
    partial = destination.with_suffix(destination.suffix + ".partial")
    with partial.open("w", encoding = "utf-8", newline = "") as handle:
        writer = csv.DictWriter(handle, fieldnames = fieldnames, lineterminator = "\n")
        writer.writeheader()
        writer.writerows(rows)
    partial.replace(destination)

def classify_disease(record, disease_area_ids, rules):
    record_namespace = namespace(record["id"])
    if (record_namespace not in set(rules["allowed_primary_namespaces"])):
        return (False, "NON_PRIMARY_ONTOLOGY_ENTITY", "namespace_not_prespecified")
    if (bool(record.get("isTherapeuticArea"))):
        return (False, "THERAPEUTIC_AREA_RECORD", "therapeutic_area_record")
    if (not set(record.get("therapeuticAreas", ())) & set(disease_area_ids)):
        return (
            False,
            "NON_DISEASE_OR_UNASSIGNED_ENTITY",
            "no_explicit_disease_therapeutic_area",
        )
    if (record_namespace in set(rules["authority_required_for_namespaces"])):
        prefixes = tuple(
            (prefix + ":" for (prefix) in (rules["accepted_authority_prefixes"]))
        )
        references = flatten_values(record.get("dbXRefs"))
        if (not any((reference.startswith(prefixes) for (reference) in (references)))):
            return (
                False,
                "EFO_WITHOUT_DISEASE_AUTHORITY",
                "no_explicit_disease_authority_cross_reference",
            )
    return (True, "PRESPECIFIED_DISEASE_ENTITY", "")

def build_disease_pair_mapping(old_records, new_records, old_release, new_release):
    old_ids = set(old_records)
    new_ids = set(new_records)
    replacement_edges = defaultdict(set)
    replacement_terms = defaultdict(set)
    cross_reference_edges = defaultdict(set)
    cross_reference_terms = defaultdict(set)
    for (old_id, record) in (old_records.items()):
        for (value) in (flatten_values(record.get("dbXRefs"))):
            matches = normalise_reference_forms(value) & new_ids
            for (new_id) in (matches):
                if (new_id != old_id):
                    cross_reference_edges[old_id].add(new_id)
                    cross_reference_terms[old_id, new_id].add(value)
    for (new_id, record) in (new_records.items()):
        for (field) in (("obsoleteTerms", "obsoleteXRefs")):
            for (value) in (flatten_values(record.get(field))):
                matches = normalise_reference_forms(value) & old_ids
                for (old_id) in (matches):
                    if (old_id != new_id):
                        replacement_edges[old_id].add(new_id)
                        replacement_terms[old_id, new_id].add(f"{field}:{value}")
        for (value) in (flatten_values(record.get("dbXRefs"))):
            matches = normalise_reference_forms(value) & old_ids
            for (old_id) in (matches):
                if (old_id != new_id):
                    cross_reference_edges[old_id].add(new_id)
                    cross_reference_terms[old_id, new_id].add(value)
    claims = {}
    evidence_kind = {}
    for (old_id) in (old_ids):
        if (old_id in new_ids):
            claims[old_id] = {old_id}
            evidence_kind[old_id] = "IDENTITY"
        elif (replacement_edges.get(old_id)):
            claims[old_id] = set(replacement_edges[old_id])
            evidence_kind[old_id] = "REPLACEMENT"
        else:
            claims[old_id] = set(cross_reference_edges.get(old_id, set()))
            evidence_kind[old_id] = "CROSS_REFERENCE" if (claims[old_id]) else "NONE"
    reverse_claims = defaultdict(set)
    for (old_id, candidates) in (claims.items()):
        for (new_id) in (candidates):
            reverse_claims[new_id].add(old_id)
    label_index = defaultdict(set)
    for (new_id, record) in (new_records.items()):
        label = normalise_label(record.get("name"))
        if (label):
            label_index[label].add(new_id)
    rows = []
    accepted = {}
    for (old_id) in (sorted(old_records)):
        candidates = sorted(claims[old_id])
        label_candidates = []
        mapped_id = ""
        if (old_id in new_ids):
            status = "EXACT_ID"
            method = "exact_identifier"
            source = "old_and_new_release_disease_metadata"
            mapped_id = old_id
        elif (len(candidates) == 1 and len(reverse_claims[candidates[0]]) == 1):
            mapped_id = candidates[0]
            if (evidence_kind[old_id] == "REPLACEMENT"):
                status = "EXPLICIT_REPLACEMENT"
                method = "forward_explicit_replacement"
                source = "new_release_obsoleteTerms_or_obsoleteXRefs"
            else:
                status = "EXPLICIT_CROSS_REFERENCE"
                method = "direct_one_to_one_cross_reference"
                source = "release_specific_dbXRefs"
        elif (len(candidates) > 1):
            status = "AMBIGUOUS_ONE_TO_MANY"
            method = "excluded_ambiguous_mapping"
            source = "release_specific_mapping_metadata"
        elif (len(candidates) == 1 and len(reverse_claims[candidates[0]]) > 1):
            status = "AMBIGUOUS_MANY_TO_ONE"
            method = "excluded_ambiguous_mapping"
            source = "release_specific_mapping_metadata"
        else:
            label = normalise_label(old_records[old_id].get("name"))
            label_candidates = sorted(label_index.get(label, set())) if (label) else []
            status = "LABEL_ONLY_CANDIDATE" if (label_candidates) else "UNRESOLVED"
            method = (
                "diagnostic_label_candidate"
                if (label_candidates)
                else "no_mapping_found"
            )
            source = (
                "normalised_label_diagnostic"
                if (label_candidates)
                else "release_specific_mapping_metadata"
            )
        if (status in ACCEPTED_DISEASE_STATUSES):
            accepted[old_id] = mapped_id
        selected_pairs = [(old_id, candidate) for (candidate) in (candidates)]
        replacement_evidence = set()
        cross_reference_evidence = set()
        for (pair) in (selected_pairs):
            replacement_evidence.update(replacement_terms.get(pair, set()))
            cross_reference_evidence.update(cross_reference_terms.get(pair, set()))
        rows.append(
            {
                "old_release": old_release,
                "new_release": new_release,
                "old_id": old_id,
                "old_label": old_records[old_id].get("name", ""),
                "mapped_id": mapped_id,
                "identifier_mapping_status": status,
                "mapping_method": method,
                "mapping_source": source,
                "ambiguity_status": ambiguity_status(status),
                "candidate_ids": semicolon(candidates),
                "replacement_or_obsolete_term_evidence": semicolon(
                    replacement_evidence
                ),
                "cross_reference_evidence": semicolon(cross_reference_evidence),
                "label_only_candidates": semicolon(label_candidates),
            }
        )
    return (rows, accepted)

def ambiguity_status(status):
    if (status in ACCEPTED_DISEASE_STATUSES):
        return "UNAMBIGUOUS"
    if (status == "AMBIGUOUS_ONE_TO_MANY"):
        return "ONE_TO_MANY"
    if (status == "AMBIGUOUS_MANY_TO_ONE"):
        return "MANY_TO_ONE"
    if (status == "LABEL_ONLY_CANDIDATE"):
        return "LABEL_ONLY"
    if (status == "EXCLUDED_ENTITY_TYPE"):
        return "NOT_APPLICABLE_ENTITY_EXCLUSION"
    return "UNRESOLVED"

def inverse_unique(mapping):
    inverse = defaultdict(list)
    for (source_id, destination_id) in (mapping.items()):
        inverse[destination_id].append(source_id)
    return {
        destination_id: source_ids[0]
        for ((destination_id, source_ids)) in (inverse.items())
        if (len(source_ids) == 1)
    }

def incoming_mapping_rows(rows):
    result = defaultdict(list)
    for (row) in (rows):
        for (candidate) in (parse_semicolon(row["candidate_ids"])):
            result[candidate].append(row)
    return result

def mapped_area_set(
    release,
    source_area_ids,
    accepted_25_03,
    accepted_03_06,
    reference_areas,
    inverse_03_06 = None,
):
    reference = set(reference_areas)
    if (release == "25.12"):
        return {
            accepted_25_03[area]
            for (area) in (source_area_ids)
            if (area in accepted_25_03 and accepted_25_03[area] in reference)
        }
    if (release == "26.03"):
        return set(source_area_ids) & reference
    inverse = (
        inverse_03_06 if (inverse_03_06 is not None) else inverse_unique(accepted_03_06)
    )
    return {
        inverse[area]
        for (area) in (source_area_ids)
        if (area in inverse and inverse[area] in reference)
    }

def mapping_exclusion_reason(status):
    values = {
        "AMBIGUOUS_ONE_TO_MANY": "ambiguous_one_to_many_mapping",
        "AMBIGUOUS_MANY_TO_ONE": "ambiguous_many_to_one_mapping",
        "LABEL_ONLY_CANDIDATE": "label_only_mapping_forbidden",
        "UNRESOLVED": "unresolved_mapping",
        "EXCLUDED_ENTITY_TYPE": "excluded_entity_type",
    }
    return values.get(status, "")

def build_disease_crosswalk(diseases, configuration):
    releases = tuple(configuration["releases"])
    reference_release = configuration["reference_release"]
    if (releases != ("25.12", "26.03", "26.06") or reference_release != "26.03"):
        raise ValueError(
            "The implementation contract requires releases 25.12, 26.03 and 26.06 with reference release 26.03"
        )
    rules = configuration["disease_mapping"]
    pair_25_03, accepted_25_03 = build_disease_pair_mapping(
        diseases["25.12"], diseases["26.03"], "25.12", "26.03"
    )
    pair_03_06, accepted_03_06 = build_disease_pair_mapping(
        diseases["26.03"], diseases["26.06"], "26.03", "26.06"
    )
    pair_25_lookup = {row["old_id"]: row for (row) in (pair_25_03)}
    pair_03_lookup = {row["old_id"]: row for (row) in (pair_03_06)}
    incoming_06 = incoming_mapping_rows(pair_03_06)
    inverse_03_06 = inverse_unique(accepted_03_06)
    inverse_25_03 = inverse_unique(accepted_25_03)
    reference_areas = tuple(rules["reference_therapeutic_areas"])
    release_area_ids = {
        "25.12": {
            inverse_25_03[area]
            for (area) in (reference_areas)
            if (area in inverse_25_03)
        },
        "26.03": set(reference_areas),
        "26.06": {
            accepted_03_06[area]
            for (area) in (reference_areas)
            if (area in accepted_03_06)
        },
    }
    reference_labels = {
        identifier: record.get("name", "")
        for ((identifier, record)) in (diseases[reference_release].items())
    }
    base_rows = []
    area_membership = {}
    for (release) in (releases):
        for (source_id) in (sorted(diseases[release])):
            record = diseases[release][source_id]
            entity_eligible, entity_type, entity_exclusion = classify_disease(
                record, release_area_ids[release], rules
            )
            if (release == "25.12"):
                mapping_row = pair_25_lookup[source_id]
                canonical_id = accepted_25_03.get(source_id, "")
            elif (release == "26.03"):
                mapping_row = {
                    "identifier_mapping_status": "EXACT_ID",
                    "mapping_method": "reference_release_identity",
                    "mapping_source": "reference_release_disease_metadata",
                    "ambiguity_status": "UNAMBIGUOUS",
                    "replacement_or_obsolete_term_evidence": "",
                    "cross_reference_evidence": "",
                    "label_only_candidates": "",
                }
                canonical_id = source_id
            else:
                canonical_id = inverse_03_06.get(source_id, "")
                if (canonical_id):
                    forward = pair_03_lookup[canonical_id]
                    mapping_row = {
                        **forward,
                        "mapping_method": "inverse_of_" + forward["mapping_method"],
                        "mapping_source": forward["mapping_source"],
                    }
                else:
                    incoming = incoming_06.get(source_id, [])
                    incoming_statuses = {
                        row["identifier_mapping_status"] for (row) in (incoming)
                    }
                    if ("AMBIGUOUS_MANY_TO_ONE" in incoming_statuses):
                        status = "AMBIGUOUS_MANY_TO_ONE"
                    elif ("AMBIGUOUS_ONE_TO_MANY" in incoming_statuses):
                        status = "AMBIGUOUS_ONE_TO_MANY"
                    else:
                        status = "UNRESOLVED"
                    mapping_row = {
                        "identifier_mapping_status": status,
                        "mapping_method": "unresolved_inverse_reference_mapping",
                        "mapping_source": "release_specific_mapping_metadata",
                        "ambiguity_status": ambiguity_status(status),
                        "replacement_or_obsolete_term_evidence": semicolon(
                            (
                                row["replacement_or_obsolete_term_evidence"]
                                for (row) in (incoming)
                            )
                        ),
                        "cross_reference_evidence": semicolon(
                            (row["cross_reference_evidence"] for (row) in (incoming))
                        ),
                        "label_only_candidates": "",
                    }
            identifier_status = mapping_row["identifier_mapping_status"]
            mapping_status = (
                identifier_status if (entity_eligible) else "EXCLUDED_ENTITY_TYPE"
            )
            canonical_areas = mapped_area_set(
                release,
                record.get("therapeuticAreas", ()),
                accepted_25_03,
                accepted_03_06,
                reference_areas,
                inverse_03_06,
            )
            area_membership[release, source_id] = canonical_areas
            primary_eligible = (
                entity_eligible
                and identifier_status in ACCEPTED_DISEASE_STATUSES
                and bool(canonical_id)
            )
            exclusion_reason = (
                entity_exclusion
                if (not entity_eligible)
                else mapping_exclusion_reason(identifier_status)
            )
            base_rows.append(
                {
                    "release": release,
                    "source_disease_id": source_id,
                    "source_label": record.get("name", ""),
                    "source_namespace": namespace(source_id),
                    "source_entity_type": entity_type,
                    "source_therapeutic_area_ids": semicolon(
                        record.get("therapeuticAreas", ())
                    ),
                    "canonical_disease_id": canonical_id,
                    "canonical_label": reference_labels.get(canonical_id, ""),
                    "canonical_therapeutic_area_ids": semicolon(canonical_areas),
                    "identifier_mapping_status": identifier_status,
                    "mapping_status": mapping_status,
                    "mapping_method": mapping_row["mapping_method"],
                    "mapping_source": mapping_row["mapping_source"],
                    "ambiguity_status": ambiguity_status(mapping_status),
                    "replacement_or_obsolete_term_evidence": mapping_row.get(
                        "replacement_or_obsolete_term_evidence", ""
                    ),
                    "cross_reference_evidence": mapping_row.get(
                        "cross_reference_evidence", ""
                    ),
                    "label_only_candidates": mapping_row.get(
                        "label_only_candidates", ""
                    ),
                    "release_native_eligible": truth_text(entity_eligible),
                    "primary_support_eligible": truth_text(primary_eligible),
                    "exclusion_reason": exclusion_reason,
                }
            )
    by_canonical = defaultdict(dict)
    for (row) in (base_rows):
        if (truth(row["primary_support_eligible"])):
            by_canonical[row["canonical_disease_id"]][row["release"]] = row
    pair_labels = (("25.12", "26.03"), ("26.03", "26.06"))
    final_rows = []
    for (row) in (base_rows):
        canonical_rows = by_canonical.get(row["canonical_disease_id"], {})
        eligible_pairs = [
            f"{old_release}_to_{new_release}"
            for ((old_release, new_release)) in (pair_labels)
            if (old_release in canonical_rows and new_release in canonical_rows)
        ]
        complete = all((release in canonical_rows for (release) in (releases)))
        area_sets = (
            [
                set(
                    parse_semicolon(
                        canonical_rows[release]["canonical_therapeutic_area_ids"]
                    )
                )
                for (release) in (releases)
            ]
            if (complete)
            else []
        )
        stable_areas = bool(
            complete and area_sets[0] and (area_sets[0] == area_sets[1] == area_sets[2])
        )
        stable_metadata = bool(
            complete
            and stable_areas
            and all(
                (
                    canonical_rows[release]["identifier_mapping_status"] == "EXACT_ID"
                    for (release) in (releases)
                )
            )
            and (
                len(
                    {
                        canonical_rows[release]["source_disease_id"]
                        for (release) in (releases)
                    }
                )
                == 1
            )
            and (
                len(
                    {
                        canonical_rows[release]["source_label"]
                        for (release) in (releases)
                    }
                )
                == 1
            )
            and (
                len(
                    {
                        canonical_rows[release]["source_entity_type"]
                        for (release) in (releases)
                    }
                )
                == 1
            )
        )
        final_rows.append(
            {
                **row,
                "pair_mapping_eligible": semicolon(eligible_pairs),
                "three_release_mapping_eligible": truth_text(complete),
                "therapeutic_area_stability": "STABLE"
                if (stable_areas)
                else "UNSTABLE"
                if (complete)
                else "NOT_COMPARABLE",
                "stable_entity_eligible": truth_text(stable_metadata),
            }
        )
    area_rows = []
    for (row) in (final_rows):
        canonical_areas = parse_semicolon(row["canonical_therapeutic_area_ids"])
        stable = row["therapeutic_area_stability"] == "STABLE"
        weight = 1 / len(canonical_areas) if (stable and canonical_areas) else ""
        source_areas = parse_semicolon(row["source_therapeutic_area_ids"])
        source_to_canonical = {}
        for (source_area) in (source_areas):
            mapped = mapped_area_set(
                row["release"],
                (source_area,),
                accepted_25_03,
                accepted_03_06,
                reference_areas,
                inverse_03_06,
            )
            for (canonical_area) in (mapped):
                source_to_canonical[source_area] = canonical_area
        for (source_area, canonical_area) in (sorted(source_to_canonical.items())):
            area_rows.append(
                {
                    "release": row["release"],
                    "source_disease_id": row["source_disease_id"],
                    "canonical_disease_id": row["canonical_disease_id"],
                    "source_therapeutic_area_id": source_area,
                    "canonical_therapeutic_area_id": canonical_area,
                    "canonical_therapeutic_area_label": reference_labels.get(
                        canonical_area, ""
                    ),
                    "therapeutic_area_stability": row["therapeutic_area_stability"],
                    "primary_area_contrast_eligible": truth_text(
                        stable and truth(row["primary_support_eligible"])
                    ),
                    "fractional_weight": weight,
                }
            )
    final_rows.sort(
        key = lambda row: (releases.index(row["release"]), row["source_disease_id"])
    )
    area_rows.sort(
        key = lambda row: (
            releases.index(row["release"]),
            row["source_disease_id"],
            row["canonical_therapeutic_area_id"],
        )
    )
    return (final_rows, area_rows, pair_25_03, pair_03_06)

def target_prior_symbol_candidates(record, reference_records):
    source_symbol = str(record.get("approvedSymbol", ""))
    if (not source_symbol):
        return ()
    candidates = []
    for (identifier, reference) in (reference_records.items()):
        values = {str(reference.get("approvedSymbol", ""))}
        values.update(flatten_values(reference.get("obsoleteSymbols")))
        if (source_symbol in values):
            candidates.append(identifier)
    return tuple(sorted(candidates))

def build_target_crosswalk(targets, configuration):
    releases = tuple(configuration["releases"])
    reference_release = configuration["reference_release"]
    reference = targets[reference_release]
    common = set.intersection(*(set(targets[release]) for (release) in (releases)))
    stable = {
        identifier
        for (identifier) in (common)
        if (
            all(
                (
                    len(
                        {
                            str(targets[release][identifier].get(field, ""))
                            for (release) in (releases)
                        }
                    )
                    == 1
                    for (field) in (("approvedSymbol", "approvedName", "biotype"))
                )
            )
        )
    }
    rows = []
    for (release) in (releases):
        for (source_id) in (sorted(targets[release])):
            record = targets[release][source_id]
            if (source_id in reference):
                canonical_id = source_id
                status = "EXACT_ID"
                method = "exact_stable_identifier"
                source = "release_and_reference_target_metadata"
                candidates = ()
                eligible = True
                exclusion_reason = ""
            else:
                canonical_id = ""
                status = "UNRESOLVED"
                method = "no_exact_stable_identifier"
                source = "release_and_reference_target_metadata"
                candidates = target_prior_symbol_candidates(record, reference)
                eligible = False
                exclusion_reason = (
                    "symbol_only_mapping_forbidden"
                    if (candidates)
                    else "unresolved_target_identifier"
                )
            rows.append(
                {
                    "release": release,
                    "source_target_id": source_id,
                    "source_approved_symbol": record.get("approvedSymbol", ""),
                    "source_approved_name": record.get("approvedName", ""),
                    "source_biotype": record.get("biotype", ""),
                    "canonical_target_id": canonical_id,
                    "canonical_symbol": reference.get(canonical_id, {}).get(
                        "approvedSymbol", ""
                    ),
                    "mapping_method": method,
                    "mapping_status": status,
                    "mapping_source": source,
                    "prior_symbol_evidence": semicolon(candidates),
                    "stable_entity_eligible": truth_text(source_id in stable),
                    "primary_support_eligible": truth_text(eligible),
                    "exclusion_reason": exclusion_reason,
                }
            )
    rows.sort(key = lambda row: (releases.index(row["release"]), row["source_target_id"]))
    return rows

def quote_identifier(value):
    if (not re.fullmatch("[A-Za-z_][A-Za-z0-9_]*", str(value))):
        raise ValueError(f"Unsafe SQL identifier: {value}")
    return '"' + str(value).replace('"', '""') + '"'

def quote_literal(value):
    return "'" + str(value).replace("'", "''") + "'"

def parquet_source(path):
    return (
        "read_parquet("
        + quote_literal(str(Path(path).resolve()).replace("\\", "/"))
        + ")"
    )

def discover_inputs(root, configuration):
    root_path = Path(root).resolve()
    grouped = {}
    for (release) in (configuration["releases"]):
        for (dataset) in (configuration["required_datasets"]):
            base = root_path / "data" / "raw" / release / dataset
            files = sorted(base.rglob("*.parquet")) if (base.is_dir()) else []
            if (not files):
                raise FileNotFoundError(
                    f"No Parquet inputs for {release}/{dataset}: {base}"
                )
            if (any(
                (not path.resolve().is_relative_to(root_path) for (path) in (files))
            )):
                raise ValueError("Input path escapes the working directory")
            grouped[release, dataset] = [
                (path.relative_to(root_path).as_posix(), path) for (path) in (files)
            ]
    return grouped

def load_disease_records(connection, files):
    source = (
        "read_parquet(["
        + ",".join(
            (quote_literal(str(path).replace("\\", "/")) for ((_, path)) in (files))
        )
        + "])"
    )
    rows = connection.execute(
        f"SELECT id, code, name, dbXRefs, parents, obsoleteTerms, obsoleteXRefs, ancestors, therapeuticAreas, ontology.isTherapeuticArea FROM {source} ORDER BY id"
    ).fetchall()
    records = {}
    for (row) in (rows):
        if (row[0] in records):
            raise RuntimeError(f"Duplicate disease identifier: {row[0]}")
        records[row[0]] = {
            "id": row[0],
            "code": row[1] or "",
            "name": row[2] or "",
            "dbXRefs": tuple(flatten_values(row[3])),
            "parents": tuple(flatten_values(row[4])),
            "obsoleteTerms": tuple(flatten_values(row[5])),
            "obsoleteXRefs": tuple(flatten_values(row[6])),
            "ancestors": tuple(flatten_values(row[7])),
            "therapeuticAreas": tuple(flatten_values(row[8])),
            "isTherapeuticArea": bool(row[9]),
        }
    return records

def load_target_records(connection, files):
    source = (
        "read_parquet(["
        + ",".join(
            (quote_literal(str(path).replace("\\", "/")) for ((_, path)) in (files))
        )
        + "])"
    )
    rows = connection.execute(
        f"SELECT id, approvedSymbol, approvedName, biotype, alternativeGenes, obsoleteSymbols, obsoleteNames FROM {source} ORDER BY id"
    ).fetchall()
    records = {}
    for (row) in (rows):
        if (row[0] in records):
            raise RuntimeError(f"Duplicate target identifier: {row[0]}")
        records[row[0]] = {
            "id": row[0],
            "approvedSymbol": row[1] or "",
            "approvedName": row[2] or "",
            "biotype": row[3] or "",
            "alternativeGenes": tuple(flatten_values(row[4])),
            "obsoleteSymbols": tuple(flatten_values(row[5])),
            "obsoleteNames": tuple(flatten_values(row[6])),
        }
    return records

def projection_sql(files, release, dataset, field_contract, endpoint_id = None):
    parts = []
    disease_field = quote_identifier(field_contract["disease_id"])
    target_field = quote_identifier(field_contract["target_id"])
    score_field = quote_identifier(field_contract["score"])
    evidence_field = quote_identifier(field_contract["evidence_count"])
    filters = []
    if ("filter_field" in field_contract):
        filters.append(
            f"{quote_identifier(field_contract['filter_field'])} = {quote_literal(field_contract['filter_value'])}"
        )
    if (endpoint_id is not None):
        filters.append(
            f"{quote_identifier(field_contract['endpoint_field'])} = {quote_literal(endpoint_id)}"
        )
    where = " WHERE " + " AND ".join(filters) if (filters) else ""
    for (relative, absolute) in (files):
        source = parquet_source(absolute)
        parts.append(
            "SELECT "
            + quote_literal(release)
            + " AS release, "
            + f"CAST({disease_field} AS VARCHAR) AS source_disease_id, "
            + f"CAST({target_field} AS VARCHAR) AS source_target_id, "
            + f"CAST({score_field} AS DOUBLE) AS score, "
            + f"CAST({evidence_field} AS BIGINT) AS evidence_count, "
            + quote_literal(dataset)
            + " AS score_source_dataset, "
            + quote_literal(relative)
            + " AS source_object_relative_path "
            + f"FROM {source}{where}"
        )
    return " UNION ALL ".join(parts)

def required_assertion_sql(files, field_contract):
    if ("assert_field" not in field_contract):
        return "SELECT 0 AS invalid_count"
    filter_field = quote_identifier(field_contract["filter_field"])
    filter_value = quote_literal(field_contract["filter_value"])
    assert_field = quote_identifier(field_contract["assert_field"])
    assert_value = quote_literal(field_contract["assert_value"])
    parts = [
        f"SELECT COUNT(*) AS invalid_count FROM {parquet_source(absolute)} WHERE {filter_field} = {filter_value} AND ({assert_field} IS NULL OR CAST({assert_field} AS VARCHAR) <> {assert_value})"
        for ((_, absolute)) in (files)
    ]
    return (
        "SELECT SUM(invalid_count) AS invalid_count FROM ("
        + " UNION ALL ".join(parts)
        + ")"
    )
