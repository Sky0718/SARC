from dataclasses import dataclass
from decimal import Decimal
from math import isclose, isfinite
from typing import Iterable

class ContractError(ValueError):
    pass

def require_text(value: str, name: str) -> str:
    if (not isinstance(value, str) or not value or value != value.strip()):
        raise ContractError(f"{name} must be a nonempty, unpadded string")
    return value

def require_number(value: object, name: str) -> float:
    if (isinstance(value, bool) or not isinstance(value, (int, float, Decimal))):
        raise ContractError(f"{name} must be numeric")
    try:
        result = float(value)
    except (OverflowError, ValueError):
        raise ContractError(f"{name} must be finite") from None
    if (not isfinite(result)):
        raise ContractError(f"{name} must be finite")
    return result

def require_count(value: object, name: str) -> int:
    if (isinstance(value, bool) or not isinstance(value, int) or value < 0):
        raise ContractError(f"{name} must be a nonnegative integer")
    return value

@dataclass(frozen = True)
class ResultContext:
    resource: str
    releases: tuple[str, ...]
    cohort: str
    support: str
    estimand: str
    unit: str
    source: str
    configuration: str
    panel_id: str | None = None

    def __post_init__(self):
        for (name) in (("resource", "cohort", "support", "estimand", "unit", "source")):
            require_text(getattr(self, name), name)
        if (not isinstance(self.releases, tuple) or not self.releases):
            raise ContractError("releases must be a nonempty tuple")
        for (release) in (self.releases):
            require_text(release, "release")
        require_text(self.configuration, "configuration")
        if (self.panel_id is not None):
            require_text(self.panel_id, "panel_id")

@dataclass(frozen = True)
class Measure:
    context: ResultContext
    kind: str
    value: object
    support_denominator: int | float
    numerator: int | None = None
    denominator: int | None = None
    reason: str | None = None
    tolerance: float = 1e-12

    def __post_init__(self):
        if (not isinstance(self.context, ResultContext)):
            raise ContractError("context must be a ResultContext")
        if (self.kind not in {
            "fraction",
            "weighted_count",
            "count",
            "scalar",
            "boolean",
            "category",
        }):
            raise ContractError("unknown measurement kind")
        support = require_number(self.support_denominator, "support_denominator")
        tolerance = require_number(self.tolerance, "tolerance")
        if (support < 0 or not 0 <= tolerance <= 1e-12):
            raise ContractError(
                "support must be nonnegative and tolerance cannot exceed the frozen 1e-12"
            )
        if (self.kind == "fraction"):
            self._validate_fraction()
            return
        if (self.numerator is not None or self.denominator is not None):
            raise ContractError("non-fraction measures cannot carry fraction fields")
        if (self.value is None):
            if (support != 0):
                raise ContractError(
                    "unavailable measures require zero eligible support"
                )
            require_text(self.reason, "unavailable reason")
            return
        if (self.reason is not None):
            raise ContractError(
                "an available measure cannot have an unavailable reason"
            )
        if (self.kind == "category"):
            require_text(self.value, "category value")
        elif (self.kind == "boolean"):
            if (not isinstance(self.value, bool)):
                raise ContractError("boolean value must be bool")
        elif (self.kind == "count"):
            require_count(self.value, "count value")
        else:
            value = require_number(self.value, "value")
            if (self.kind == "weighted_count" and value < 0):
                raise ContractError("weighted counts cannot be negative")

    def _validate_fraction(self):
        require_count(self.support_denominator, "fraction support_denominator")
        fields = (self.value, self.numerator, self.denominator)
        if (all((value is None for (value) in (fields)))):
            if (self.support_denominator != 0):
                raise ContractError("unavailable fraction has nonzero support")
            require_text(self.reason, "unavailable fraction reason")
            return
        if (any((value is None for (value) in (fields)))):
            raise ContractError("incomplete fraction fields")
        if (self.reason is not None):
            raise ContractError(
                "an estimable fraction cannot have an unavailable reason"
            )
        numerator = require_count(self.numerator, "numerator")
        denominator = require_count(self.denominator, "denominator")
        value = require_number(self.value, "fraction value")
        if (denominator <= 0 or numerator > denominator or (not 0 <= value <= 1)):
            raise ContractError("fraction is outside its closed domain")
        if (self.support_denominator != denominator):
            raise ContractError("fraction denominator does not match eligible support")
        if (not isclose(
            value, numerator / denominator, rel_tol = 0, abs_tol = self.tolerance
        )):
            raise ContractError("fraction value does not match numerator / denominator")

@dataclass(frozen = True)
class ScoreObservation:
    raw_score: float | None
    state: str
    passes_threshold: bool | None

    def __post_init__(self):
        if (self.passes_threshold is not None and (
            not isinstance(self.passes_threshold, bool)
        )):
            raise ContractError("passes_threshold must be bool or null")
        if (self.state in {"ABSENT_OR_CENSORED", "PRESENT_NULL"}):
            if (self.raw_score is not None or self.passes_threshold is not None):
                raise ContractError(
                    "unobserved scores cannot carry values or threshold results"
                )
        elif (self.state in {"PRESENT_ZERO", "PRESENT_POSITIVE"}):
            value = require_number(self.raw_score, "raw_score")
            if (
                self.state == "PRESENT_ZERO"
                and value != 0
                or (self.state == "PRESENT_POSITIVE" and value <= 0)
            ):
                raise ContractError("score state contradicts its observed value")
        else:
            raise ContractError("unknown score observation state")

def observe_score(
    value: object,
    *,
    present: bool = True,
    lower: float = 0,
    upper: float = 1,
    threshold: float | None = None,
    integral: bool = False,
) -> ScoreObservation:
    if (not isinstance(present, bool) or not isinstance(integral, bool)):
        raise ContractError("present and integral must be bool")
    low = require_number(lower, "lower bound")
    high = require_number(upper, "upper bound")
    if (low > high):
        raise ContractError("score bounds are reversed")
    cutoff = None if (threshold is None) else require_number(threshold, "threshold")
    if (cutoff is not None and (not low <= cutoff <= high)):
        raise ContractError("threshold is outside score bounds")
    if (not present):
        if (value is not None):
            raise ContractError("an absent record cannot carry a score")
        return ScoreObservation(None, "ABSENT_OR_CENSORED", None)
    if (value is None):
        return ScoreObservation(None, "PRESENT_NULL", None)
    score = require_number(value, "score")
    if (not low <= score <= high or (integral and (not score.is_integer()))):
        raise ContractError("score is outside the declared domain")
    state = "PRESENT_ZERO" if (score == 0) else "PRESENT_POSITIVE"
    return ScoreObservation(score, state, None if (cutoff is None) else score >= cutoff)

@dataclass(frozen = True)
class MappingProposal:
    baseline_id: str
    followup_id: str
    kind: str

    def __post_init__(self):
        require_text(self.baseline_id, "baseline_id")
        require_text(self.followup_id, "followup_id")
        if (self.kind not in {
            "EXPLICIT_ONE_TO_ONE_REPLACEMENT",
            "EXPLICIT_ONE_TO_ONE_CROSS_REFERENCE",
        }):
            raise ContractError(
                "mapping proposals require an explicit supported mapping kind"
            )

@dataclass(frozen = True)
class IdentityAlignment:
    matches: tuple[tuple[str, str, str], ...]
    unmatched_baseline: frozenset[str]
    unmatched_followup: frozenset[str]
    ambiguous_baseline: frozenset[str]
    ambiguous_followup: frozenset[str]

def unique_identifiers(values: Iterable[str], name: str) -> frozenset[str]:
    ordered = tuple(values)
    for (value) in (ordered):
        require_text(value, name)
    if (len(set(ordered)) != len(ordered)):
        raise ContractError(f"duplicate {name}")
    return frozenset(ordered)

def align_identities(
    baseline: Iterable[str],
    followup: Iterable[str],
    proposals: Iterable[MappingProposal] = (),
) -> IdentityAlignment:
    left = unique_identifiers(baseline, "baseline identifier")
    right = unique_identifiers(followup, "followup identifier")
    exact = left & right
    candidates = {}
    for (proposal) in (proposals):
        if (not isinstance(proposal, MappingProposal)):
            raise ContractError("invalid mapping proposal")
        (source, target) = (proposal.baseline_id, proposal.followup_id)
        if (source not in left or target not in right):
            raise ContractError(
                "mapping endpoint is outside its declared identity universe"
            )
        if (source in exact or target in exact):
            if (source == target):
                continue
            raise ContractError("mapping collides with a reserved exact match")
        candidates.setdefault(source, {}).setdefault(proposal.kind, set()).add(target)
    selected = {}
    for (source, by_kind) in (candidates.items()):
        kind = (
            "EXPLICIT_ONE_TO_ONE_REPLACEMENT"
            if ("EXPLICIT_ONE_TO_ONE_REPLACEMENT" in by_kind)
            else "EXPLICIT_ONE_TO_ONE_CROSS_REFERENCE"
        )
        selected[source] = (kind, by_kind[kind])
    incoming = {}
    for (source, (_, targets)) in (selected.items()):
        for (target) in (targets):
            incoming.setdefault(target, set()).add(source)
    bad_left = {
        source for ((source, (_, targets))) in (selected.items()) if (len(targets) != 1)
    }
    bad_right = {
        target for ((target, sources)) in (incoming.items()) if (len(sources) != 1)
    }
    bad_left.update(
        (source for (target) in (bad_right) for (source) in (incoming[target]))
    )
    bad_right.update(
        (target for (source) in (bad_left) for (target) in (selected[source][1]))
    )
    matches = [(value, value, "EXACT_ID") for (value) in (exact)]
    for (source, (kind, targets)) in (selected.items()):
        if (source not in bad_left):
            target = next(iter(targets))
            if (target not in bad_right):
                matches.append((source, target, kind))
    matched_left = {source for ((source, _, _)) in (matches)}
    matched_right = {target for ((_, target, _)) in (matches)}
    return IdentityAlignment(
        tuple(sorted(matches)),
        left - matched_left,
        right - matched_right,
        frozenset(bad_left),
        frozenset(bad_right),
    )

@dataclass(frozen = True)
class MembershipPartition:
    persistent: frozenset[str]
    entering: frozenset[str]
    exiting: frozenset[str]

def partition_membership(
    baseline: Iterable[str], followup: Iterable[str], *, identity_decidable: bool = True
) -> MembershipPartition:
    if (identity_decidable is not True):
        raise ContractError("membership requires resolved identity")
    left = unique_identifiers(baseline, "baseline member")
    right = unique_identifiers(followup, "followup member")
    return MembershipPartition(left & right, right - left, left - right)

@dataclass(frozen = True)
class PanelSnapshot:
    resource: str
    release: str
    panel_id: str
    members: tuple[tuple[str, float], ...]
    identity_decidable: bool = True

    def __post_init__(self):
        if (self.resource not in {"OpenTargets", "STRING"}):
            raise ContractError("unsupported panel resource")
        require_text(self.release, "release")
        require_text(self.panel_id, "panel_id")
        if (not isinstance(self.members, tuple) or not isinstance(
            self.identity_decidable, bool
        )):
            raise ContractError(
                "members must be a tuple and identity_decidable must be bool"
            )
        unique_identifiers((member for ((member, _)) in (self.members)), "panel member")
        for (_, value) in (self.members):
            score = observe_score(
                value,
                upper = 1000 if (self.resource == "STRING") else 1,
                integral = self.resource == "STRING",
            )
            if (score.state != "PRESENT_POSITIVE"):
                raise ContractError(
                    "panel members must already be eligible observed positive records"
                )

@dataclass(frozen = True)
class Eligibility:
    estimable: bool
    reason: str | None
    observed_members: int

def panel_eligibility(
    panel: PanelSnapshot, endpoint: str, *, minimum_roster: int = 30, top_k: int = 10
) -> Eligibility:
    if (
        require_count(minimum_roster, "minimum_roster") == 0
        or require_count(top_k, "top_k") == 0
    ):
        raise ContractError("roster minimum and top_k must be positive")
    minimums = {
        "baseline_cohort": minimum_roster,
        "native_leader": 1,
        "native_top_k": top_k,
        "fixed_leader": minimum_roster,
        "fixed_top_k": max(top_k, minimum_roster),
        "fixed_rank_correlation": minimum_roster,
        "fixed_score_revision": 1,
    }
    if (endpoint not in minimums):
        raise ContractError("unknown panel endpoint")
    count = len(panel.members)
    if (not panel.identity_decidable):
        return Eligibility(False, "IDENTITY_UNRESOLVED", count)
    if (count < minimums[endpoint]):
        return Eligibility(
            False, "EMPTY_SUPPORT" if (count == 0) else "INSUFFICIENT_MEMBERS", count
        )
    if (
        endpoint == "fixed_rank_correlation"
        and len({value for ((_, value)) in (panel.members)}) < 2
    ):
        return Eligibility(False, "CONSTANT_SCORES", count)
    return Eligibility(True, None, count)

def leading_members(panel: PanelSnapshot) -> frozenset[str] | None:
    if (not panel_eligibility(panel, "native_leader").estimable):
        return None
    maximum = max((value for ((_, value)) in (panel.members)))
    return frozenset(
        (member for ((member, value)) in (panel.members) if (value == maximum))
    )

def top_members(panel: PanelSnapshot, top_k: int = 10) -> frozenset[str] | None:
    if (not panel_eligibility(panel, "native_top_k", top_k = top_k).estimable):
        return None
    cutoff = sorted((value for ((_, value)) in (panel.members)), reverse = True)[top_k - 1]
    return frozenset(
        (member for ((member, value)) in (panel.members) if (value >= cutoff))
    )
