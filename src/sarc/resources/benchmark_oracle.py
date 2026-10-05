import math
from collections import Counter

def oracle_identity(baseline_rows, followup_rows, mappings, identity_decidable):
    left = dict(baseline_rows)
    right_ids = [row[1] for (row) in (followup_rows)]
    if (len(set(right_ids)) != len(right_ids) or not identity_decidable):
        return (left, None, False)
    right = {member: score for ((_, member, score)) in (followup_rows)}
    exact = set(left).intersection(right)
    source_targets = {}
    target_sources = {}
    for (source, target) in (mappings):
        if (
            source not in left
            or target not in right
            or source in exact
            or (target in exact)
        ):
            return (left, None, False)
        source_targets.setdefault(source, set()).add(target)
        target_sources.setdefault(target, set()).add(source)
    if (any((len(values) != 1 for (values) in (source_targets.values()))) or any(
        (len(values) != 1 for (values) in (target_sources.values()))
    )):
        return (left, None, False)
    replacements = {target: source for ((source, target)) in (mappings)}
    right = {
        replacements.get(member, member): score for ((member, score)) in (right.items())
    }
    return (left, right, True)

def oracle_set(left, right, top_k, minimum):
    if (len(left) < max(top_k, minimum) or len(right) < max(top_k, minimum)):
        return {"estimable": False, "value": None, "support_denominator": 0}
    ranked_left = sorted(left, key = lambda member: left[member], reverse = True)
    ranked_right = sorted(right, key = lambda member: right[member], reverse = True)
    (left_value, right_value) = (
        left[ranked_left[top_k - 1]],
        right[ranked_right[top_k - 1]],
    )
    first = {member for (member) in (ranked_left) if (left[member] >= left_value)}
    second = {member for (member) in (ranked_right) if (right[member] >= right_value)}
    return {
        "estimable": True,
        "value": first != second,
        "support_denominator": 1,
        "baseline_set": sorted(first),
        "followup_set": sorted(second),
        "jaccard_numerator": len(first.intersection(second)),
        "jaccard_denominator": len(first.union(second)),
    }

def tie_pair_count(values):
    return sum((count * (count - 1) // 2 for (count) in (Counter(values).values())))

def oracle_tau(first, second):
    if (len(set(first)) < 2 or len(set(second)) < 2):
        return None
    pairs = sorted(zip(first, second))
    positions = {
        value: index + 1 for ((index, value)) in (enumerate(sorted(set(second))))
    }
    tree = [0] * (len(positions) + 1)
    discordant = 0
    inserted = 0
    offset = 0
    while (offset < len(pairs)):
        stop = offset + 1
        while (stop < len(pairs) and pairs[stop][0] == pairs[offset][0]):
            stop += 1
        for (_, value) in (pairs[offset:stop]):
            position = positions[value]
            smaller_or_equal = 0
            while (position > 0):
                smaller_or_equal += tree[position]
                position -= position & -position
            discordant += inserted - smaller_or_equal
        for (_, value) in (pairs[offset:stop]):
            position = positions[value]
            while (position < len(tree)):
                tree[position] += 1
                position += position & -position
            inserted += 1
        offset = stop
    pair_count = len(first) * (len(first) - 1) // 2
    (first_ties, second_ties) = (tie_pair_count(first), tie_pair_count(second))
    both_ties = tie_pair_count(list(zip(first, second)))
    numerator = pair_count - first_ties - second_ties + both_ties - 2 * discordant
    return numerator / math.sqrt((pair_count - first_ties) * (pair_count - second_ties))

def average_ranks(values):
    ordered = sorted(range(len(values)), key = lambda index: values[index])
    ranks = [0.0] * len(values)
    start = 0
    while (start < len(values)):
        stop = start + 1
        while (stop < len(values) and values[ordered[start]] == values[ordered[stop]]):
            stop += 1
        rank = (start + 1 + stop) / 2
        for (index) in (ordered[start:stop]):
            ranks[index] = rank
        start = stop
    return ranks

def oracle_spearman(first, second):
    if (len(set(first)) < 2 or len(set(second)) < 2):
        return None
    (first, second) = (average_ranks(first), average_ranks(second))
    centre = (len(first) + 1) / 2
    numerator = sum(((x - centre) * (y - centre) for ((x, y)) in (zip(first, second))))
    denominator = math.sqrt(
        sum(((x - centre) ** 2 for (x) in (first)))
        * sum(((y - centre) ** 2 for (y) in (second)))
    )
    return numerator / denominator

def expected_answers(scenario, minimum_roster = 30, top_ks = (5, 10, 20)):
    (left, right, decidable) = oracle_identity(
        scenario.baseline,
        scenario.followup,
        scenario.mappings,
        scenario.identity_decidable,
    )
    answers = {
        "mapping_ambiguity_requires_abstention": {
            "estimable": True,
            "value": not decidable,
            "support_denominator": 1,
        }
    }
    if (not decidable):
        absent = {"estimable": False, "value": None, "support_denominator": 0}
        answers.update(
            {
                name: dict(absent)
                for (name) in (
                    (
                        "fixed_leading_set_change",
                        "native_leading_set_change",
                        "persistent_versus_roster_membership_description",
                    )
                )
            }
        )
        return {
            "tasks": answers,
            "identity_decidable": False,
            "fixed_count": None,
            "top_k": {},
            "rank": None,
        }
    persistent = set(left).intersection(right)
    fixed_left = {member: left[member] for (member) in (persistent)}
    fixed_right = {member: right[member] for (member) in (persistent)}
    answers["fixed_leading_set_change"] = oracle_set(
        fixed_left, fixed_right, 1, minimum_roster
    )
    answers["native_leading_set_change"] = oracle_set(left, right, 1, 1)
    answers["persistent_versus_roster_membership_description"] = {
        "estimable": True,
        "value": {
            "persistent": sorted(persistent),
            "entering": sorted(set(right).difference(left)),
            "exiting": sorted(set(left).difference(right)),
        },
        "support_denominator": 1,
    }
    top_k = {
        str(count): {
            "native": oracle_set(left, right, count, 1),
            "fixed": oracle_set(fixed_left, fixed_right, count, minimum_roster),
        }
        for (count) in (top_ks)
    }
    members = sorted(persistent)
    (first, second) = (
        [left[member] for (member) in (members)],
        [right[member] for (member) in (members)],
    )
    rank = (
        None
        if (len(members) < minimum_roster)
        else {
            "kendall_tau_b": oracle_tau(first, second),
            "spearman": oracle_spearman(first, second),
        }
    )
    return {
        "tasks": answers,
        "identity_decidable": True,
        "fixed_count": len(persistent),
        "top_k": top_k,
        "rank": rank,
    }
