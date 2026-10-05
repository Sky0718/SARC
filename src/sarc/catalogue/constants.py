CONTEXT = "COAD_TCGA"
METHODS = ("DawnRank", "PRODIGY", "PersonaDrive")
SEEDS = (104729, 130363, 155921, 196613, 228017)
TRANSITIONS = ("11_0_to_11_5", "11_5_to_12_0")
NATIVE = ("native_11_0", "native_11_5", "native_12_0")
PERSISTENT = tuple(
    transition + "__" + state
    for (transition) in (TRANSITIONS)
    for (state) in (("persistent_old", "persistent_new"))
)
HYBRIDS = tuple(
    transition + "__" + state
    for (transition) in (TRANSITIONS)
    for (state) in (("confidence_first", "topology_first"))
)
PARAMETERS = {"DawnRank": "mu_3", "PRODIGY": "alpha_0.05", "PersonaDrive": "original"}
REFERENCES = (
    "NCG6_primary_all",
    "NCG6_known_subgroup",
    "NCG6_candidate_subgroup",
    "CancerMine2019_secondary",
)
SCALES = (10, 1, 5, 20)
FACTORS = ("context", "method", "network_id", "parameter_id")
ENDPOINTS = tuple(
    (support, reference, k)
    for (support) in (("common", "native"))
    for (reference) in (REFERENCES)
    for (k) in (SCALES)
)

def network_states(method):
    return NATIVE + PERSISTENT + (HYBRIDS if (method == "PRODIGY") else ())

def group_key(method, network):
    return CONTEXT, method, network, PARAMETERS[method]
