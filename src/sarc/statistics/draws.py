import random

from .endpoint import require

def colorectal_draws(donor_ids):
    require(
        (donor_ids == sorted(donor_ids) and len(set(donor_ids)) == len(donor_ids) == 83),
        "Full ordered 83-donor axis required",
    )
    generator = random.Random(20261004)
    return [[generator.randrange(83) for (unused) in (range(83))] for (replicate) in (range(10000))]
