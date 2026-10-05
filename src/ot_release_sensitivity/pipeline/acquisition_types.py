from __future__ import annotations

ENDPOINT = "https://open-targets-public-data-releases.s3.eu-west-1.amazonaws.com"
RELEASES = ("25.12", "26.03", "26.06")
DATASETS = (
    "association_overall_direct",
    "association_by_datatype_direct",
    "association_by_datasource_direct",
    "disease",
    "target",
)
NAMESPACE = {"s3": "http://s3.amazonaws.com/doc/2006-03-01/"}
