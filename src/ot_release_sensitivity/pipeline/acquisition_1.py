from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.parse import quote
from xml.etree import ElementTree
import os
import time
import requests
from .acquisition_types import DATASETS, ENDPOINT, NAMESPACE, RELEASES

def request(method, url, **kwargs):
    last_error = None
    for (attempt) in (range(7)):
        try:
            response = requests.request(method, url, timeout = (30, 300), **kwargs)
            response.raise_for_status()
            return response
        except requests.RequestException as error:
            last_error = error
            if (attempt < 6):
                time.sleep(min(2**attempt, 30))
    raise last_error

def list_response(prefix, delimiter = None, continuation_token = None):
    params = {"list-type": "2", "prefix": prefix, "max-keys": "1000"}
    if (delimiter):
        params["delimiter"] = delimiter
    if (continuation_token):
        params["continuation-token"] = continuation_token
    return request("GET", ENDPOINT, params = params)

def parse_listing(content):
    root = ElementTree.fromstring(content)
    objects = []
    prefixes = []
    for (item) in (root.findall("s3:Contents", NAMESPACE)):
        objects.append(
            {
                "key": item.findtext("s3:Key", namespaces = NAMESPACE),
                "last_modified": item.findtext("s3:LastModified", namespaces = NAMESPACE),
                "etag": item.findtext("s3:ETag", namespaces = NAMESPACE).strip('"'),
                "size": int(item.findtext("s3:Size", namespaces = NAMESPACE)),
                "storage_class": item.findtext(
                    "s3:StorageClass", default = "", namespaces = NAMESPACE
                ),
            }
        )
    for (item) in (root.findall("s3:CommonPrefixes", NAMESPACE)):
        prefixes.append(item.findtext("s3:Prefix", namespaces = NAMESPACE))
    truncated = (
        root.findtext("s3:IsTruncated", default = "false", namespaces = NAMESPACE).lower()
        == "true"
    )
    token = root.findtext("s3:NextContinuationToken", namespaces = NAMESPACE)
    return (objects, prefixes, truncated, token)

def list_objects(prefix):
    objects = []
    token = None
    while (True):
        response = list_response(prefix, continuation_token = token)
        page_objects, _, truncated, token = parse_listing(response.content)
        objects.extend(page_objects)
        if (not truncated):
            return objects

def download_object(item, release, dataset, root):
    prefix = f"platform/{release}/output/{dataset}/"
    relative = item["key"][len(prefix) :]
    destination = Path(root) / "data" / "raw" / release / dataset / relative
    if (not destination.resolve().is_relative_to(Path(root).resolve())):
        raise ValueError("Source object path escapes the working directory")
    destination.parent.mkdir(parents = True, exist_ok = True)
    if (destination.exists() and destination.stat().st_size == item["size"]):
        return destination
    partial = destination.with_suffix(destination.suffix + ".partial")
    existing = partial.stat().st_size if (partial.exists()) else 0
    headers = {"Range": f"bytes={existing}-"} if (0 < existing < item["size"]) else {}
    mode = "ab" if (headers) else "wb"
    url = f"{ENDPOINT}/{quote(item['key'], safe='/')}"
    with request("GET", url, headers = headers, stream = True) as response:
        if (headers and response.status_code != 206):
            mode = "wb"
        with partial.open(mode) as handle:
            for (block) in (response.iter_content(chunk_size = 8 * 1024 * 1024)):
                if (block):
                    handle.write(block)
    if (partial.stat().st_size != item["size"]):
        raise RuntimeError(f"Incomplete source object: {item['key']}")
    os.replace(partial, destination)
    return destination

def run(root, workers = 8):
    objects = []
    for (release) in (RELEASES):
        for (dataset) in (DATASETS):
            prefix = f"platform/{release}/output/{dataset}/"
            entries = list_objects(prefix)
            if (not any((item["key"].endswith(".parquet") for (item) in (entries)))):
                raise FileNotFoundError(f"No source partitions for {release}/{dataset}")
            objects.extend(((item, release, dataset) for (item) in (entries)))
    with ThreadPoolExecutor(max_workers = workers) as executor:
        pending = [
            executor.submit(download_object, item, release, dataset, root)
            for ((item, release, dataset)) in (objects)
        ]
        return [future.result() for (future) in (as_completed(pending))]
