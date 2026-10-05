import json
import subprocess
from pathlib import Path

def read_rds(requests, output_directory, rscript = "Rscript", base_directory = "."):
    base = Path(base_directory).resolve()
    output = (base / output_directory).resolve()
    if (output.exists()):
        raise FileExistsError("Use a new typed-reader output directory")
    identities = [row["source_id"] for (row) in (requests)]
    if (len(identities) != len(set(identities)) or any(
        not isinstance(value, str) or not value for (value) in (identities)
    )):
        raise ValueError("Unique explicit source identities required")
    lines, files = [], []
    for (index, row) in (enumerate(requests)):
        if (row["mode"] not in ("checkpoint", "result", "pathways")):
            raise ValueError("Unsupported scientific RDS projection")
        source = (base / row["path"]).resolve()
        if (not source.is_file()):
            raise FileNotFoundError(source)
        destination = output / (str(index).zfill(5) + ".json")
        values = (row["mode"], str(source), str(destination))
        if (any("\t" in value or "\n" in value or "\r" in value for (value) in (values))):
            raise ValueError("Unsupported control character in RDS request")
        lines.append("\t".join(values))
        files.append(destination)
    output.mkdir(parents = True)
    request_path = output / "requests.tsv"
    request_path.write_text("\n".join(lines) + "\n", encoding = "utf-8")
    helper = Path(__file__).parent / "r" / "pro_hypothesis_reader.R"
    executable = (
        str((base / rscript).resolve())
        if ("/" in rscript or "\\" in rscript)
        else rscript
    )
    result = subprocess.run(
        [executable, str(helper), str(request_path)],
        capture_output = True,
        text = True,
        check = False,
    )
    (output / "reader_stdout.txt").write_text(result.stdout, encoding = "utf-8")
    (output / "reader_stderr.txt").write_text(result.stderr, encoding = "utf-8")
    if (result.returncode):
        raise RuntimeError(
            "RDS reader failed; partial outputs and diagnostics are retained"
        )
    return {
        "payloads": {
            identity: json.loads(path.read_text(encoding = "utf-8"))
            for ((identity, path)) in (zip(identities, files))
        },
        "projection_modes": {row["source_id"]: row["mode"] for (row) in (requests)},
    }
