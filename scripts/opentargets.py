from pathlib import Path
import argparse
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ot_release_sensitivity.workflow import COMMANDS

parser = argparse.ArgumentParser(
    description = "Compare OpenTargets release scores, supports and priorities"
)
parser.add_argument("command", choices = tuple(COMMANDS))
parser.add_argument("--root", type = Path, required = True)
parser.add_argument("--config-dir", type = Path)
parser.add_argument("--threads", type = int, default = 8)
parser.add_argument("--memory-limit", default = "8GB")
parser.add_argument("--workers", type = int, default = 8)
arguments = parser.parse_args()
COMMANDS[arguments.command](
    arguments.root,
    arguments.config_dir,
    arguments.threads,
    arguments.memory_limit,
    arguments.workers,
)
