import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
# The lab oracle (phase1_lab/common/oracle.py) is THE spec; never copied, always imported.
ORACLE_DIR = Path(os.environ.get("POLYCORE_ORACLE_DIR", HERE.parents[1] / "common"))
sys.path.insert(0, str(ORACLE_DIR))
sys.path.insert(0, str(HERE))
sys.dont_write_bytecode = True


def env_int(name, default):
    return int(os.environ.get(name, default))
