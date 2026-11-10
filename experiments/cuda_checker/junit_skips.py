"""Fail the GPU check on any skipped test (Implementation.md D8): on the box a skip hides a missing input or a CPU
build, so `run_box.sh --gpu-check` runs `python -I -m experiments.cuda_checker.junit_skips <gpu.xml>` after pytest.
Prints each skipped test with its reason; exit 1 when there is one."""
import argparse
import sys
import xml.etree.ElementTree as ET


def skipped(path: str) -> list[str]:
    return [f"{case.get('classname')}::{case.get('name')}: {s.get('message')}"
            for case in ET.parse(path).getroot().iter("testcase") for s in case.findall("skipped")]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("junit", help="the pytest --junitxml file of the GPU check")
    bad = skipped(ap.parse_args().junit)
    print(*(f"skipped on the box: {b}" for b in bad), sep="\n", file=sys.stderr)
    sys.exit(1 if bad else 0)
