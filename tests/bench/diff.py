"""Pure diff of two sets of per-testbed BMF files on the accuracy-digits measure."""
from __future__ import annotations

import glob
import json
import os
from dataclasses import dataclass, field


@dataclass
class Change:
    testbed: str
    name: str
    base_digits: float
    pr_digits: float
    delta_digits: float
    pct: float


@dataclass
class DiffResult:
    improvements: list
    regressions: list
    unchanged_count: int
    added: list = field(default_factory=list)
    removed: list = field(default_factory=list)
    by_testbed: dict = field(default_factory=dict)


def digits(bmf, name):
    return float(bmf[name]["accuracy-digits"]["value"])


def load_bmf_tree(directory):
    tree = {}
    for path in sorted(glob.glob(os.path.join(directory, "*.bmf.json"))):
        testbed = os.path.basename(path)[: -len(".bmf.json")]
        with open(path, encoding="utf-8") as f:
            tree[testbed] = json.load(f)
    return tree


def diff(pr_tree, base_tree, gate=0.5):
    improvements, regressions, added, removed = [], [], [], []
    unchanged = 0
    by_testbed = {}
    for testbed in sorted(set(pr_tree) | set(base_tree)):
        pr_bmf = pr_tree.get(testbed, {})
        base_bmf = base_tree.get(testbed, {})
        counts = {"improved": 0, "regressed": 0, "unchanged": 0}
        for name in sorted(set(pr_bmf) | set(base_bmf)):
            in_pr, in_base = name in pr_bmf, name in base_bmf
            if in_pr and not in_base:
                added.append(f"{testbed}::{name}")
                continue
            if in_base and not in_pr:
                removed.append(f"{testbed}::{name}")
                continue
            b, p = digits(base_bmf, name), digits(pr_bmf, name)
            delta = p - b
            pct = (100.0 * delta / b) if b != 0 else 0.0
            if abs(delta) >= gate:
                change = Change(testbed, name, b, p, delta, pct)
                (improvements if delta > 0 else regressions).append(change)
                counts["improved" if delta > 0 else "regressed"] += 1
            else:
                unchanged += 1
                counts["unchanged"] += 1
        by_testbed[testbed] = counts
    improvements.sort(key=lambda c: abs(c.delta_digits), reverse=True)
    regressions.sort(key=lambda c: abs(c.delta_digits), reverse=True)
    return DiffResult(improvements, regressions, unchanged, added, removed, by_testbed)
