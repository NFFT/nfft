#!/usr/bin/env python3
# Copyright (c) 2026 The NFFT Authors
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, see <https://www.gnu.org/licenses/>.

# /// script
# requires-python = ">=3.9"
# dependencies = []
# ///

"""Check or update the project copyright notices.

    uv run support/update_copyright.py          # report files that differ
    uv run support/update_copyright.py --fix    # rewrite them

Run from anywhere inside the repository. Only tracked files are touched.

The years in each notice are the years in which the file had a nontrivial
change, taken from git history with renames and copies followed. Files in
the initial import also get the years from FIRST_YEAR. A change is trivial if
it only touches copyright or license text, blank lines, or whitespace, so
earlier notice updates and reformats do not add years. Commits listed in
.git-blame-ignore-revs are skipped entirely. Consecutive years are written as
a range; README.md states that a range is inclusive.

Every tracked file with a project notice is updated. A notice is added to
source files of the types in STYLES that have no copyright notice at all and
at least MIN_LINES lines. Files with someone else's notice are left alone.
"""

import argparse
import datetime
import os
import re
import subprocess
import sys

HOLDER = "The NFFT Authors"
# The GPLv2 notice from the FSF's "How to Apply These Terms" appendix.
LICENSE = """\

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, see <https://www.gnu.org/licenses/>."""
BODY = LICENSE.split("\n")


def _words(lines):
    return " ".join(" ".join(lines).split())


# License texts that are replaced by BODY: the current one, and the ones used
# until 2026, which end with the FSF postal address instead of the URL.
_ADDRESS = (
    "write to the Free Software Foundation, Inc., "
    "51 Franklin Street, Fifth Floor, Boston, MA 02110-1301{} USA."
)
_URL = "see <https://www.gnu.org/licenses/>."
LICENSES = {_words(BODY)} | {
    _words(BODY).replace(_URL, _ADDRESS.format(c)) for c in ("", ",")
}

# The history starts with the CVS import of 2005; files in it date from 2002.
FIRST_YEAR = 2002

MIN_LINES = 15

# Comment style for added notices: "c" block comment, or a line prefix.
STYLES = {
    ".c": "c",
    ".h": "c",
    ".cpp": "c",
    ".dox": "c",
    ".py": "#",
    ".jl": "#",
    ".sh": "#",
    ".m4": "#",
    ".ac": "#",
    ".cmake": "#",
    "CMakeLists.txt": "#",
    ".m": "%",
}

# Files that get no added notice: generated, stubs of foreign APIs, agent
# tooling, code based on third-party code, and MATLAB Contents.m, whose first
# comment block is the toolbox description.
NO_ADD = re.compile(
    r"^(tests/data/generated/|kernel/util/kbpoly_data\.h$"
    r"|cmake/test-fixtures/|\.claude/|.*/lgwt\.m$"
    r"|.*/Contents\.m$)"
)

# A generated file carries the notice its generator writes, so both share
# their years.
GROUPS = [
    ("gen/besselgen/generate.py", "kernel/util/bessel_i0_data.h"),
    ("gen/sincgen/generate.py", "kernel/util/sinc_data.h"),
]

# Notices for the whole package take the years of all files.
PACKAGE = ["README.md"]

NOTICE = re.compile(
    r"^(?P<lead>\s*(?:\*|#|%|//)?\s*)Copyright \(c\) "
    r"(?P<years>[0-9][0-9, -]*?) "
    r"(?:Jens Keiner, Stefan Kunis, Daniel Potts|" + re.escape(HOLDER) + r")\s*$"
)
FOREIGN = re.compile(
    r"(?i)copyright\s*(\(c\)|\xa9|[0-9])|@license"
    r"|SPDX-License-Identifier"
)

# Words of all license and notice variants in the history.
VOCAB = set(
    re.findall(
        r"[a-z]+",
        " ".join(BODY).lower()
        + """
    copyright c jens keiner stefan kunis daniel potts the nfft authors and
    all rights reserved temple place suite mass ave cambridge see https http
    www gnu org licenses write to free software foundation inc franklin street
    fifth floor boston ma usa""".lower(),
    )
)
SVN_KEYWORD = re.compile(
    r"\$(Id|Rev|Revision|Date|Author|Header|HeadURL|URL"
    r"|LastChanged\w*)\b"
)


def git(*args):
    return subprocess.run(
        ["git", "-c", "core.quotepath=off", *args],
        check=True,
        capture_output=True,
        encoding="utf-8",
        errors="replace",
    ).stdout


def is_notice_line(s):
    """True if s is blank or reads as notice or license text."""
    s = s.strip().strip("/*%#!").strip()
    if not s or SVN_KEYWORD.search(s):
        return True
    words = re.findall(r"[a-z]+", s.lower())
    return (
        bool(words)
        and set(words) <= VOCAB
        and re.fullmatch(r"[A-Za-z0-9 ,.;:()<>@/'_\-]*", s) is not None
    )


def is_trivial(removed, added):
    """True if a diff only touches notice text or whitespace."""

    def key(lines):
        return sorted("".join(s.split()) for s in lines if not is_notice_line(s))

    return key(removed) == key(added)


def parse_diff(text):
    """Yield (commit, date, old_path, new_path, copy, trivial) per file diff.

    Commits start with a \\x01 line as written by --format=%x01%H %cs.
    old_path is None for an added file, new_path None for a deleted one.
    """
    commit = date = None
    cur = None

    def done():
        if cur is not None:
            return (
                commit,
                date,
                cur["old"],
                cur["new"],
                cur["copy"],
                is_trivial(cur["rem"], cur["add"]),
            )

    for line in text.split("\n"):
        if line.startswith("\x01"):
            if r := done():
                yield r
            cur = None
            commit, date = line[1:].split(" ")
        elif line.startswith("diff --git "):
            if r := done():
                yield r
            a, b = line[len("diff --git a/") :].split(" b/", 1)
            cur = {
                "old": a,
                "new": b,
                "rem": [],
                "add": [],
                "hunk": False,
                "copy": False,
            }
        elif cur is None:
            continue
        elif line.startswith("@@"):
            cur["hunk"] = True
        elif not cur["hunk"]:
            if line.startswith("new file mode"):
                cur["old"] = None
            elif line.startswith("deleted file mode"):
                cur["new"] = None
            elif line.startswith("rename from "):
                cur["old"] = line[len("rename from ") :]
            elif line.startswith("rename to "):
                cur["new"] = line[len("rename to ") :]
            elif line.startswith("copy from "):
                cur["old"] = line[len("copy from ") :]
                cur["copy"] = True
            elif line.startswith("copy to "):
                cur["new"] = line[len("copy to ") :]
        elif line.startswith("-"):
            cur["rem"].append(line[1:])
        elif line.startswith("+"):
            cur["add"].append(line[1:])
    if r := done():
        yield r


def file_years(tracked):
    """Map each tracked path to the set of years with nontrivial changes."""
    try:
        with open(".git-blame-ignore-revs") as f:
            ignore = {l.split()[0] for l in f if l.strip() and not l.startswith("#")}
    except FileNotFoundError:
        ignore = set()
    root = set(git("rev-list", "--max-parents=0", "HEAD").split())
    years = {p: set() for p in tracked}
    # Path at the current point in history -> tracked paths whose content
    # descends from it. Renames and copies extend the set; a path that later
    # belonged to a different file maps to the empty set.
    alias = {}
    log = git(
        "log",
        "--no-merges",
        "--no-color",
        "--no-ext-diff",
        "-M50%",
        "-C75%",
        "--find-copies-harder",
        "-U0",
        "-p",
        "--format=%x01%H %cs",
        "HEAD",
    )
    for commit, date, old, new, copy, trivial in parse_diff(log):
        if new is None:
            alias[old] = set()
            continue
        paths = alias.get(new, {new} & years.keys())
        if commit not in ignore:
            for path in paths:
                if not trivial:
                    years[path].add(int(date[:4]))
                if commit in root:
                    years[path].update(range(FIRST_YEAR, int(date[:4])))
        if old != new:
            if old is not None:
                if copy:
                    paths = paths | alias.get(old, {old} & years.keys())
                alias[old] = paths
            alias[new] = set()
    this_year = datetime.date.today().year
    for _, _, old, new, _, trivial in parse_diff(
        "\x01- -\n" + git("diff", "--no-color", "--no-ext-diff", "-M", "-U0", "HEAD")
    ):
        if new in years and not trivial:
            years[new].add(this_year)
    for group in GROUPS:
        union = set().union(*(years.get(p, set()) for p in group))
        for p in group:
            if p in years:
                years[p] = set(union)
    package = set().union(*years.values())
    for p in PACKAGE:
        years[p] = set(package)
    return years


def format_years(years):
    """Join a set of years into '2002-2005, 2007, 2009-2011'."""
    out = []
    for y in sorted(years):
        if out and out[-1][1] == y - 1:
            out[-1][1] = y
        else:
            out.append([y, y])
    return ", ".join(str(a) if a == b else f"{a}-{b}" for a, b in out)


def render(lead, years):
    blank = lead.rstrip()
    return [lead + "Copyright (c) " + years + " " + HOLDER] + [
        lead + l if l else blank for l in BODY
    ]


def update(lines, years):
    """Rewrite every project notice in lines. Return (lines, problems)."""
    out, problems, i = [], [], 0
    while i < len(lines):
        m = NOTICE.match(lines[i].rstrip("\r"))
        if not m:
            out.append(lines[i])
            i += 1
            continue
        lead = m.group("lead")
        blank = lead.rstrip()
        j = i + 1
        while (
            j < len(lines)
            and blank
            and lines[j].startswith(blank)
            and not lines[j].lstrip().startswith("*/")
            and is_notice_line(lines[j])
        ):
            j += 1
        while j > i + 1 and not lines[j - 1].strip().strip("*%#/"):
            j -= 1
        eol = "\r" if lines[i].endswith("\r") else ""
        if j == i + 1:
            out.append(lead + "Copyright (c) " + years + " " + HOLDER + eol)
        else:
            body = [l.rstrip("\r")[len(blank) :].strip() for l in lines[i + 1 : j]]
            if _words(body) not in LICENSES:
                problems.append(f"line {i + 1}: license text differs")
                out.extend(lines[i:j])
            else:
                out.extend(l + eol for l in render(lead, years))
        i = j
    return out, problems


def add(path, lines, years):
    """Insert a notice into a file of a type in STYLES."""
    name = os.path.basename(path)
    style = STYLES.get(name) or STYLES.get(
        os.path.splitext(path.removesuffix(".in"))[1]
    )
    if style == "c":
        return ["/*"] + render(" * ", years) + [" */", ""] + lines
    if style == "#":
        i = 0
        while (
            i < len(lines)
            and i < 2
            and (lines[i].startswith("#!") or "coding" in lines[i])
        ):
            i += 1
        return lines[:i] + render("# ", years) + [""] + lines[i:]
    # MATLAB: after the function line and its help text, which must stay the
    # first comment block.
    i = 1 if lines and re.match(r"\s*(function|classdef)\b", lines[0]) else 0
    j = i
    while j < len(lines) and lines[j].lstrip().startswith("%"):
        j += 1
    if j == i:
        return render("% ", years) + [""] + lines
    return lines[:j] + [""] + render("% ", years) + lines[j:]


def wants_notice(path, text):
    name = os.path.basename(path)
    base = path.removesuffix(".in")
    return (
        (name in STYLES or os.path.splitext(base)[1] in STYLES)
        and not NO_ADD.match(path)
        and not FOREIGN.search(text)
        and "::Package::" not in text
        and text.count("\n") >= MIN_LINES
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fix", action="store_true", help="rewrite files")
    args = ap.parse_args()
    os.chdir(git("rev-parse", "--show-toplevel").strip())
    tracked = [
        p
        for p in git("ls-files", "-z").split("\0")
        if p and os.path.isfile(p) and not os.path.islink(p)
    ]
    years = file_years(tracked)
    changed = problems = 0
    for path in tracked:
        with open(path, "rb") as f:
            raw = f.read()
        if b"\0" in raw[:8192]:
            continue
        text = raw.decode("utf-8", "surrogateescape")
        lines = text.split("\n")
        y = format_years(years[path])
        new, issues = update(lines, y)
        if new == lines and wants_notice(path, text):
            new = add(path, lines, y)
        for issue in issues:
            print(f"{path}: {issue}", file=sys.stderr)
            problems += 1
        if new != lines:
            changed += 1
            print(path)
            if args.fix:
                with open(path, "wb") as f:
                    f.write("\n".join(new).encode("utf-8", "surrogateescape"))
    verb = "updated" if args.fix else "differ"
    print(f"{changed} file(s) {verb}, {problems} problem(s)", file=sys.stderr)
    return 1 if problems or (changed and not args.fix) else 0


if __name__ == "__main__":
    sys.exit(main())
