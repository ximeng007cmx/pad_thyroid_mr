#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
relocate.py — point the analysis scripts at *your* local checkout.

Background
----------
The scripts in this release were written for a fixed working directory and each
contains exactly one hard-coded project root, e.g.

    ROOT = r"D:/WorkSpace_Study/pad_thyroid_mr"

The scripts are shipped unmodified so that the archived code is byte-identical
to the code that produced the published results. Running this script once
rewrites that root to the location where you cloned/extracted this repository.

Usage
-----
    python relocate.py            # target = this file's directory (recommended)
    python relocate.py /path/to/clone
    python relocate.py --dry-run  # show what would change, write nothing
    python relocate.py --restore  # revert to the original root

The script only rewrites the project-root literal. Nothing else is touched.
"""

import argparse
import re
import sys
from pathlib import Path

# The literal that appears in every script, in either slash style.
PATTERN = re.compile(r"[A-Za-z]:[\\/]WorkSpace_Study[\\/]pad_thyroid_mr")
ORIGINAL = "D:/WorkSpace_Study/pad_thyroid_mr"


def iter_scripts(root: Path):
    for p in sorted(root.rglob("*.py")):
        if p.name == Path(__file__).name:
            continue
        yield p


def main() -> int:
    ap = argparse.ArgumentParser(description="Rewrite the hard-coded project root.")
    ap.add_argument("target", nargs="?", default=None,
                    help="new project root (default: the directory containing this script)")
    ap.add_argument("--dry-run", action="store_true", help="show changes without writing")
    ap.add_argument("--restore", action="store_true",
                    help="revert to the original root (D:/WorkSpace_Study/pad_thyroid_mr)")
    args = ap.parse_args()

    here = Path(__file__).resolve().parent

    if args.restore:
        new_root = ORIGINAL
    else:
        new_root = str(Path(args.target).resolve()) if args.target else str(here)
        if sys.platform == "win32":
            new_root = new_root.replace("\\", "/")

    if not args.restore and not Path(new_root).is_dir():
        print(f"ERROR: target directory does not exist: {new_root}", file=sys.stderr)
        return 2

    print(f"project root  ->  {new_root}\n")

    changed_files = 0
    changed_lines = 0
    for script in iter_scripts(here):
        text = script.read_text(encoding="utf-8", errors="surrogateescape")
        # lambda avoids backslash-escaping pitfalls inside re.sub replacements
        new_text, n = PATTERN.subn(lambda _m: new_root, text)
        rel = script.relative_to(here)
        if n:
            changed_files += 1
            changed_lines += n
            print(f"  [{n}] {rel}")
            if not args.dry_run:
                script.write_text(new_text, encoding="utf-8", errors="surrogateescape")

    verb = "would rewrite" if args.dry_run else "rewrote"
    print(f"\n{verb} {changed_lines} line(s) across {changed_files} file(s).")
    if changed_files == 0:
        print("Nothing matched — the root may already point elsewhere, "
              "or the pattern was edited.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
