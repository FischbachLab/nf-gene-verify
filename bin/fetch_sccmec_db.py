#!/usr/bin/env python3
"""Deprecated shim -- use fetch_predefined_db.py --db sccmec.

Kept so existing commands keep working. The SCCmec and IS256 builders share
the same output format and manifest, so they live in one script.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from fetch_predefined_db import main  # noqa: E402

if __name__ == "__main__":
    print("NOTE: fetch_sccmec_db.py is a shim; "
          "use fetch_predefined_db.py --db sccmec", file=sys.stderr)
    sys.exit(main(["--db", "sccmec"] + sys.argv[1:]))
