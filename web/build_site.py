"""
Assemble the AbTaxa web version into a folder that any static host (GitHub Pages) can serve.

The site runs the repository's own abtaxa.py, resolve_orthologs.py and structure_compare.py,
copied byte for byte, so the web version cannot drift from the source code: every number it
shows comes from the same functions as the command line.

    python web/build_site.py                 # -> _site/
    python -m http.server -d _site 8000      # then open http://localhost:8000
"""
import argparse
import json
import os
import re
import shutil
from pathlib import Path

PYODIDE_VERSION = "314.0.7"
PYODIDE_CDN = f"https://cdn.jsdelivr.net/pyodide/v{PYODIDE_VERSION}/full/"
PY_FILES = ["abtaxa.py", "resolve_orthologs.py", "structure_compare.py"]
WEB_FILES = ["index.html", "worker.js"]


def build(out, pyodide_url=PYODIDE_CDN):
    web = Path(__file__).resolve().parent
    root = web.parent
    out = Path(out)
    if out.exists():
        shutil.rmtree(out)
    (out / "py").mkdir(parents=True)

    for f in PY_FILES:
        shutil.copyfile(root / f, out / "py" / f)
    for f in WEB_FILES:
        text = (web / f).read_text()
        if f == "worker.js":
            if not pyodide_url.endswith("/"):
                pyodide_url += "/"
            text = text.replace("__PYODIDE_URL__", pyodide_url)
        (out / f).write_text(text)
    (out / ".nojekyll").write_text("")

    version = re.search(r'^__version__ = "([^"]+)"', (root / "abtaxa.py").read_text(), re.M)
    info = {"abtaxa": version.group(1) if version else "?",
            "pyodide": PYODIDE_VERSION,
            "commit": os.environ.get("GITHUB_SHA", "")[:12]}
    (out / "version.json").write_text(json.dumps(info, indent=2) + "\n")
    return info


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Build the AbTaxa web version into a folder.")
    ap.add_argument("--out", default="_site", help="output folder (default: _site)")
    ap.add_argument("--pyodide-url", default=PYODIDE_CDN,
                    help="where Pyodide is served from (default: the jsDelivr CDN)")
    a = ap.parse_args()
    info = build(a.out, a.pyodide_url)
    print(f"built {a.out}/ with AbTaxa {info['abtaxa']} on Pyodide {info['pyodide']}")
