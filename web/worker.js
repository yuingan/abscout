// AbTaxa web version: runs the unchanged AbTaxa Python source (abtaxa.py,
// resolve_orthologs.py, structure_compare.py) inside the browser with Pyodide.
// Nothing is sent to any AbTaxa server; queries go straight from this browser to
// Ensembl and UniProt, exactly as the command-line version makes them.
//
// build_site.py fills in PYODIDE_URL and copies the Python files into py/.

const PYODIDE_URL = "__PYODIDE_URL__";
const PY_FILES = ["structure_compare.py", "resolve_orthologs.py", "abtaxa.py"];
const PACKAGES = ["numpy", "biopython", "requests"];

let pyodide = null;
let ready = null;          // the boot promise, shared by every message
let booted = false;

function post(type, extra) {
  self.postMessage(Object.assign({ type }, extra || {}));
}

async function boot() {
  post("status", { text: "Loading the Python runtime (first visit only: about 15 MB)..." });
  const { loadPyodide } = await import(PYODIDE_URL + "pyodide.mjs");
  pyodide = await loadPyodide({ indexURL: PYODIDE_URL });
  // AbTaxa writes progress to stderr; show it line by line while a run is going
  pyodide.setStderr({ batched: (line) => post("log", { line }) });
  pyodide.setStdout({ batched: () => {} });

  post("status", { text: "Loading numpy and Biopython..." });
  await loadPackages();

  post("status", { text: "Loading AbTaxa..." });
  pyodide.FS.mkdirTree("/home/pyodide/abtaxa");
  for (const f of PY_FILES) {
    const r = await fetch("py/" + f, { cache: "no-cache" });
    if (!r.ok) throw new Error("could not load " + f + " (" + r.status + ")");
    pyodide.FS.writeFile("/home/pyodide/abtaxa/" + f, await r.text());
  }
  const version = pyodide.runPython(`
import sys
sys.path.insert(0, "/home/pyodide/abtaxa")
import abtaxa
abtaxa.__version__`);
  booted = true;
  post("ready", { version });
}

// Pyodide's loadPackage does not throw when a package fails to download or install: it only
// logs to the console, and the failure surfaced later as "No module named 'Bio'". Collect its
// messages, retry once (CDN hiccups are usually transient), and stop with the reason if a
// package is still missing.
async function loadPackages() {
  const problems = [];
  const opts = { messageCallback: () => {}, errorCallback: (m) => problems.push(m) };
  const missing = () => PACKAGES.filter((p) => !pyodide.loadedPackages[p]);
  for (let attempt = 1; attempt <= 2 && missing().length; attempt++) {
    if (attempt > 1) post("status", { text: "Retrying " + missing().join(", ") + "..." });
    try {
      await pyodide.loadPackage(missing(), opts);
    } catch (e) {
      problems.push(String(e && e.message ? e.message : e));
    }
  }
  if (missing().length) {
    throw new Error(
      "Your browser could not load " + missing().join(", ") + ", part of the Python runtime " +
      "the page downloads from " + new URL(PYODIDE_URL, self.location).host + ". Reload the " +
      "page to try again. If it keeps failing, a network filter or browser extension may be " +
      "blocking that download; please report it at https://github.com/yuingan/abtaxa/issues " +
      "with the details below." + (problems.length ? "\n\nDetails:\n" + problems.join("\n") : ""));
  }
}

// One cheap request to each service before the first run. If the browser cannot reach
// them (offline, firewall, an extension blocking requests), say so plainly instead of
// letting it surface later as "no ortholog found".
let reachable = false;
async function checkServices() {
  if (reachable) return;
  const probes = [
    ["Ensembl", "https://rest.ensembl.org/info/ping?content-type=application/json"],
    ["UniProt", "https://rest.uniprot.org/uniprotkb/P07766.fasta"],
  ];
  for (const [name, url] of probes) {
    try {
      const r = await fetch(url);
      if (!r.ok) throw new Error("HTTP " + r.status);
    } catch (e) {
      throw new Error(
        "This browser could not reach " + name + " (" + e.message + "). AbTaxa needs " +
        "Ensembl and UniProt; check the network connection, or whether a firewall or " +
        "browser extension blocks rest.ensembl.org / rest.uniprot.org, then try again.");
    }
  }
  reachable = true;
}

const RUN_PY = `
import json, abtaxa
a = dict(ARGS)
try:
    ep = abtaxa.parse_epitope_coor(a["epitope_coor"]) if a.get("epitope_coor") else None
    out = abtaxa.run_single(a["gene"], int(a["taxon"]), ep,
                             a.get("epitope_ref") or "human", a.get("epitope_seq") or None)
    res = {"ok": True, "html": out["html"], "result": json.dumps(out["result"], indent=2),
           "report": out["report"]}
except (ValueError, abtaxa.AbTaxaError) as exc:
    res = {"ok": False, "error": str(exc)}
json.dumps(res)
`;

self.onmessage = async (ev) => {
  const msg = ev.data || {};
  try {
    if (!ready) ready = boot();
    await ready;
    if (msg.type !== "run") return;
    post("status", { text: "Checking that Ensembl and UniProt can be reached..." });
    await checkServices();
    post("status", { text: "Running AbTaxa..." });
    pyodide.globals.set("ARGS", pyodide.toPy(msg.args));
    const out = JSON.parse(pyodide.runPython(RUN_PY));
    post("done", out);
  } catch (e) {
    if (!booted) ready = null;                   // a failed load can be retried
    post("done", { ok: false, error: String(e && e.message ? e.message : e),
                   unexpected: booted && !(e && e.message && e.message.startsWith("This browser")) });
  }
};
