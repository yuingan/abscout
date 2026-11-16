#!/usr/bin/env python3
"""
abscout.py — AbScout: cross-species antibody cross-reactivity by membrane region.

One command. Given a human gene symbol and a target species, AbScout:
  1. resolves orthologs in the antibody-rich model species + the target (Ensembl Compara)
  2. determines each MODEL protein's membrane topology (curated UniProt where annotated,
     otherwise TMbed with --out-format 4)
  3. compares the target against each model per region (structure_compare.py)
  4. reports per-region sequence identity, and -- when an antibody's immunogen is supplied --
     that epitope's conservation, coordinates, indels and alignment across species.

The report is deliberately neutral: it gives percentages and coordinates for the user to
interpret, and makes no claim about whether a given antibody will bind.

Chains: resolve_orthologs.py  +  TMbed  +  structure_compare.py  (all must be importable/installed).

USAGE
  python abscout.py --gene EGFR --target-taxon 7955 --json egfr_abscout.json
  python abscout.py --batch cases.csv --batch-out results.csv
  # taxids: 7955 zebrafish | 10181 naked mole-rat | 9541 cynomolgus | 9031 chicken | 9823 pig

SCOPE: identity is a proxy for epitope conservation. It cannot see conformational epitopes,
differential glycosylation, or single hot-spot contact residues, so treat the output as a
prioritisation screen rather than a binding prediction.
"""
import argparse
import csv
import datetime
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import requests

import resolve_orthologs as ro
import structure_compare as sc

MODEL_ORDER = ["human", "mouse", "rat", "rabbit"]
__version__ = "0.8.0"

# ----------------------------------------------------------------------------------
# On-disk cache. Sequences and topology for a given accession never change, so caching
# them makes repeat and batch runs near-instant and avoids hammering the public APIs.
# ----------------------------------------------------------------------------------
CACHE_DIR = Path(os.environ.get("ABSCOUT_CACHE", Path.home() / ".abscout_cache"))
CACHE_ENABLED = True
# What AbScout *derives* from a UniProt entry (topology, anchor type, N-terminal processing,
# the reviewed-entry search) depends on its parsing code, so it is cached per version: a
# release that parses better must not be served an answer an older release worked out.
# Human THY1 kept a pre-GPI "no topology" answer this way and was sent to TMbed for good.
# Sequences and TMbed predictions are raw data and are kept across versions.
DERIVED_KINDS = {"topo_uniprot", "topokind_uniprot", "nterm_uniprot", "reviewed_alts"}


def _cache_path(kind, key):
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in str(key))
    if kind in DERIVED_KINDS:
        kind = f"{kind}@{__version__}"
    return CACHE_DIR / kind / f"{safe}.txt"


def cache_get(kind, key):
    if not CACHE_ENABLED:
        return None
    p = _cache_path(kind, key)
    try:
        return p.read_text() if p.exists() else None
    except OSError:
        return None


def cache_put(kind, key, value):
    if not CACHE_ENABLED or value is None:
        return
    p = _cache_path(kind, key)
    try:
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(value)
    except OSError:
        pass                                    # a cache failure must never break a run


MIN_TMBED_LEN = 20          # TMbed cannot batch very short sequences (crashes in make_batches)
IN_BROWSER = sys.platform == "emscripten"       # the web version runs this file under Pyodide


def tmbed_available():
    """TMbed is an optional external tool (installed with conda). The web version never has
    it, and a pip-only install may not either; topology then comes from UniProt alone."""
    import importlib.util
    try:
        return importlib.util.find_spec("tmbed") is not None
    except (ImportError, ValueError):
        return False


def _tmbed_hint():
    return ("the source-code version of AbScout predicts it with TMbed" if IN_BROWSER
            else "install TMbed (see README) to predict it")


def run_tmbed(seq, work, tag, acc=None):
    """Run TMbed on a sequence and return its per-residue topology string (o/i/H/S...).
    Results are cached by accession: TMbed is the slowest step (~20-40 s per protein)."""
    seq = "".join(seq.split()).upper()
    if len(seq) < MIN_TMBED_LEN:
        raise ValueError(f"sequence too short for TMbed ({len(seq)} aa; likely a fragment)")
    if acc:
        cached = cache_get("topo_tmbed", acc)
        if cached:
            return cached
    fasta = Path(work) / f"{tag}.fasta"
    pred = Path(work) / f"{tag}.pred"
    fasta.write_text(f">{tag}\n{seq}\n")
    proc = subprocess.run([sys.executable, "-m", "tmbed", "predict",
                           "-f", str(fasta), "-p", str(pred), "--out-format", "4"],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        tail = (proc.stderr or "").strip().splitlines()
        raise RuntimeError("TMbed failed: " + (tail[-1] if tail else "unknown error"))
    preds = sc.parse_tmbed_pred(str(pred))
    topo = next(iter(preds.values())) if preds else None
    if acc and topo:
        cache_put("topo_tmbed", acc, topo)
    return topo


def _uniprot_get(url, tries=4):
    """GET a UniProt URL with polite backoff on rate limits / transient errors."""
    r = None
    for i in range(1, tries + 1):
        try:
            r = requests.get(url, timeout=30)
        except requests.RequestException:
            time.sleep(1.5 * i)
            continue
        if r.status_code == 200:
            return r
        if r.status_code in (429, 500, 502, 503, 504) and i < tries:
            time.sleep(min(float(r.headers.get("Retry-After", 1.5 * i)), 20))
            continue
        break
    if r is not None:
        r.raise_for_status()
    raise requests.RequestException(f"failed to GET {url}")


def fetch_uniprot_seq(acc):
    """Fetch the canonical protein sequence for a UniProt accession (no structure needed)."""
    cached = cache_get("seq", acc)
    if cached:
        return cached
    r = _uniprot_get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta")
    seq = "".join(ln.strip() for ln in r.text.splitlines() if ln and not ln.startswith(">"))
    if not seq:
        # UniProt answers 200 with an empty FASTA for an entry it has withdrawn, and Ensembl
        # xrefs can still point at one (naked mole-rat THY1 -> A0A0P6J0R8, deleted as "not
        # part of a reference proteome"). Carrying on with a 0 aa protein made every species
        # fail with "zero length"; say what actually happened instead.
        raise UniProtUnavailable(_inactive_reason(acc) or "UniProt returned no sequence")
    cache_put("seq", acc, seq)
    return seq


class UniProtUnavailable(Exception):
    """A UniProt accession that exists as a reference but has no usable sequence."""


def _inactive_reason(acc):
    """Why UniProt holds no sequence for an accession ("deleted from UniProtKB (...)"), or ""."""
    try:
        entry = _uniprot_get(f"https://rest.uniprot.org/uniprotkb/{acc}.json").json()
    except Exception:
        return ""
    if entry.get("entryType") != "Inactive":
        return ""
    why = entry.get("inactiveReason") or {}
    kind = (why.get("inactiveReasonType") or "withdrawn").lower()
    detail = why.get("deletedReason") or ""
    return f"{kind} from UniProtKB" + (f" ({detail})" if detail else "")


def _topology_string_from_features(features, seq_len):
    """Build a per-residue topology string (o/i/H/S/.) from UniProt features.
    Returns None unless the entry annotates a transmembrane segment AND BOTH inside and
    outside domains. An incomplete annotation (e.g. a TM with no sides, common in TrEMBL)
    is worse than a full TMbed prediction, so in that case we fall back to TMbed."""
    labels = ["."] * seq_len
    has_tm = has_extra = has_intra = False
    for f in features:
        loc = f.get("location", {})
        s = loc.get("start", {}).get("value")
        e = loc.get("end", {}).get("value")
        if not s or not e:
            continue
        ftype = f.get("type", "")
        if ftype in ("Transmembrane", "Intramembrane"):
            has_tm, code = True, "H"
        elif ftype == "Signal":
            code = "S"
        elif ftype == "Topological domain":
            desc = (f.get("description") or "").lower()
            if "extracellular" in desc or "lumen" in desc:
                code, has_extra = "o", True
            elif "cytoplasm" in desc:
                code, has_intra = "i", True
            else:
                continue
        else:
            continue
        for i in range(s - 1, min(e, seq_len)):
            labels[i] = code
    if has_tm and has_extra and has_intra:
        return "".join(labels)
    return None


def _gpi_topology_from_features(features, seq_len):
    """
    Topology for a GPI-anchored protein, or None if the entry has no GPI anchor.

    A GPI-anchored protein is tethered to the outer leaflet by a lipid rather than by a
    transmembrane helix, so TMbed correctly finds no TM and the protein would otherwise be
    reported as soluble -- badly misleading for CD90, CD59, CD73, PrP and other common
    surface targets, which are exactly the proteins people check before buying a
    flow-cytometry antibody. The mature chain of such a protein lies wholly outside the
    cell, so every residue of it is extracellular.

    UniProt marks these with a Lipidation feature whose description names a GPI anchor
    (THY1/P04216: Signal 1-19, Chain 20-130, Propeptide 131-161 "Removed in mature form",
    Lipidation 130 "GPI-anchor amidated cysteine"). The C-terminal propeptide is the GPI
    attachment signal and is cleaved off, so it is left unassigned rather than called
    extracellular.
    """
    def span(f):
        loc = f.get("location", {})
        return loc.get("start", {}).get("value"), loc.get("end", {}).get("value")

    has_gpi = any(f.get("type") == "Lipidation"
                  and "gpi-anchor" in (f.get("description") or "").lower()
                  for f in features)
    if not has_gpi:
        return None
    # a GPI anchor and a transmembrane helix are different attachments; if the entry
    # annotates a helix as well, the normal TM path describes it better
    if any(f.get("type") in ("Transmembrane", "Intramembrane") for f in features):
        return None

    labels = ["."] * seq_len
    chains = [span(f) for f in features if f.get("type") == "Chain"]
    chains = [(s, e) for s, e in chains if s and e]
    if not chains:
        return None
    for s, e in chains:                             # mature chain: entirely extracellular
        for i in range(s - 1, min(e, seq_len)):
            labels[i] = "o"
    for f in features:                              # signal peptide, if any
        if f.get("type") == "Signal":
            s, e = span(f)
            if s and e:
                for i in range(s - 1, min(e, seq_len)):
                    labels[i] = "S"
    for f in features:                              # cleaved GPI attachment signal
        if f.get("type") == "Propeptide":
            s, e = span(f)
            if s and e and s > min(c[0] for c in chains):
                for i in range(s - 1, min(e, seq_len)):
                    labels[i] = "."
    return "".join(labels) if "o" in labels else None


def _nterm_processing_from_features(features):
    """
    N-terminal processing from UniProt features, for the coordinate-numbering check.
    Returns {"signal_len", "mature_start", "parts"} where `mature_start` is the 1-based first
    residue of the mature chain: after the signal peptide and any propeptide that directly
    follows it (E-cadherin: signal 1-22, propeptide 23-154, mature chain from 155), or after
    a removed initiator methionine. `parts` names what precedes the mature chain.
    """
    def span(f):
        loc = f.get("location", {})
        return loc.get("start", {}).get("value"), loc.get("end", {}).get("value")

    sig = 0
    for f in features:
        if f.get("type") == "Signal":
            s, e = span(f)
            if s == 1 and e:
                sig = max(sig, e)
    parts = [f"{sig}-residue signal peptide"] if sig else []
    pos = sig
    if not sig:
        for f in features:
            if f.get("type") == "Initiator methionine" and span(f) == (1, 1):
                pos = 1
                parts.append("initiator methionine removed")
                break
    moved = True
    while moved:                                # chain through contiguous propeptides
        moved = False
        for f in features:
            if f.get("type") == "Propeptide":
                s, e = span(f)
                if s == pos + 1 and e and e > pos:
                    pos, moved = e, True
                    if "propeptide" not in parts:
                        parts.append("propeptide")
    return {"signal_len": sig, "mature_start": pos + 1, "parts": parts}


def _fetch_uniprot_annotations(acc, seq_len):
    """Fetch an entry once and cache everything derived from it: the topology string, how
    the protein is anchored, and the N-terminal processing. Returns
    (topology_or_None, kind, nterm_dict) where kind is "tm", "gpi" or "" ; raises on a
    transient network failure so nothing misleading is cached."""
    r = _uniprot_get(f"https://rest.uniprot.org/uniprotkb/{acc}.json")
    features = r.json().get("features", [])
    topo = _topology_string_from_features(features, seq_len)
    kind = "tm" if topo else ""
    if topo is None:                            # no helix annotated: is it GPI-anchored?
        topo = _gpi_topology_from_features(features, seq_len)
        kind = "gpi" if topo else ""
    nterm = _nterm_processing_from_features(features)
    cache_put("topo_uniprot", acc, topo or "")
    cache_put("topokind_uniprot", acc, kind)
    cache_put("nterm_uniprot", acc, json.dumps(nterm))
    return topo, kind, nterm


def uniprot_topology(acc, seq_len):
    """Curated topology for an accession as (topology_or_None, kind). kind is "tm" for a
    transmembrane protein, "gpi" for a GPI-anchored one, "" when nothing is annotated."""
    cached = cache_get("topo_uniprot", acc)
    if cached is not None:
        return (cached or None), (cache_get("topokind_uniprot", acc) or "")
    try:
        topo, kind, _ = _fetch_uniprot_annotations(acc, seq_len)
        return topo, kind
    except requests.RequestException:
        return None, ""                         # transient: don't cache a failure


def uniprot_topology_string(acc, seq_len):
    """Curated UniProt topology for a reference accession, or None if unannotated."""
    return uniprot_topology(acc, seq_len)[0]


def uniprot_nterm(acc, seq_len):
    """Cached N-terminal processing (signal peptide / propeptide) for an accession, or None
    if it cannot be determined right now."""
    cached = cache_get("nterm_uniprot", acc)
    if cached:
        try:
            return json.loads(cached)
        except ValueError:
            pass
    try:
        return _fetch_uniprot_annotations(acc, seq_len)[2]
    except requests.RequestException:
        return None


# ----------------------------------------------------------------------------------
# N-linked glycosylation sequons
#
# A glycan is bulky and can physically block antibody access, so a sequon GAINED in the
# target (or LOST relative to the immunogen) can change reactivity even when sequence
# identity is high. Identity cannot see this, which is why it is reported separately.
# ----------------------------------------------------------------------------------
SEQUON_RE = re.compile(r"N(?=[^P][ST])")        # N-X-[S/T], X != P; lookahead finds overlaps


def sequon_positions(seq):
    """0-based positions of the asparagine in each N-linked glycosylation sequon."""
    return {m.start() for m in SEQUON_RE.finditer(seq)}


def glyco_delta(ref_seq, tgt_seq, cols, r2t):
    """
    Compare N-glycosylation sequons between the reference and the target over the epitope
    columns `cols` (0-based reference positions), using the reference->target index map.
    Positions are reported in 1-based REFERENCE numbering so they line up with the epitope
    table and the alignment.
    """
    ref_sq, tgt_sq = sequon_positions(ref_seq), sequon_positions(tgt_seq)
    gained, lost = [], []
    for a in cols:
        b = r2t.get(a)
        in_ref = a in ref_sq
        in_tgt = b is not None and b in tgt_sq
        if in_tgt and not in_ref:
            gained.append(a + 1)
        elif in_ref and not in_tgt:
            lost.append(a + 1)
    return {"gained_in_target": gained, "lost_in_target": lost}


def signal_peptide_len(topology):
    """
    Length of the leading signal peptide implied by a topology string (count of leading
    'S'). Zero when there is none. Used to warn about coordinate numbering: UniProt counts
    from Met1, while datasheets sometimes count from the mature protein.
    """
    if not topology:
        return 0
    n = 0
    for c in topology:
        if c == "S":
            n += 1
        else:
            break
    return n


_REF_MAPS = {}


def _ref_map(ref_seq, seq):
    """Reference index -> species index from a global alignment, memoised: a protein with
    several regions asks for the same alignment once per segment."""
    key = (ref_seq, seq)
    if key not in _REF_MAPS:
        _REF_MAPS[key] = dict(sc.correspondence_seq({"seq": ref_seq}, {"seq": seq})[0])
    return _REF_MAPS[key]


REGION_LABELS = {"extracellular": "Extracellular", "TM": "Transmembrane",
                 "intracellular": "Intracellular"}


def topology_segments(topology):
    """Contiguous runs of each region in a topology string, as
    {region: [(start, end), ...]} in 1-based numbering. Signal peptide and unassigned
    residues are not regions and are left out."""
    segs = {k: [] for k in REGION_LABELS}
    if not topology:
        return segs
    cur, start = None, None
    for i, c in enumerate(topology + "\0"):
        reg = sc.code_to_region(c) if c != "\0" else None
        reg = reg if reg in REGION_LABELS else None
        if reg != cur:
            if cur is not None:
                segs[cur].append((start + 1, i))
            cur, start = reg, i
    return segs


def region_alignments(model_seqs, target_seq, ref, topology):
    """
    For each membrane region, the segments defined on the reference species' topology, each
    aligned across every species exactly like the epitope: per-species coordinates, length,
    indels, identity and similarity to the target, and the residue-by-residue alignment.
    Identities here are reference-anchored, so they can differ slightly from the region
    table, which measures each species against the target with its own topology.
    """
    out = {"ref": ref, "regions": []}
    for reg, segs in topology_segments(topology).items():
        blocks = []
        for s_, e_ in segs:
            try:
                ep = epitope_across_species(model_seqs, target_seq, ref, (s_, e_))
            except ValueError:
                continue
            blocks.append({"range": [s_, e_], "rows": ep["rows"],
                           "alignment": ep["alignment"],
                           "substitutions": ep["substitutions"]})
        if blocks:
            out["regions"].append({"region": reg, "label": REGION_LABELS[reg],
                                   "segments": blocks})
    return out


def epitope_across_species(model_seqs, target_seq, epitope_ref, epitope_range):
    """
    Map an immunogen region (defined on epitope_ref, 1-based) onto every species and the
    target via pairwise alignment. Returns {'rows': [...], 'alignment': {species: string}}:
      rows       -- per-species coordinates, length, identity to the target, the number of
                    indel (gap) columns vs the target, and epitope coverage
      alignment  -- reference-anchored aligned residues of the epitope (one char per
                    reference column, '-' for a deletion) so differences are visible.

    Identity is computed over columns where BOTH species have a residue, which is the
    standard convention but means an insertion or deletion is excluded from the ratio
    rather than counted against it. An indel inside an epitope reshapes the local backbone
    and is usually more disruptive than a substitution, so `gaps_vs_target` and
    `epitope_coverage_pct` are reported alongside identity and must be read with it.

    Note the identity is reference-anchored (each species is projected onto the reference's
    epitope columns), not a direct pairwise model-vs-target alignment; where two species
    carry different indels this projection is what makes their columns comparable.
    """
    es, ee = epitope_range
    ref_seq = next((s for c, s in model_seqs if c == epitope_ref), None)
    if ref_seq is None:
        # Never fall back to another species: the coordinates would then be read against the
        # wrong sequence while the report still named the species the user asked for.
        have = ", ".join(c for c, _ in model_seqs) or "none"
        raise ValueError(
            f"the epitope reference '{epitope_ref}' is not among the species compared in this "
            f"run ({have}). It may be the target species itself, or its ortholog may have been "
            "skipped as a partial sequence. Choose a reference from that list, or pass the "
            "immunogen peptide with --epitope-seq and let AbScout find the species.")

    # Coordinates are only meaningful against a specific sequence. Silently returning an
    # empty epitope for an out-of-range span would look like "no conservation" rather than
    # "wrong numbering", which is the likely cause: datasheets sometimes number the mature
    # protein (signal peptide removed) while UniProt numbers from Met1.
    n_ref = len(ref_seq)
    if es < 1 or es > n_ref:
        raise ValueError(
            f"epitope start {es} is outside {epitope_ref} (1-{n_ref}). Check the numbering "
            "convention -- datasheets sometimes count from the mature protein (signal peptide "
            "removed) while UniProt counts from Met1. Passing the immunogen peptide with "
            "--epitope-seq avoids this entirely.")
    if ee > n_ref:
        sys.stderr.write(f"[warn] epitope end {ee} is past the end of {epitope_ref} "
                         f"({n_ref} aa); truncating to {n_ref}.\n")
        ee = n_ref
    if ee < es:
        raise ValueError(f"epitope range {es}-{ee} is empty (end before start).")
    cols = list(range(es - 1, ee))

    species = list(model_seqs) + [("target", target_seq)]
    dupes = {c for c, _ in species if [x for x, _ in species].count(c) > 1}
    if dupes:
        # Everything below is keyed by species name, so a repeated name would make one
        # entry's residues stand in for another's and report a wrong identity for both.
        raise ValueError(
            f"more than one sequence was supplied for {', '.join(sorted(dupes))}; the epitope "
            "tables are per species, so the caller must choose one ortholog per species first.")
    alignment, idx_map, maps = {}, {}, {}
    for common, seq in species:
        r2s = _ref_map(ref_seq, seq)
        maps[common] = r2s
        alignment[common] = "".join(seq[r2s[a]] if a in r2s else "-" for a in cols)
        idx_map[common] = [r2s[a] for a in cols if a in r2s]

    tgt = alignment.get("target", "")
    n_cols = len(cols)
    rows = []
    for common, _ in species:
        s_idx = idx_map[common]
        coords = f"{min(s_idx) + 1}-{max(s_idx) + 1}" if s_idx else "-"
        if common == "target":
            rows.append({"common": common, "coordinates": coords, "length": len(s_idx),
                         "identity_to_target": None, "similarity_to_target": None,
                         "gaps_vs_target": None, "radical_subs": None,
                         "epitope_coverage_pct": round(100 * len(s_idx) / max(n_cols, 1), 1)})
            continue
        gaps = sum(1 for a, b in zip(alignment[common], tgt)
                   if (a == "-") != (b == "-"))
        sub = sc.score_substitutions(list(zip(alignment[common], tgt)))
        rows.append({"common": common, "coordinates": coords, "length": len(s_idx),
                     "absent": not s_idx,     # no residue of this ortholog aligns to the epitope
                     "identity_to_target": sub["seq_identity_pct"],
                     "similarity_to_target": sub["similarity_pct"],
                     "radical_subs": sub["radical_subs"],
                     "mean_grantham_subs": sub["mean_grantham_subs"],
                     "gaps_vs_target": gaps,
                     "epitope_coverage_pct": round(100 * sub["compared"] / max(n_cols, 1), 1)})

    # residue-level detail for the comparison that actually decides binding: the species the
    # antibody was raised in, versus the target
    ref_aln = alignment.get(epitope_ref, "")
    subs = []
    for i, (a, b) in enumerate(zip(ref_aln, tgt)):
        if a == "-" or b == "-" or a == b:
            continue
        d = sc.grantham_distance(a, b)
        subs.append({"position": cols[i] + 1, "from": a, "to": b,
                     "grantham": (round(d) if d == d else None),
                     "class": sc.substitution_class(d),
                     "conservative": sc.is_conservative(a, b)})

    glyco = glyco_delta(ref_seq, target_seq, cols, maps.get("target", {}))
    return {"rows": rows, "alignment": alignment, "substitutions": subs, "glyco": glyco,
            "range_used": [es, ee]}             # may be shorter than requested (overhang cut)


def locate_peptide(peptide, ref_seq):
    """
    Find an immunogen peptide within the reference sequence. Returns
    (start, end, pct_identity) 1-based inclusive (exact match => 100.0), or None if the
    peptide is empty or longer than the reference. Falls back to the best-matching window
    when there's no exact hit (handles minor datasheet/isoform differences).
    """
    peptide = "".join(peptide.split()).upper()
    L = len(peptide)
    if L == 0 or L > len(ref_seq):
        return None
    i = ref_seq.find(peptide)
    if i >= 0:
        return (i + 1, i + L, 100.0)
    best_i, best_m = 0, -1
    for i in range(len(ref_seq) - L + 1):
        m = sum(a == b for a, b in zip(peptide, ref_seq[i:i + L]))
        if m > best_m:
            best_m, best_i = m, i
    return (best_i + 1, best_i + L, round(100.0 * best_m / L, 1))


MIN_PEPTIDE_MATCH = 70.0    # below this, a "best window" is noise, not a located immunogen

# ----------------------------------------------------------------------------------
# Is this sequence really the ortholog?
#
# UniProt's gene-symbol search matches synonyms. Asking for human "CDH1" returned Q9UM11,
# which is FZR1 -- a cell-cycle protein that lists CDH1 among its synonyms -- not
# E-cadherin. Any step that picks a sequence by gene symbol can therefore pull in an
# unrelated protein, so every such step is checked against the other species' sequences.
#
# The score is the fraction of another species' sequence that this one reproduces
# (identical aligned residues / that sequence's length), taking the best over the other
# species in the run. A real ortholog scores well against at least one other ortholog; a
# protein that merely shares a symbol scores low against all of them; a fragment scores
# low because it reproduces only part of the full-length sequence. Calibrated on UniProt:
# real orthologs here 0.45-0.62 (CD3E evolves fast: human-mouse is 0.59), reviewed dog and
# rabbit CD3E 0.60 and 0.62 against 0.45 and 0.51 for the TrEMBL entries the xref gave,
# FZR1 against E-cadherin 0.14.
# ----------------------------------------------------------------------------------
ORTHOLOG_FLOOR = 0.25       # below this against every other species: probably not the ortholog
SWAP_TOLERANCE = 0.02       # a reviewed entry may be this much less similar and still win


def _identical_aligned(a, b, memo):
    """Identical residues in the global alignment of a and b (memoised, order-independent)."""
    key = (a, b) if a <= b else (b, a)
    if key not in memo:
        x, y = key
        pairs, _ = sc.correspondence_seq({"seq": x}, {"seq": y})
        memo[key] = sum(x[i] == y[j] for i, j in pairs)
    return memo[key]


def ortholog_score(seq, others, memo):
    """Best fraction of another species' sequence that `seq` reproduces (0-1)."""
    best = 0.0
    for o in others:
        if o and seq:
            best = max(best, _identical_aligned(seq, o, memo) / len(o))
    return best


def reviewed_alternatives(gene, taxid):
    """Reviewed (Swiss-Prot) accessions for gene+species. The symbol search matches synonyms,
    so callers must check that a result is the same protein before using it."""
    key = f"{gene}_{taxid}"
    cached = cache_get("reviewed_alts", key)
    if cached is not None:
        return [a for a in cached.split(",") if a]
    results = ro._uniprot_search({
        "query": f"(gene_exact:{gene}) AND (organism_id:{taxid}) AND (reviewed:true)",
        "fields": "accession", "format": "json", "size": "10"})
    accs = sorted({r.get("primaryAccession") for r in results if r.get("primaryAccession")})
    cache_put("reviewed_alts", key, ",".join(accs))
    return accs


class AbScoutError(Exception):
    """A run cannot proceed (e.g. no ortholog found). An ordinary exception rather than
    SystemExit so that batch mode can record it for one case and carry on with the rest."""


def detect_epitope_from_peptide(model_seqs, peptide, min_match=MIN_PEPTIDE_MATCH, notes=None):
    """
    Locate an immunogen peptide across all model species and return the best match as
    (species, start, end, match_pct). The best-matching species is the immunogen species,
    so the user does not need to know it.

    locate_peptide() always returns its best-scoring window, so without a floor an absent
    or mistyped peptide would still yield a confident-looking position (e.g. 5% match).
    Returns None when the best match falls below `min_match`. The reason is written to
    stderr and, when a `notes` list is passed, appended to it so it travels with the result.
    """
    def say(level, short, detail=""):
        if notes is not None:
            notes.append(short)
        sys.stderr.write(f"[{level}] {short}{detail}\n")

    best = None
    for common, seq in model_seqs:
        loc = locate_peptide(peptide, seq)
        if loc and (best is None or loc[2] > best[3]):
            best = (common, loc[0], loc[1], loc[2])
    if best is None:
        say("warn", "immunogen peptide not located: it is empty or longer than every model "
                    "sequence; no epitope reported")
        return None
    if best[3] < min_match:
        say("warn", f"immunogen peptide not located: best match {best[0]} {best[3]:.0f}%, below "
                    f"the {min_match:.0f}% floor; no epitope reported",
            ". Check for a typo, a truncated paste, tag/linker residues, or an immunogen from a "
            "species outside human/mouse/rat/rabbit.")
        return None
    if best[3] < 100.0:
        say("note", f"immunogen peptide matches {best[0]} at {best[3]:.0f}%, not exactly; "
                    "located coordinates are approximate",
            " -- it may come from a different isoform or carry non-native residues "
            "(for example a linker cysteine).")
    return best


def analyze(gene, target_taxid, work, epitope=None, epitope_ref="human", epitope_seq=None):
    resolved = ro.resolve(gene, target_taxid)
    orths = resolved["orthologs"]
    targets = [o for o in orths if o["is_target"] and o["uniprot"]]
    models = [o for o in orths if not o["is_target"] and o["uniprot"]]
    if not targets:
        raise AbScoutError(
            f"No Compara ortholog for {gene} in taxid {target_taxid} -- please check it "
            f"manually. AbScout does not fall back to matching on gene name, because a gene "
            f"with the same symbol is not necessarily the ortholog and a wrong sequence here "
            f"would silently distort every identity in the run. If Ensembl timed out above, "
            f"the family may just be slow: retry with --timeout 300.")
    notes = []                  # everything that qualifies the numbers; travels to JSON / CSV

    # ---- the target sequence --------------------------------------------------------------
    # Every identity in the run is measured against this one sequence, so taking whichever
    # ortholog Ensembl happened to list first was the weakest link in the pipeline: mouse
    # CDH1 gave both an 884 aa entry and a 493 aa fragment, and the fragment won by position
    # alone. Rank the candidates exactly as the models are ranked.
    def rank_against(longest):
        # Completeness outranks curation. Swiss-Prot carries short entries -- mouse CDH1 is
        # reviewed at 493 aa while the full 884 aa entry is unreviewed -- so ranking on
        # review status first would pick the fragment and deflate every identity in the run.
        def rank(item):                         # higher is better; fully deterministic
            o, seq = item
            return (len(seq) >= 0.6 * longest,          # not truncated vs the best available
                    o.get("uniprot_status") == "reviewed",
                    o.get("orthology_type") == "ortholog_one2one",
                    len(seq),
                    o["uniprot"])                       # stable tie-break
        return rank

    target_cands, unusable = [], []
    for t in targets:
        try:
            target_cands.append((t, fetch_uniprot_seq(t["uniprot"])))
        except Exception as exc:
            sys.stderr.write(f"  [warn] could not fetch target {t['uniprot']}: {exc}\n")
            unusable.append(f"{t['uniprot']} ({t['species']}): {exc}")
    if not target_cands:
        raise AbScoutError(
            f"No usable target sequence for {gene} in taxid {target_taxid}. Ensembl Compara "
            f"found the ortholog, but its UniProt entry cannot be used -- "
            + "; ".join(unusable) + ". Please check it manually.")
    target_cands.sort(key=rank_against(max(len(x) for _, x in target_cands)), reverse=True)
    target, target_seq = target_cands[0]
    if len(target_cands) > 1:
        others = ", ".join(f"{o['uniprot']} ({len(s)} aa)" for o, s in target_cands[1:])
        sys.stderr.write(f"  [note] target: {len(target_cands)} orthologs; using "
                         f"{target['uniprot']} ({len(target_seq)} aa), not {others}\n")
        notes.append(f"target: {len(target_cands)} orthologs returned, used "
                     f"{target['uniprot']} ({len(target_seq)} aa) over {others}")

    sys.stderr.write(f"[abscout] target {target['uniprot']} ({target['species']}); "
                     f"{len(models)} model orthologs to compare\n")

    # ---- one ortholog per model species -------------------------------------------------
    # Compara can return several orthologs for one species (CD3E gave two rabbit entries,
    # CDH1 two mouse). Keeping them all is wrong three ways: "will a rabbit antibody work"
    # has one answer per species, the duplicates run TMbed twice, and downstream the epitope
    # tables are keyed by species name so the second entry silently overwrote the first.
    # Pick the best representative deterministically and record what was set aside.
    with_seq = []
    for m in models:
        try:
            with_seq.append((m, fetch_uniprot_seq(m["uniprot"])))
        except Exception as exc:
            sys.stderr.write(f"  [warn] could not fetch {m['common']} ({m['uniprot']}): {exc}\n")
            notes.append(f"{m['common']} ({m['uniprot']}) skipped: {str(exc)[:100]}")

    by_species, chosen = {}, []
    for m, seq in with_seq:
        by_species.setdefault(m["common"], []).append((m, seq))
    for common, cands in by_species.items():
        # Judged against the best entry for this same species, not against the target: if the
        # target is itself short, every candidate would pass a target-relative test and the
        # ranking would fall through to review status and pick a fragment.
        cands.sort(key=rank_against(max(len(x) for _, x in cands)), reverse=True)
        chosen.append(cands[0])
        if len(cands) > 1:
            others = ", ".join(f"{o['uniprot']} ({len(s)} aa)" for o, s in cands[1:])
            best = cands[0][0]["uniprot"]
            sys.stderr.write(f"  [note] {common}: {len(cands)} orthologs; using {best} "
                             f"({len(cands[0][1])} aa), not {others}\n")
            notes.append(f"{common}: {len(cands)} orthologs returned, used {best} "
                         f"({len(cands[0][1])} aa) over {others}")

    memo = {}                   # alignment results shared by every identity check below

    def others_than(common):
        """Sequences of every other species in the run (models and target)."""
        out = [x for mm, x in chosen if mm["common"] != common]
        if common != "__target__":
            out.append(target_seq)
        return out

    # ---- rescue partial models before anything is derived from them ----------------------
    # A model much shorter than the target is a partial entry: its regions would be missing
    # and its identity deflated. Swiss-Prot carries short entries too (human CDH1 resolved to
    # a 496 aa reviewed accession against an 882 aa protein), so before giving up, look for a
    # longer entry for this gene and species. Only reached when the length already looks
    # wrong, and a reviewed accession is only ever swapped for another reviewed one.
    rescued = []
    for m, mseq in chosen:
        if len(mseq) < 0.6 * len(target_seq):
            try:
                better, bstatus, _ = ro.upgrade_partial(
                    m["uniprot"], m.get("uniprot_status"), gene, m["taxid"])
                if better and better != m["uniprot"]:
                    bseq = fetch_uniprot_seq(better)
                    oth = others_than(m["common"])
                    same_protein = (ortholog_score(bseq, oth, memo) >= ORTHOLOG_FLOOR and
                                    ortholog_score(bseq, oth, memo) >
                                    ortholog_score(mseq, oth, memo))
                    if not same_protein:
                        sys.stderr.write(f"  [{m['common']}] longer entry {better} found by gene "
                                         f"symbol does not resemble the other species; not "
                                         f"used\n")
                        notes.append(f"{m['common']}: longer entry {better} found by gene symbol "
                                     "was not used -- it does not resemble the other species")
                    elif len(bseq) > len(mseq):
                        sys.stderr.write(f"  [{m['common']}] {m['uniprot']} is {len(mseq)} aa "
                                         f"vs target {len(target_seq)} aa -> using {better} "
                                         f"({len(bseq)} aa)\n")
                        notes.append(f"{m['common']}: {m['uniprot']} was partial ({len(mseq)} aa "
                                     f"vs target {len(target_seq)} aa), used {better} "
                                     f"({len(bseq)} aa)")
                        m, mseq = dict(m, uniprot=better, uniprot_status=bstatus), bseq
            except Exception as exc:
                sys.stderr.write(f"  [warn] could not check {m['common']} for a longer "
                                 f"entry: {exc}\n")
        rescued.append((m, mseq))
    chosen = rescued

    # ---- is the target itself a fragment? -------------------------------------------------
    # Judged against the orthologs just fetched rather than an extra lookup, so this costs
    # nothing unless the target really does look short. A truncated target is the most
    # damaging failure there is: it deflates every identity, drops whole regions, and makes
    # full-length models look partial by comparison.
    if chosen:
        typical = sorted(len(s) for _, s in chosen)[len(chosen) // 2]       # median model
        if len(target_seq) < 0.6 * typical:
            better, bstatus, _ = ro.upgrade_partial(
                target["uniprot"], target.get("uniprot_status"), gene, target["taxid"])
            if better and better != target["uniprot"]:
                bseq = fetch_uniprot_seq(better)
                oth = others_than("__target__")
                if (ortholog_score(bseq, oth, memo) < ORTHOLOG_FLOOR or
                        ortholog_score(bseq, oth, memo) <=
                        ortholog_score(target_seq, oth, memo)):
                    notes.append(f"target: longer entry {better} found by gene symbol was not "
                                 "used -- it does not resemble the other species")
                    bseq = ""                    # disqualify
                if bseq and len(bseq) > len(target_seq):
                    sys.stderr.write(f"  [target] {target['uniprot']} is {len(target_seq)} aa "
                                     f"vs {typical} aa for the orthologs -> using {better} "
                                     f"({len(bseq)} aa)\n")
                    notes.append(f"target: {target['uniprot']} was partial ({len(target_seq)} aa "
                                 f"vs {typical} aa for the orthologs), used {better} "
                                 f"({len(bseq)} aa)")
                    target = dict(target, uniprot=better, uniprot_status=bstatus)
                    target_seq = bseq
            if len(target_seq) < 0.6 * typical:
                notes.append(f"target {target['uniprot']} is {len(target_seq)} aa against "
                             f"{typical} aa for the orthologs and no longer entry was found: "
                             "identities are measured against a possibly partial sequence")
                sys.stderr.write(f"  [warn] target {target['uniprot']} looks partial "
                                 f"({len(target_seq)} aa vs {typical} aa); no longer entry "
                                 "found -- identities may be understated\n")

    # ---- prefer curated (Swiss-Prot) entries --------------------------------------------
    # Compara's cross-reference sometimes points at an unreviewed TrEMBL model even when a
    # curated entry exists for the same gene and species. Those models can be wrong exactly
    # where it matters: dog CD3E A0A8I3MGF0 has a mispredicted C-terminus that drops the
    # conserved ERPPPVPNPDYEP motif (30.8% epitope identity) where the reviewed P27597 has it
    # (92.3%); rabbit A0A5F9CJ85 likewise. A reviewed entry replaces the unreviewed one only
    # if it resembles the other species at least as well -- the gene-symbol search can return
    # a different protein entirely, and that check is what keeps it out.
    def prefer_reviewed(entry, seq, common):
        if entry.get("uniprot_status") == "reviewed":
            return entry, seq
        try:
            alts = [a for a in reviewed_alternatives(gene, entry["taxid"])
                    if a != entry["uniprot"]]
        except Exception as exc:
            sys.stderr.write(f"  [warn] could not look up reviewed entries for {common}: "
                             f"{exc}\n")
            return entry, seq
        if not alts:
            return entry, seq
        oth = others_than(common)
        cur = ortholog_score(seq, oth, memo)
        best = None
        for acc in alts:
            try:
                aseq = fetch_uniprot_seq(acc)
            except Exception:
                continue
            sc_ = ortholog_score(aseq, oth, memo)
            if best is None or sc_ > best[2]:
                best = (acc, aseq, sc_)
        if best and best[2] >= ORTHOLOG_FLOOR and best[2] >= cur - SWAP_TOLERANCE:
            label = "target" if common == "__target__" else common
            sys.stderr.write(f"  [{label}] using reviewed {best[0]} ({len(best[1])} aa) instead "
                             f"of unreviewed {entry['uniprot']} ({len(seq)} aa)\n")
            notes.append(f"{label}: used reviewed {best[0]} instead of unreviewed "
                         f"{entry['uniprot']} (similarity to other species {best[2]:.2f} vs "
                         f"{cur:.2f})")
            return dict(entry, uniprot=best[0], uniprot_status="reviewed"), best[1]
        if best:
            label = "target" if common == "__target__" else common
            notes.append(f"{label}: reviewed {best[0]} exists but resembles the other species "
                         f"less ({best[2]:.2f} vs {cur:.2f}); kept {entry['uniprot']}")
        return entry, seq

    target, target_seq = prefer_reviewed(target, target_seq, "__target__")
    chosen = [prefer_reviewed(m, x, m["common"]) for m, x in chosen]

    # ---- is every sequence plausibly the ortholog? --------------------------------------
    # A last check that does not depend on how a sequence was found: one that resembles
    # none of the other species is flagged for manual checking rather than trusted.
    suspect = set()
    for m, x in chosen:
        if ortholog_score(x, others_than(m["common"]), memo) < ORTHOLOG_FLOOR:
            suspect.add(m["common"])
            notes.append(f"{m['common']} ({m['uniprot']}) resembles none of the other species "
                         f"(best {ortholog_score(x, others_than(m['common']), memo):.2f}) -- it "
                         "may not be the ortholog; please check it manually")
    if chosen and ortholog_score(target_seq, others_than("__target__"), memo) < ORTHOLOG_FLOOR:
        notes.append(f"target ({target['uniprot']}) resembles none of the model species "
                     f"(best {ortholog_score(target_seq, others_than('__target__'), memo):.2f}) "
                     "-- it may not be the ortholog; please check it manually")
    for name in sorted(suspect):
        sys.stderr.write(f"  [warn] {name} may not be the ortholog -- please check it manually\n")

    # ---- topology ------------------------------------------------------------------------
    # Decide membrane-or-not from the best evidence available, not from whichever ortholog
    # came back first. Curated annotation outranks prediction, and TMbed's "no TM" is only
    # allowed to settle the question when it ran on the longest model: a fragment that
    # happens to exclude the helix previously marked the whole protein soluble.
    curated, anchor = {}, {}
    for m, mseq in chosen:
        curated[m["uniprot"]], anchor[m["uniprot"]] = uniprot_topology(m["uniprot"], len(mseq))
    # A curated helix or a curated GPI anchor both settle the question: the protein sits in
    # a membrane, and TMbed must not be allowed to call it soluble.
    has_curated_tm = any(k in ("tm", "gpi") for k in anchor.values())
    tmbed_ok = tmbed_available()
    # Without TMbed and without a curated topology anywhere, nothing can say where the
    # membrane is. That is "unknown", never "soluble".
    unresolved = not has_curated_tm and not tmbed_ok
    no_topology = []            # species compared on whole-protein identity only
    tmbed_cache, no_tm_seen = {}, False
    if not has_curated_tm and tmbed_ok:
        needs = [(m, s) for m, s in chosen if curated.get(m["uniprot"]) is None]
        if needs:
            m0, s0 = max(needs, key=lambda it: len(it[1]))       # longest, not first
            try:
                t0 = run_tmbed(s0, work, m0["common"], acc=m0["uniprot"])
                tmbed_cache[m0["uniprot"]] = t0
                if t0 and not any(c in t0 for c in "BbHh"):
                    no_tm_seen = True
            except Exception as exc:
                sys.stderr.write(f"  [warn] TMbed failed on {m0['common']}: {exc}\n")

    model_results, model_seqs = [], []
    topo_by_common = {}         # per-species topology, kept for the signal-peptide check
    for m, mseq in chosen:
        try:
            if len(mseq) < 0.6 * len(target_seq):
                sys.stderr.write(
                    f"  [warn] {m['common']} ({m['uniprot']}) is {len(mseq)} aa vs target "
                    f"{len(target_seq)} aa -- looks like a partial sequence, and no longer "
                    f"entry was found; skipping\n")
                notes.append(f"{m['common']} ({m['uniprot']}) skipped: partial sequence, "
                             f"{len(mseq)} aa vs target {len(target_seq)} aa, no longer entry")
                continue
            if m["uniprot"] not in curated:        # e.g. an accession swapped in after the pass
                curated[m["uniprot"]], anchor[m["uniprot"]] = uniprot_topology(
                    m["uniprot"], len(mseq))
            topo = curated[m["uniprot"]]                              # curated first
            kind = anchor.get(m["uniprot"], "")
            source = "UniProt (GPI-anchored)" if kind == "gpi" else "UniProt"
            if topo is None and m["uniprot"] in tmbed_cache:
                topo, source = tmbed_cache[m["uniprot"]], "TMbed"
            elif topo is None and no_tm_seen:
                source = "none (soluble)"                             # settled above
            elif topo is None and not tmbed_ok:                       # nothing to predict with
                source = "not annotated in UniProt"
                no_topology.append(f"{m['common']} ({m['uniprot']})")
            elif topo is None:                                        # fall back to prediction
                topo = run_tmbed(mseq, work, m["common"], acc=m["uniprot"])
                source = "TMbed"
            if topo is None and not tmbed_ok and source != "none (soluble)":
                sys.stderr.write(f"  [{m['common']}] no curated topology in UniProt and TMbed "
                                 f"is not available: whole-protein identity only\n")
            else:
                sys.stderr.write(f"  [{m['common']}] reference topology from {source}\n")
            topo_by_common[m["common"]] = topo
            res = sc.compare_sequences(mseq, target_seq, ref_topology=topo)
            regions = {k: res["regions"].get(k, {}).get("seq_identity_pct")
                       for k in ("extracellular", "TM", "intracellular")}
            model_results.append({"species": m["species"], "common": m["common"],
                                  "uniprot": m["uniprot"], "topology_source": source,
                                  "regions": regions,
                                  "anchor": anchor.get(m["uniprot"], ""),
                                  "region_detail": res["regions"],
                                  "overall": res["global"]["seq_identity_pct"],
                                  "overall_similarity": res["global"]["similarity_pct"],
                                  "radical_subs": res["global"]["radical_subs"],
                                  "coverage_pct": res["global"]["coverage_pct"],
                                  "length": len(mseq)})
            model_seqs.append((m["common"], mseq))
        except Exception as exc:                       # skip a model rather than abort the run
            sys.stderr.write(f"  [warn] skipping {m['common']} ({m['uniprot']}): {exc}\n")
            notes.append(f"{m['common']} ({m['uniprot']}) skipped: {str(exc)[:120]}")

    # resolve the epitope range: from a pasted immunogen peptide (auto-located) or explicit coords
    epitope_range, located = epitope, None
    if epitope_seq and model_seqs:
        best = detect_epitope_from_peptide(model_seqs, epitope_seq, notes=notes)
        if best:
            common, s, e, pct = best
            if epitope_ref and epitope_ref != "human" and epitope_ref != common:
                sys.stderr.write(
                    f"[note] --epitope-ref {epitope_ref} was given, but the peptide matches "
                    f"{common} better ({pct:.0f}%); using {common} as the epitope reference.\n")
                notes.append(f"--epitope-ref {epitope_ref} overridden: peptide matches {common}")
            epitope_ref = common                       # detected immunogen species = reference
            epitope_range = (s, e)
            located = {"match_identity": pct, "detected_species": common,
                       "peptide_len": len("".join(epitope_seq.split()))}
            sys.stderr.write(f"[abscout] immunogen peptide best matches {common} at residues "
                             f"{s}-{e} ({pct:.0f}% match); using {common} as the epitope reference\n")

    # A protein with no transmembrane segment in ANY model is soluble (e.g. ZO-1/TJP1, a
    # cytoplasmic scaffold). The extracellular/TM/intracellular split describes which side of
    # the bilayer a residue sits on, so it is meaningless here: report overall identity only.
    is_membrane = any(m["regions"].get("TM") is not None or m.get("anchor") == "gpi"
                      for m in model_results)
    if model_results and not is_membrane and not unresolved:
        for m in model_results:
            m["regions"] = {"extracellular": None, "TM": None, "intracellular": None}
            m["soluble"] = True
        notes.append("no transmembrane segment detected: reported as whole-protein identity")
    elif model_results and unresolved:
        notes.append("membrane topology unknown: no species has a curated topology in UniProt "
                     "and TMbed is not available, so only whole-protein identity is shown (the "
                     "protein is not assumed to be soluble); " + _tmbed_hint())


    epitope_info = None
    if epitope_range and model_seqs:
        ep = epitope_across_species(model_seqs, target_seq, epitope_ref, epitope_range)
        es, ee = ep["range_used"]
        epitope_info = {"ref": epitope_ref, "range": [es, ee],
                        "rows": ep["rows"], "alignment": ep["alignment"],
                        "substitutions": ep["substitutions"], "glyco": ep["glyco"]}
        if [es, ee] != list(epitope_range):
            epitope_info["requested_range"] = list(epitope_range)
            notes.append(f"epitope range {epitope_range[0]}-{epitope_range[1]} runs past the end "
                         f"of {epitope_ref}; cut to {es}-{ee}")
        if located:
            epitope_info["located"] = located
        else:
            # coordinates were supplied by hand: surface the numbering convention, because
            # UniProt counts from Met1 while datasheets sometimes count the mature protein
            ref_entry = next((m for m in model_results if m["common"] == epitope_ref), None)
            ref_seq = next((s for c, s in model_seqs if c == epitope_ref), "")
            n_ref = len(ref_seq)
            nterm = uniprot_nterm(ref_entry["uniprot"], n_ref) if ref_entry else None
            if not nterm or nterm["mature_start"] == 1:
                # nothing curated: fall back to a signal peptide predicted in the topology
                sp = signal_peptide_len(topo_by_common.get(epitope_ref))
                if sp:
                    nterm = {"signal_len": sp, "mature_start": sp + 1,
                             "parts": [f"{sp}-residue signal peptide (predicted)"]}
                elif not nterm:
                    nterm = {"signal_len": 0, "mature_start": 1, "parts": []}
            off = nterm["mature_start"] - 1

            def clip(s, n=60):
                return s if len(s) <= n else s[:n] + "..."

            numbering = {
                "signal_peptide_len": nterm["signal_len"],
                "mature_start": nterm["mature_start"],
                "processing": nterm["parts"],
                "residues_at_range": clip(ref_seq[es - 1:ee]),
                "mature_equivalent": None, "mature_fits": None, "mature_residues": None,
                "reference_length": n_ref,
                "reference_uniprot": (ref_entry or {}).get("uniprot"),
            }
            if off > 0:
                a, b = es + off, ee + off
                numbering["mature_equivalent"] = [a, b]
                numbering["mature_fits"] = b <= n_ref
                if b <= n_ref:
                    numbering["mature_residues"] = clip(ref_seq[a - 1:b])
            epitope_info["numbering"] = numbering
            if numbering["mature_fits"]:
                a, b = numbering["mature_equivalent"]
                notes.append(f"numbering: {es}-{ee} read as UniProt numbering "
                             f"({numbering['residues_at_range']}); if the datasheet numbers the "
                             f"mature protein the range is {a}-{b} "
                             f"({numbering['mature_residues']})")

    low = [m for m in model_results if (m.get("coverage_pct") or 100) < 80]
    if low:
        notes.append("low alignment coverage: " +
                     ", ".join(f"{m['common']} {m['coverage_pct']:.0f}%" for m in low))
    region_info = None
    if model_results and is_membrane:
        anchor_sp = next((c for c in ["human"] + [m["common"] for m in model_results]
                          if topo_by_common.get(c) and c in dict(model_seqs)), None)
        if anchor_sp:
            region_info = region_alignments(model_seqs, target_seq, anchor_sp,
                                            topo_by_common[anchor_sp])
    resolved["region_alignments"] = region_info
    if no_topology and not unresolved:
        notes.append("no curated membrane topology in UniProt for " + ", ".join(no_topology)
                     + ": whole-protein identity only for "
                     + ("this species" if len(no_topology) == 1 else "these species")
                     + (" (its residues are still shown in the region-by-region alignments)"
                        if region_info else "") + "; " + _tmbed_hint())

    resolved["notes"] = notes
    return resolved, target, targets, model_results, epitope_info


def _fmt(x):
    return f"{x:.1f}" if x is not None else "  -"


def render_html(gene, target, model_results, epitope_info=None, notes=None, regions=None):
    """Self-contained HTML report (no server, opens in a browser)."""
    order = {c: i for i, c in enumerate(MODEL_ORDER)}
    rows = sorted(model_results, key=lambda m: order.get(m["common"], 99))
    soluble = bool(model_results) and all(m.get("soluble") for m in model_results)

    def cell(pct):
        if pct is None:
            return '<td style="background:#eef0f2;color:#8a919a">n/a</td>'
        c = "#1b9e57" if pct >= 85 else "#e0972a" if pct >= 70 else "#cf3b3b"
        return f'<td style="background:{c};color:#fff;font-weight:600">{pct:.1f}%</td>'

    if soluble:
        body = ""
        for m in rows:
            body += (f'<tr><td class="lbl">{m["common"].capitalize()}'
                     f'<span class="acc">{m["uniprot"]}</span></td>{cell(m["overall"])}</tr>')
        note = ('<p class="sub">No transmembrane segment detected &mdash; this appears to be a '
                'soluble (non-membrane) protein, so the extracellular / transmembrane / '
                'intracellular split does not apply.</p>')
        main_table = (f'{note}<table><tr><th>Model (antibody host)</th>'
                      f'<th>Whole-protein identity</th></tr>{body}</table>')
    else:
        body = ""
        for m in rows:
            r = m["regions"]
            body += (f'<tr><td class="lbl">{m["common"].capitalize()}'
                     f'<span class="acc">{m["uniprot"]} · topology: {m.get("topology_source","?")}</span></td>'
                     f'{cell(r.get("extracellular"))}{cell(r.get("TM"))}'
                     f'{cell(r.get("intracellular"))}{cell(m["overall"])}</tr>')
        main_table = ('<table><tr><th>Model (antibody host)</th><th>Extracellular</th>'
                      '<th>Transmembrane</th><th>Intracellular</th><th>Overall</th></tr>'
                      f'{body}</table>')

    ep_html = ""
    if epitope_info:
        er = epitope_info["range"]
        erows = ""
        for row in epitope_info["rows"]:
            ident = row["identity_to_target"]
            id_cell = (cell(ident) if ident is not None
                       else '<td style="background:#eef0f2;color:#8a919a">&mdash;</td>')
            if row.get("absent") and row["common"] != "target":
                id_cell = ('<td style="background:#eef0f2;color:#6b7280" title="No residue of this '
                           'ortholog aligns to the epitope; the entry is probably truncated '
                           'there. This is not 0%.">absent</td>')
            name = "Target" if row["common"] == "target" else row["common"].capitalize()
            gaps = row.get("gaps_vs_target")
            gcell = ('<td style="background:#eef0f2;color:#8a919a">&mdash;</td>' if gaps is None
                     else (f'<td style="color:#b91c1c;font-weight:700">{gaps}</td>' if gaps
                           else '<td>0</td>'))
            sim = row.get("similarity_to_target")
            scell = (cell(sim) if sim is not None
                     else '<td style="background:#eef0f2;color:#8a919a">&mdash;</td>')
            rad = row.get("radical_subs")
            rcell = ('<td style="background:#eef0f2;color:#8a919a">&mdash;</td>' if rad is None
                     else (f'<td style="color:#b91c1c;font-weight:700">{rad}</td>' if rad
                           else '<td>0</td>'))
            erows += (f'<tr><td class="lbl">{name}</td><td>{row["coordinates"]}</td>'
                      f'<td>{row["length"]}</td>{gcell}{id_cell}{scell}{rcell}</tr>')
        loc = epitope_info.get("located")
        loc_note = (f' &middot; located from a pasted peptide ({loc["match_identity"]:.0f}% match)'
                    if loc else '')
        rq = epitope_info.get("requested_range")
        if rq:
            loc_note += (f' &middot; requested {rq[0]}&ndash;{rq[1]}, cut to the '
                         f'protein&rsquo;s length')

        subs = epitope_info.get("substitutions") or []
        subs_html = ""
        if subs:
            srows = ""
            for s in subs:
                g = "&mdash;" if s["grantham"] is None else f'{s["grantham"]:.0f}'
                colour = {"radical": "#b91c1c", "moderate": "#b45309",
                          "conservative": "#15803d"}.get(s["class"], "#6b7280")
                srows += (f'<tr><td>{s["position"]}</td>'
                          f'<td><b>{s["from"]}&rarr;{s["to"]}</b></td><td>{g}</td>'
                          f'<td style="color:{colour};font-weight:600">{s["class"]}</td></tr>')
            subs_html = (f'<h2 style="margin-top:16px">Substitutions '
                         f'({epitope_info["ref"]} &rarr; target)</h2>'
                         f'<table><tr><th>Position</th><th>Change</th><th>Grantham</th>'
                         f'<th>Class</th></tr>{srows}</table>')

        gl = epitope_info.get("glyco") or {}
        glyco_html = ""
        if gl.get("gained_in_target") or gl.get("lost_in_target"):
            bits = []
            if gl.get("gained_in_target"):
                bits.append("<b>gained</b> in target at " +
                            ", ".join(str(p) for p in gl["gained_in_target"]))
            if gl.get("lost_in_target"):
                bits.append("<b>lost</b> in target at " +
                            ", ".join(str(p) for p in gl["lost_in_target"]))
            glyco_html = (f'<p class="sub" style="color:#b45309">N-glycosylation sequon '
                          f'(N-X-S/T) {"; ".join(bits)} ({epitope_info["ref"]} numbering). '
                          f'A glycan can physically block antibody access, so this can change '
                          f'reactivity independently of sequence identity.</p>')

        nb = epitope_info.get("numbering")
        numbering_html = ""
        if nb:
            ref = epitope_info["ref"]
            if nb.get("mature_equivalent"):
                why = " + ".join(nb.get("processing") or []) or "N-terminal processing"
                a, b = nb["mature_equivalent"]
                if nb.get("mature_fits"):
                    mature = (f' {ref.capitalize()}&rsquo;s mature chain starts at residue '
                              f'{nb["mature_start"]} ({why}). If the datasheet numbered the '
                              f'<b>mature</b> protein, the same range is UniProt {a}&ndash;{b}: '
                              f'<code>{nb["mature_residues"]}</code>. Compare both with the '
                              f'datasheet; if the second matches, re-run with that range.')
                else:
                    mature = (f' {ref.capitalize()}&rsquo;s mature chain starts at residue '
                              f'{nb["mature_start"]} ({why}). Read as mature-protein numbering '
                              f'the range would be UniProt {a}&ndash;{b}, which runs past the '
                              f'end of the protein, so it can only be UniProt numbering.')
            else:
                mature = (f' {ref.capitalize()} has no signal peptide or propeptide annotated, '
                          f'so mature-protein and UniProt numbering coincide here.')
            box = ' numbox' if nb.get("mature_fits") else ''
            numbering_html = (
                f'<p class="sub{box}"><b>Numbering check.</b> {er[0]}&ndash;{er[1]} read as '
                f'{ref} UniProt numbering (from Met1) of {nb["reference_uniprot"]}, '
                f'{nb["reference_length"]} aa. Residues there: '
                f'<code>{nb["residues_at_range"]}</code>.{mature}</p>')
        aln = epitope_info.get("alignment") or {}
        aln_html = ""
        if aln:
            tgt = aln.get("target", "")
            names = ([n for n in MODEL_ORDER if n in aln]
                     + [n for n in aln if n not in MODEL_ORDER and n != "target"]
                     + (["target"] if "target" in aln else []))
            arows = ""
            for n in names:
                s = aln[n]
                spans = ""
                for i, ch in enumerate(s):
                    diff = (n != "target" and i < len(tgt) and ch != "-"
                            and tgt[i] != "-" and ch != tgt[i])
                    spans += f'<span class="d">{ch}</span>' if diff else ch
                label = "target" if n == "target" else n
                arows += f'<div><span class="al">{label}</span>{spans}</div>'
            aln_html = (f'<h2 style="margin-top:16px">Epitope alignment</h2>'
                        f'<p class="sub">Red = residue differs from the target</p>'
                        f'<div class="aln">{arows}</div>')
        ep_html = (f'<h2 style="margin-top:20px">Antibody epitope region</h2>'
                   f'<p class="sub">Defined on {epitope_info["ref"]}, residues {er[0]}&ndash;{er[1]}'
                   f'{loc_note} &middot; identity is each species vs the target over this region</p>'
                   f'{numbering_html}'
                   f'<table><tr><th>Species</th><th>Coordinates</th><th>Length</th>'
                   f'<th>Indels</th><th>Identity</th><th>Similarity</th><th>Radical subs</th>'
                   f'</tr>{erows}</table>'
                   f'<p class="sub">Identity = exact matches. Similarity = identical or '
                   f'conservative (BLOSUM62 &gt; 0). Radical = substitutions with Grantham '
                   f'&gt; 100, i.e. a substantial change in charge, size or polarity.</p>'
                   f'{subs_html}{glyco_html}{aln_html}')

    import html as _html

    def wrapped_alignment(aln, ref_name, ref_start, width=60):
        """Alignment in blocks of `width` columns, numbered on the reference, with each residue
        that differs from the target highlighted (same convention as the epitope view)."""
        tgt = aln.get("target", "")
        names = ([n for n in MODEL_ORDER if n in aln]
                 + [n for n in aln if n not in MODEL_ORDER and n != "target"]
                 + (["target"] if "target" in aln else []))
        n = len(next(iter(aln.values()), ""))
        ref_aln = aln.get(ref_name, "")
        out, ref_pos = "", ref_start
        for b in range(0, n, width):
            out += f'<div class="pos">{ref_name} {ref_pos}</div>'
            for nm in names:
                seg = aln[nm][b:b + width]
                cells = ""
                for k, ch in enumerate(seg):
                    i = b + k
                    diff = (nm != "target" and ch != "-" and i < len(tgt)
                            and tgt[i] != "-" and ch != tgt[i])
                    cells += f'<span class="d">{ch}</span>' if diff else ch
                out += f'<div><span class="al">{nm}</span>{cells}</div>'
            ref_pos += sum(1 for ch in ref_aln[b:b + width] if ch != "-")
        return f'<div class="aln">{out}</div>'

    regions_html = ""
    if regions and regions.get("regions"):
        ref = regions["ref"]
        parts = []
        for reg in regions["regions"]:
            segs = reg["segments"]
            blocks = []
            for k, seg in enumerate(segs, 1):
                a, b = seg["range"]
                trows = ""
                for row in seg["rows"]:
                    nm = "Target" if row["common"] == "target" else row["common"].capitalize()
                    ident = row.get("identity_to_target")
                    gaps = row.get("gaps_vs_target")
                    trows += (f'<tr><td class="lbl">{nm}</td><td>{row["coordinates"]}</td>'
                              f'<td>{row["length"]}</td>'
                              f'<td>{"&mdash;" if gaps is None else gaps}</td>'
                              f'{cell(ident) if ident is not None else "<td>&mdash;</td>"}</tr>')
                title = (f'{reg["label"]} {k} of {len(segs)} &middot; {ref} {a}&ndash;{b} '
                         f'({b - a + 1} aa)') if len(segs) > 1 else \
                        (f'{reg["label"]} &middot; {ref} {a}&ndash;{b} ({b - a + 1} aa)')
                opened = " open" if (b - a + 1) <= 120 else ""
                blocks.append(
                    f'<details class="seg"{opened}><summary>{title}</summary>'
                    f'<table><tr><th>Species</th><th>Coordinates</th><th>Length</th>'
                    f'<th>Indels</th><th>Identity to target</th></tr>{trows}</table>'
                    f'{wrapped_alignment(seg["alignment"], ref, a)}</details>')
            parts.append(f'<h3>{reg["label"]}</h3>{"".join(blocks)}')
        regions_html = (
            '<h2 style="margin-top:22px">Membrane regions, residue by residue</h2>'
            f'<p class="sub">Each region as defined on the {ref} topology, aligned across species. '
            'Red = residue differs from the target. Long segments are collapsed: click to open. '
            'Identities here are measured on these reference-defined segments, so they can '
            'differ slightly from the table above, which uses each species&rsquo; own '
            'topology.</p>' + "".join(parts))

    notes_html = ""
    if notes:
        items = "".join(f"<li>{_html.escape(n)}</li>" for n in notes)
        notes_html = (f'<div class="notes"><b>Notes on this run</b> &mdash; read these before '
                      f'the numbers<ul>{items}</ul></div>')

    return f'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>AbScout &mdash; {gene} in {target["species"]}</title>
<style>
 body{{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;color:#222;max-width:820px;margin:28px auto;padding:0 18px;line-height:1.5}}
 h1{{font-size:22px;margin:0 0 2px}} .sub{{color:#6b7280;margin:0 0 16px;font-size:14px}}
 h2{{font-size:16px;margin:0 0 6px}}
 table{{border-collapse:collapse;width:100%;margin:16px 0}}
 th,td{{border:1px solid #e5e7eb;padding:9px 10px;text-align:center;font-size:14px}}
 th{{background:#f7f8fa;font-weight:600}} td.lbl,th:first-child{{text-align:left}}
 .acc{{display:block;font-weight:400;color:#9aa0a6;font-size:12px}}
 .legend span{{display:inline-block;padding:2px 8px;border-radius:4px;color:#fff;font-size:12px;margin:0 4px 4px 0}}
 .aln{{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:13px;line-height:1.7;overflow-x:auto;white-space:pre;background:#fafafa;border:1px solid #eee;border-radius:6px;padding:8px}}
 .aln .al{{display:inline-block;width:64px;color:#6b7280}}
 .aln .d{{background:#fde2e2;color:#b91c1c;font-weight:700}}
 .numbox{{background:#fffbeb;border:1px solid #f3d38a;border-radius:6px;padding:9px 12px;color:#4b5563}}
 code{{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:13px;color:#111}}
 h3{{font-size:15px;margin:18px 0 6px}}
 details.seg{{border:1px solid #e5e7eb;border-radius:6px;padding:6px 12px;margin:8px 0}}
 details.seg summary{{cursor:pointer;font-weight:600;font-size:14px;padding:2px 0}}
 .aln .pos{{color:#9aa0a6;font-size:11px;margin-top:6px}}
 .notes{{background:#f8fafc;border:1px solid #cbd5e1;border-radius:6px;padding:8px 14px;margin:14px 0;font-size:13px;color:#334155}}
 .notes ul{{margin:4px 0 2px 18px;padding:0}} .notes li{{margin:3px 0}}
 .caveat{{color:#6b7280;font-size:12px;border-top:1px solid #e5e7eb;padding-top:12px;margin-top:22px}}
</style></head><body>
<h1>AbScout &mdash; {gene} in <i>{target["species"]}</i></h1>
<p class="sub">Target protein {target["uniprot"]}</p>
{notes_html}
<div class="legend">Sequence identity to target:
 <span style="background:#1b9e57">&ge;85%</span>
 <span style="background:#e0972a">70&ndash;85%</span>
 <span style="background:#cf3b3b">&lt;70%</span>
 <span style="background:#eef0f2;color:#8a919a">n/a</span></div>
{main_table}
{ep_html}
{regions_html}
<p class="caveat">These are percent identities for your interpretation &mdash; not a guarantee
of binding (conformational epitopes, glycosylation, and single hot-spot residues can still break
binding). Generated by AbScout {__version__}{" (web version)" if IN_BROWSER else ""}.</p>
</body></html>'''


def rank_and_report(gene, target, targets, model_results, epitope_info=None, notes=None):
    order = {c: i for i, c in enumerate(MODEL_ORDER)}
    rows = sorted(model_results, key=lambda m: order.get(m["common"], 99))
    soluble = bool(model_results) and all(m.get("soluble") for m in model_results)

    L = [f"AbScout -- {gene} in {target['species']} (taxid {target['taxid']})",
         f"Target protein: {target['uniprot']}"]
    if len(targets) > 1:
        L.append(f"(note: {len(targets)} target paralogs found; analysed {target['uniprot']})")

    if soluble:
        L += ["", "No transmembrane segment detected -- this appears to be a soluble "
                  "(non-membrane) protein, so the extracellular / transmembrane / intracellular "
                  "split does not apply.", "", "Whole-protein sequence identity to target (%):", ""]
        L.append(f"  {'model':<10}{'overall':>10}")
        for m in rows:
            L.append(f"  {m['common']:<10}{_fmt(m['overall']):>10}")
    else:
        L += ["", "Sequence identity to target, by membrane region (%):", ""]
        L.append(f"  {'model':<10}{'extracellular':>15}{'transmembrane':>15}"
                 f"{'intracellular':>15}{'overall':>10}")
        for m in rows:
            r = m["regions"]
            L.append(f"  {m['common']:<10}{_fmt(r.get('extracellular')):>15}{_fmt(r.get('TM')):>15}"
                     f"{_fmt(r.get('intracellular')):>15}{_fmt(m['overall']):>10}")

    if epitope_info:
        er = epitope_info["range"]
        loc = epitope_info.get("located")
        loc_note = (f", auto-located from the pasted peptide at {loc['match_identity']:.0f}% match"
                    if loc else "")
        rq = epitope_info.get("requested_range")
        if rq:
            loc_note += f"; requested {rq[0]}-{rq[1]}, cut to the protein's length"
        L += ["", f"Antibody epitope region (defined on {epitope_info['ref']}, residues "
                  f"{er[0]}-{er[1]}{loc_note}) -- identity to target (%):", ""]
        L.append(f"  {'species':<10}{'coordinates':>14}{'len':>6}{'indels':>8}"
                 f"{'identity':>10}{'similarity':>12}{'radical':>9}")
        for row in epitope_info["rows"]:
            ident = _fmt(row["identity_to_target"]) if row["identity_to_target"] is not None else "-"
            if row.get("absent") and row["common"] != "target":
                ident = "absent"
            simil = (_fmt(row.get("similarity_to_target"))
                     if row.get("similarity_to_target") is not None else "-")
            gaps = row.get("gaps_vs_target")
            rad = row.get("radical_subs")
            L.append(f"  {row['common']:<10}{row['coordinates']:>14}{row['length']:>6}"
                     f"{('-' if gaps is None else gaps):>8}{ident:>10}{simil:>12}"
                     f"{('-' if rad is None else rad):>9}")
        L.append("  identity = exact matches; similarity = identical or conservative "
                 "(BLOSUM62 > 0); radical = substitutions with Grantham > 100")

        gone = [r["common"] for r in epitope_info["rows"]
                if r.get("absent") and r["common"] != "target"]
        if gone:
            L.append(f"\n  [note] no residue of {', '.join(gone)} aligns to the epitope region -- "
                     "the entry is probably truncated there, so no identity can be given "
                     "(this is not 0%).")
        indel_rows = [r for r in epitope_info["rows"]
                      if (r.get("gaps_vs_target") or 0) > 0 and not r.get("absent")]
        if indel_rows:
            who = ", ".join(f"{r['common']} ({r['gaps_vs_target']})" for r in indel_rows)
            L.append(f"\n  [note] insertion/deletion inside the epitope vs the target: {who}. "
                     "Identity is computed over columns present in both species, so indels are "
                     "excluded from it -- read the indel count alongside the percentage.")

        subs = epitope_info.get("substitutions") or []
        if subs:
            L += ["", f"  Substitutions, {epitope_info['ref']} -> target "
                      f"({len(subs)} over the epitope):", ""]
            L.append(f"    {'pos':>7}  {'change':<8}{'Grantham':>9}  class")
            for s in subs:
                g = "-" if s["grantham"] is None else f"{s['grantham']:.0f}"
                L.append(f"    {s['position']:>7}  {s['from']}->{s['to']:<5}{g:>9}  {s['class']}")
            rad = [s for s in subs if s["class"] == "radical"]
            if rad:
                where = ", ".join(f"{s['from']}{s['position']}{s['to']}" for s in rad)
                L.append(f"\n  [note] physicochemically radical substitution(s) in the epitope: "
                         f"{where}. These change charge, size or polarity substantially and are "
                         "the changes most likely to affect binding.")

        gl = epitope_info.get("glyco") or {}
        if gl.get("gained_in_target") or gl.get("lost_in_target"):
            parts = []
            if gl.get("gained_in_target"):
                parts.append("gained in target at " +
                             ", ".join(str(p) for p in gl["gained_in_target"]))
            if gl.get("lost_in_target"):
                parts.append("lost in target at " +
                             ", ".join(str(p) for p in gl["lost_in_target"]))
            L.append(f"\n  [note] N-glycosylation sequon (N-X-S/T) {'; '.join(parts)} "
                     f"({epitope_info['ref']} numbering). A glycan can block antibody access, "
                     "so this can change reactivity independently of identity.")

        nb = epitope_info.get("numbering")
        if nb:
            ref = epitope_info["ref"]
            L.append(f"\n  Numbering check -- range {er[0]}-{er[1]} read as {ref} UniProt "
                     f"numbering (from Met1) of {nb['reference_uniprot']}, "
                     f"{nb['reference_length']} aa.")
            L.append(f"    residues at {er[0]}-{er[1]}: {nb['residues_at_range']}")
            if nb.get("mature_equivalent"):
                why = " + ".join(nb.get("processing") or []) or "N-terminal processing"
                a, b = nb["mature_equivalent"]
                if nb.get("mature_fits"):
                    L.append(f"    {ref}'s mature chain starts at residue {nb['mature_start']} "
                             f"({why}). If the datasheet numbered the MATURE protein, the same "
                             f"range is UniProt {a}-{b}:")
                    L.append(f"    residues at {a}-{b}: {nb['mature_residues']}")
                    L.append(f"    Compare both with the datasheet; if the second matches, re-run "
                             f"with --epitope-coor {a}-{b}. --epitope-seq avoids the ambiguity "
                             "entirely.")
                else:
                    L.append(f"    {ref}'s mature chain starts at residue {nb['mature_start']} "
                             f"({why}). Read as mature-protein numbering the range would be "
                             f"UniProt {a}-{b}, which runs past the end of the protein, so it "
                             "can only be UniProt numbering.")
            else:
                L.append(f"    {ref} has no signal peptide or propeptide annotated, so "
                         "mature-protein and UniProt numbering coincide here.")

        aln = epitope_info.get("alignment")
        if aln:
            L += ["", "  Epitope alignment (^ = position varies across species):", ""]
            names = ([n for n in MODEL_ORDER if n in aln]
                     + [n for n in aln if n not in MODEL_ORDER and n != "target"]
                     + (["target"] if "target" in aln else []))
            lblw = max(len(n) for n in names)
            alen = len(next(iter(aln.values())))
            for start in range(0, alen, 60):
                for n in names:
                    L.append(f"  {n:<{lblw}}  {aln[n][start:start + 60]}")
                seqs = [aln[n][start:start + 60] for n in names]
                marker = "".join("^" if len({s[i] for s in seqs if s[i] != '-'}) > 1 else " "
                                 for i in range(len(seqs[0])))
                L += [f"  {'':<{lblw}}  {marker}", ""]

    low_cov = [m for m in model_results
               if m.get("coverage_pct") is not None and m["coverage_pct"] < 80]
    if low_cov:
        names = ", ".join(f"{m['common']} {m['coverage_pct']:.0f}%" for m in low_cov)
        L.append(f"\n  [note] low alignment coverage ({names}) -- identity there describes only "
                 "the aligned part of the protein; check for a partial sequence or a "
                 "different isoform.")

    if notes:
        L += ["", "  Notes on this run (read these before the numbers):"]
        for n in notes:
            L.append(f"    - {n}")

    L.append("\n  Percent identities are provided for your interpretation -- "
             "not a guarantee of binding.")

    result = {"abscout_version": __version__,
              "generated": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
              "gene": gene, "target": target, "n_target_paralogs": len(targets),
              "models": model_results}
    if epitope_info:
        result["epitope"] = epitope_info
    return "\n".join(L), result


BATCH_COLUMNS = ["case", "gene", "target_taxon", "model", "uniprot", "length",
                 "extracellular", "transmembrane", "intracellular", "overall",
                 "coverage_pct", "topology_source", "epitope_range", "epitope_ref",
                 "epitope_residues", "epitope_alt_numbering",
                 "epitope_coordinates", "epitope_length", "epitope_indels",
                 "epitope_identity", "epitope_similarity", "epitope_radical_subs",
                 "epitope_glyco_gained", "epitope_glyco_lost",
                 "overall_similarity", "radical_subs", "note", "error"]


def run_batch(path, out_csv, timeout):
    """
    Run many cases from a CSV and write one tidy results table -- one row per
    (case, model species). Input columns: gene, target_taxon are required;
    case, epitope_coor, epitope_seq, epitope_ref are optional.
    """
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        in_cols = list(reader.fieldnames or [])
        cases = [r for r in reader if r.get("gene") and r.get("target_taxon")]
    if not cases:
        raise SystemExit(f"No usable rows in {path} (need at least: gene,target_taxon)")

    # any extra input columns (e.g. known_outcome, evidence) ride along so the results CSV
    # can be analysed directly against ground truth without a join
    known = {"case", "gene", "target_taxon", "epitope_coor", "epitope_seq", "epitope_ref"}
    # also exclude names AbScout writes itself, or the output would get a duplicate header
    # and the user's value would be clobbered by the computed one
    reserved = known | set(BATCH_COLUMNS)
    passthrough = [c for c in in_cols if c not in reserved]
    clashes = [c for c in in_cols if c in set(BATCH_COLUMNS) - known]
    if clashes:
        sys.stderr.write(f"[warn] input column(s) {', '.join(clashes)} clash with AbScout's own "
                         "output columns and were not carried through; rename them to keep them.\n")

    ro.ENSEMBL_TIMEOUT = timeout
    n_ok = n_err = 0

    def clip(s, n=60):
        return s if len(s) <= n else s[:n] + "..."

    # Rows are written as each case finishes, so an interrupted run keeps what it has done.
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=BATCH_COLUMNS + passthrough, extrasaction="ignore")
        w.writeheader()
        for i, c in enumerate(cases, 1):
            name = (c.get("case") or f"{c['gene']}_{c['target_taxon']}").strip()
            sys.stderr.write(f"\n[batch {i}/{len(cases)}] {name}\n")
            base = {"case": name, "gene": c["gene"].strip(),
                    "target_taxon": c["target_taxon"].strip()}
            base.update({k: c.get(k, "") for k in passthrough})
            rows = []
            try:
                epitope = None
                coor = (c.get("epitope_coor") or "").strip()
                if coor:
                    a, b = coor.split("-")
                    epitope = (int(a), int(b))
                with tempfile.TemporaryDirectory() as work:
                    resolved, target, targets, model_results, epitope_info = analyze(
                        c["gene"].strip(), int(c["target_taxon"]), work, epitope,
                        (c.get("epitope_ref") or "human").strip(),
                        (c.get("epitope_seq") or "").strip() or None)
                case_note = " | ".join((resolved or {}).get("notes", []))
                if not model_results:
                    rows.append({**base, "note": case_note,
                                 "error": "no model orthologs compared"})
                ei = epitope_info or {}
                ep_by = {r["common"]: r for r in ei.get("rows", [])}
                nb = ei.get("numbering") or {}
                ep_res = clip((ei.get("alignment") or {}).get(ei.get("ref"), ""))
                alt = ""
                if nb.get("mature_fits"):
                    alt = (f"{nb['mature_equivalent'][0]}-{nb['mature_equivalent'][1]}:"
                           f"{nb['mature_residues']}")
                for m in model_results:
                    r = m["regions"]
                    ep = ep_by.get(m["common"], {})
                    rows.append({**base, "model": m["common"], "uniprot": m["uniprot"],
                                 "length": m.get("length"),
                                 "extracellular": r.get("extracellular"),
                                 "transmembrane": r.get("TM"),
                                 "intracellular": r.get("intracellular"),
                                 "overall": m["overall"],
                                 "coverage_pct": m.get("coverage_pct"),
                                 "topology_source": m.get("topology_source"),
                                 "epitope_range": ("-".join(map(str, ei["range"]))
                                                   if ei else ""),
                                 "epitope_ref": ei.get("ref", ""),
                                 "epitope_residues": ep_res,
                                 "epitope_alt_numbering": alt,
                                 "epitope_coordinates": ep.get("coordinates", ""),
                                 "epitope_length": ep.get("length", ""),
                                 "epitope_indels": ep.get("gaps_vs_target", ""),
                                 "epitope_identity": ep.get("identity_to_target", ""),
                                 "epitope_similarity": ep.get("similarity_to_target", ""),
                                 "epitope_radical_subs": ep.get("radical_subs", ""),
                                 "epitope_glyco_gained": (
                                     len(ei.get("glyco", {}).get("gained_in_target", []))
                                     if ei else ""),
                                 "epitope_glyco_lost": (
                                     len(ei.get("glyco", {}).get("lost_in_target", []))
                                     if ei else ""),
                                 "overall_similarity": m.get("overall_similarity"),
                                 "radical_subs": m.get("radical_subs"),
                                 "note": case_note,
                                 "error": ""})
            except (Exception, SystemExit) as exc:      # one bad case must not sink the batch
                sys.stderr.write(f"  [error] {name}: {exc}\n")
                rows = [{**base, "error": str(exc)[:300]}]
            for row in rows:
                w.writerow(row)
                if row.get("error"):
                    n_err += 1
                else:
                    n_ok += 1
            fh.flush()

    sys.stderr.write(f"\n[wrote] {out_csv}  ({n_ok} result rows, "
                     f"{n_err} errors, {len(cases)} cases)\n")


def parse_epitope_coor(text):
    """'300-350' -> (300, 350)."""
    try:
        a, b = str(text).strip().split("-")
        return int(a), int(b)
    except ValueError:
        raise ValueError("epitope coordinates must look like START-END, e.g. 300-350") from None


def run_single(gene, target_taxon, epitope=None, epitope_ref="human", epitope_seq=None):
    """
    One case end to end. The command line and the web page both call this, so the two give
    the same numbers from the same code. Returns {"report": terminal text, "result": the
    JSON-ready dict, "html": the report page, "has_regions": bool}; raises AbScoutError or
    ValueError when the case cannot run.
    """
    with tempfile.TemporaryDirectory() as work:
        resolved, target, targets, model_results, epitope_info = analyze(
            gene, target_taxon, work, epitope, epitope_ref, epitope_seq)
    if not model_results:
        raise AbScoutError("No model orthologs could be compared (no sequences/topology).")
    run_notes = resolved.get("notes", [])
    report, result = rank_and_report(gene, target, targets, model_results, epitope_info,
                                     notes=run_notes)
    result["notes"] = run_notes
    result["abscout_version"] = __version__ + (" (web)" if IN_BROWSER else "")
    result["generated"] = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    regions = resolved.get("region_alignments")
    if regions:
        result["region_alignments"] = regions
    html = render_html(gene, target, model_results, epitope_info, notes=run_notes,
                       regions=regions)
    return {"report": report, "result": result, "html": html, "has_regions": bool(regions)}


def main():
    ap = argparse.ArgumentParser(
        description="AbScout: cross-species antibody cross-reactivity by membrane region.")
    ap.add_argument("--gene", help="human gene symbol, e.g. EGFR")
    ap.add_argument("--target-taxon", type=int, help="NCBI taxid of the target species")
    ap.add_argument("--json", help="optional path to write the full result as JSON")
    ap.add_argument("--html", help="optional path to write a self-contained HTML report")
    ap.add_argument("--epitope-coor", "--epitope", dest="epitope_coor",
                    help="antibody immunogen residue range on the reference, e.g. 300-350")
    ap.add_argument("--epitope-seq", help="antibody immunogen PEPTIDE sequence from the datasheet; "
                                          "AbScout locates it in the reference automatically")
    ap.add_argument("--epitope-ref", default="human",
                    help="which model the epitope refers to (default: human)")
    ap.add_argument("--timeout", type=int, default=120,
                    help="Ensembl request timeout in seconds; raise for large gene families "
                         "such as TJP1 or OCLN (default: 120, escalates on retry)")
    ap.add_argument("--batch", help="CSV of cases (columns: gene,target_taxon[,case,"
                                    "epitope_coor,epitope_seq,epitope_ref]) to run in one go")
    ap.add_argument("--batch-out", default="abscout_batch.csv",
                    help="where to write batch results (default: abscout_batch.csv)")
    ap.add_argument("--no-cache", action="store_true",
                    help="ignore the on-disk cache of sequences and topology")
    ap.add_argument("--version", action="version", version=f"AbScout {__version__}")
    args = ap.parse_args()

    global CACHE_ENABLED
    CACHE_ENABLED = not args.no_cache

    if args.batch:
        run_batch(args.batch, args.batch_out, args.timeout)
        return
    if not args.gene or not args.target_taxon:
        ap.error("--gene and --target-taxon are required (or use --batch)")

    ro.ENSEMBL_TIMEOUT = args.timeout

    epitope = None
    if args.epitope_coor:
        try:
            epitope = parse_epitope_coor(args.epitope_coor)
        except ValueError:
            ap.error("--epitope-coor must look like START-END, e.g. 300-350")

    try:
        out = run_single(args.gene, args.target_taxon, epitope, args.epitope_ref,
                         args.epitope_seq)
    except (ValueError, AbScoutError) as exc:
        # a problem with what was asked for (bad coordinates, unavailable reference species)
        # or with what could be found (no ortholog): say so plainly, no traceback
        raise SystemExit(f"abscout: {exc}")
    print(out["report"])
    if out["has_regions"] and not args.html:
        sys.stderr.write("\n  Residue-by-residue alignments of the extracellular, transmembrane and "
                         "intracellular regions are in the HTML report (--html) and JSON (--json).\n")
    if args.json:
        Path(args.json).write_text(json.dumps(out["result"], indent=2))
        sys.stderr.write(f"\n[wrote] {args.json}\n")
    if args.html:
        Path(args.html).write_text(out["html"])
        sys.stderr.write(f"[wrote] {args.html}\n")


if __name__ == "__main__":
    main()
