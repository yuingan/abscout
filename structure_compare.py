#!/usr/bin/env python3
"""
structure_compare.py
====================
Per-region structural comparison of a protein between a TARGET species and a
MODEL (reference) species, for the cross-species antibody tool.

It fetches both AlphaFold DB structures, aligns them, and reports structural
similarity GLOBALLY and -- the whole point -- SEPARATELY for the extracellular,
transmembrane, and intracellular regions. It is pLDDT-aware, so you can tell
whether a regional difference is real biology or just low-confidence AlphaFold.

  * extracellular numbers  -> relevant for surface / flow / IHC antibodies
  * intracellular numbers  -> relevant for western / lysate antibodies

Region definitions come from the REFERENCE protein's TMbed topology (so the
regions are anchored on the well-characterised model-species protein).

Dependencies : numpy, biopython, requests   (already in your env)
Optional     : tmtools  -> true TM-align alignment + official TM-score.
               Without it, falls back to a sequence-alignment correspondence
               (fine for orthologs; weaker for divergent pairs).

USAGE
  # fetch both from AlphaFold DB and compare per region using a TMbed .pred for the reference
  python structure_compare.py --reference P00533 --target <ACC> --ref-topology egfr.pred

  # global-only (no per-region) if no topology is given
  python structure_compare.py --reference P00533 --target <ACC>

NOTE: the AlphaFold-DB fetch could not be exercised where this was written (no
network to ebi.ac.uk); test it on P00533 first. The alignment/scoring logic is
unit-tested on synthetic structures (see test_structure_compare.py).
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import requests
from Bio.PDB import MMCIFParser, PDBParser
from Bio.Align import PairwiseAligner, substitution_matrices

try:
    from tmtools import tm_align as _tm_align
    HAVE_TMTOOLS = True
except Exception:
    HAVE_TMTOOLS = False

AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

REGION_ORDER = ["extracellular", "TM", "intracellular", "signal", "periplasmic", "other"]


def code_to_region(c):
    """Map a TMbed (or DeepTMHMM) per-residue label to a region name."""
    if c in "BbHhM":
        return "TM"
    if c == "S":
        return "signal"
    if c in "iI":
        return "intracellular"
    if c in "oO":
        return "extracellular"
    if c in "pP":
        return "periplasmic"
    return "other"


# ----------------------------------------------------------------------------------
# Fetch + load structures
# ----------------------------------------------------------------------------------
def fetch_afdb(acc, outdir, prefer="cif"):
    """Download the AlphaFold DB model for a UniProt accession. Returns local path."""
    Path(outdir).mkdir(parents=True, exist_ok=True)
    meta = requests.get(AFDB_API.format(acc=acc), timeout=30)
    meta.raise_for_status()
    js = meta.json()
    if not js:
        raise ValueError(f"No AlphaFold model found for {acc}")
    entry = js[0]
    url = entry.get("cifUrl" if prefer == "cif" else "pdbUrl") or \
        entry.get("cifUrl") or entry.get("pdbUrl")
    if not url:
        raise ValueError(f"No model file URL in AlphaFold record for {acc}")
    ext = ".cif" if url.endswith(".cif") else ".pdb"
    path = Path(outdir) / f"AF-{acc}{ext}"
    path.write_bytes(requests.get(url, timeout=120).content)
    return str(path)


def load_structure(path):
    """Return {ca: Nx3 float, resnums: [int], plddt: Nx float, seq: str} for chain 1."""
    parser = MMCIFParser(QUIET=True) if str(path).endswith(".cif") else PDBParser(QUIET=True)
    structure = parser.get_structure("x", path)
    model = next(structure.get_models())
    chain = next(model.get_chains())
    ca, resnums, plddt, seq = [], [], [], []
    for res in chain:
        if res.id[0] != " " or "CA" not in res:   # skip hetero/water/missing-CA
            continue
        atom = res["CA"]
        ca.append(atom.get_coord())
        resnums.append(res.id[1])
        plddt.append(atom.get_bfactor())           # B-factor column = pLDDT in AFDB
        seq.append(THREE_TO_ONE.get(res.resname, "X"))
    if not ca:
        raise ValueError(f"No CA atoms parsed from {path}")
    return {"ca": np.asarray(ca, float), "resnums": resnums,
            "plddt": np.asarray(plddt, float), "seq": "".join(seq)}


# ----------------------------------------------------------------------------------
# Superposition + metrics (numpy only)
# ----------------------------------------------------------------------------------
def superpose(P, Q):
    """Optimally superpose P onto Q (both Nx3). Return (rmsd, per_residue_distances)."""
    if len(P) < 3:
        return float("nan"), np.full(len(P), np.nan)
    Pc = P - P.mean(0)
    Qc = Q - Q.mean(0)
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    Pr = Pc @ R.T
    dist = np.sqrt(((Pr - Qc) ** 2).sum(1))
    return float(np.sqrt((dist ** 2).mean())), dist


def tm_score_from_dist(dist, l_ref):
    """Approximate TM-score from aligned-pair distances, normalised by reference length."""
    dist = dist[~np.isnan(dist)]
    if l_ref <= 0 or len(dist) == 0:
        return float("nan")
    d0 = max(1.24 * (max(l_ref - 15, 1)) ** (1.0 / 3.0) - 1.8, 0.5)
    return float((1.0 / (1.0 + (dist / d0) ** 2)).sum() / l_ref)


def lddt_per_residue(A, B, pairs, inclusion_radius=15.0, thresholds=(0.5, 1.0, 2.0, 4.0)):
    """
    CA-lDDT per reference residue (superposition-FREE; robust to domain/hinge motion,
    which is why it is more reliable than RMSD for flexible multi-domain proteins).
    For each reference residue it measures how well its local Cα-Cα distances are
    preserved in the target. Returns {reference_residue_index: lddt_fraction_0_to_1}.
    """
    if len(pairs) < 2:
        return {}
    a2b = dict(pairs)
    aidx = np.array(sorted(a2b))
    bidx = np.array([a2b[a] for a in aidx])
    PA = A["ca"][aidx]
    PB = B["ca"][bidx]
    DA = np.linalg.norm(PA[:, None, :] - PA[None, :, :], axis=2)
    DB = np.linalg.norm(PB[:, None, :] - PB[None, :, :], axis=2)
    mask = (DA < inclusion_radius) & (DA > 0.0)          # in-radius reference pairs, exclude self
    diff = np.abs(DA - DB)
    preserved = np.zeros_like(DA)
    for t in thresholds:
        preserved += (diff < t)
    preserved /= len(thresholds)
    num = (preserved * mask).sum(axis=1)
    den = mask.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        per = np.where(den > 0, num / den, np.nan)
    return {int(a): float(v) for a, v in zip(aidx, per)}


# ----------------------------------------------------------------------------------
# Residue correspondence
# ----------------------------------------------------------------------------------
def global_tm_score(A, B):
    """
    Global TM-score from a true TM-align structural superposition (normalised by the
    reference length). Used ONLY as a supplementary global number -- never for residue
    correspondence, because a single rigid superposition misaligns the domains of flexible
    multi-domain proteins and would corrupt per-region sequence identity.
    """
    try:
        res = _tm_align(A["ca"], B["ca"], A["seq"], B["seq"])
        return float(res.tm_norm_chain1)
    except Exception:
        return float("nan")


_BLOSUM62 = substitution_matrices.load("BLOSUM62")


# ----------------------------------------------------------------------------------
# Conservation-weighted substitution scoring
#
# Percent identity treats every mismatch alike, but an antibody does not: Leu->Ile is
# nearly invisible to it while Arg->Glu reverses the charge. Two epitopes at the same
# identity can therefore differ in whether binding survives. Two complementary measures:
#
#   BLOSUM62  -- evolutionary exchangeability. A score >= 0 is the standard definition of
#                a "conservative" substitution (BLAST's "positives").
#   Grantham  -- physicochemical distance, computed from Grantham (1974) Science 185:862
#                via his own formula over composition, polarity and molecular volume, so
#                the values are derived rather than transcribed. Range ~5-215.
# ----------------------------------------------------------------------------------
# amino acid -> (composition c, polarity p, molecular volume v)   [Grantham 1974, Table 1]
_GRANTHAM_PROPS = {
    "A": (0.00, 8.1, 31.0),   "R": (0.65, 10.5, 124.0), "N": (1.33, 11.6, 56.0),
    "D": (1.38, 13.0, 54.0),  "C": (2.75, 5.5, 55.0),   "Q": (0.89, 10.5, 85.0),
    "E": (0.92, 12.3, 83.0),  "G": (0.74, 9.0, 3.0),    "H": (0.58, 10.4, 96.0),
    "I": (0.00, 5.2, 111.0),  "L": (0.00, 4.9, 111.0),  "K": (0.33, 11.3, 119.0),
    "M": (0.00, 5.7, 105.0),  "F": (0.00, 5.2, 132.0),  "P": (0.39, 8.0, 32.5),
    "S": (1.42, 9.2, 32.0),   "T": (0.71, 8.6, 61.0),   "W": (0.13, 5.4, 170.0),
    "Y": (0.20, 6.2, 136.0),  "V": (0.00, 5.9, 84.0),
}
_G_ALPHA, _G_BETA, _G_GAMMA, _G_RHO = 1.833, 0.1018, 0.000399, 50.723

# Grantham's own conventional bands (his "mean chemical difference" classes)
GRANTHAM_CONSERVATIVE = 50      # <= 50  : conservative
GRANTHAM_RADICAL = 100          # >  100 : radical (moderately radical upward)


def grantham_distance(a, b):
    """Grantham physicochemical distance between two residues (nan if either is unknown)."""
    pa, pb = _GRANTHAM_PROPS.get(a), _GRANTHAM_PROPS.get(b)
    if pa is None or pb is None:
        return float("nan")
    dc, dp, dv = pa[0] - pb[0], pa[1] - pb[1], pa[2] - pb[2]
    return _G_RHO * math.sqrt(_G_ALPHA * dc * dc + _G_BETA * dp * dp + _G_GAMMA * dv * dv)


def substitution_class(d):
    """Band a Grantham distance: conservative / moderate / radical."""
    if d != d:                                  # nan
        return "unknown"
    if d <= GRANTHAM_CONSERVATIVE:
        return "conservative"
    if d <= GRANTHAM_RADICAL:
        return "moderate"
    return "radical"


def blosum62_score(a, b):
    """BLOSUM62 score for a residue pair, or None for unknown residues."""
    try:
        return float(_BLOSUM62[a, b])
    except (KeyError, IndexError):
        return None


def is_conservative(a, b):
    """
    True if the substitution is conservative by BLOSUM62, using BLAST's "positives"
    convention of a strictly positive score. Score 0 is deliberately excluded: Arg->Glu
    scores 0 yet reverses the charge, which is not conservative for a surface contact.
    """
    s = blosum62_score(a, b)
    return bool(s is not None and s > 0)


def score_substitutions(pairs_ab):
    """
    Summarise a list of (ref_residue, target_residue) pairs.
    Returns identity/similarity percentages plus how disruptive the differences are:
      similarity_pct      -- identical OR conservative by BLOSUM62 ("positives")
      radical_subs        -- substitutions with Grantham > 100
      mean_grantham_subs  -- mean Grantham over SUBSTITUTED positions only, i.e. how
                             disruptive the changes that do exist are (nan if none)
    """
    n = same = similar = radical = 0
    gsum = gn = 0.0
    for a, b in pairs_ab:
        if a == "-" or b == "-":
            continue
        n += 1
        if a == b:
            same += 1
            similar += 1
            continue
        if is_conservative(a, b):
            similar += 1
        d = grantham_distance(a, b)
        if d == d:
            gsum += d
            gn += 1
            if d > GRANTHAM_RADICAL:
                radical += 1
    if not n:
        return {"compared": 0, "seq_identity_pct": None, "similarity_pct": None,
                "radical_subs": 0, "mean_grantham_subs": None}
    return {"compared": n,
            "seq_identity_pct": round(100 * same / n, 1),
            "similarity_pct": round(100 * similar / n, 1),
            "radical_subs": radical,
            "mean_grantham_subs": (round(gsum / gn, 1) if gn else None)}


def correspondence_seq(A, B):
    """Sequence-alignment-based aligned index pairs (fallback). Returns (pairs, info)."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = _BLOSUM62
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aln = aligner.align(A["seq"], B["seq"])[0]
    pairs = []
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        for k in range(a1 - a0):
            pairs.append((a0 + k, b0 + k))
    return pairs, {"engine": "biopython/sequence"}


# ----------------------------------------------------------------------------------
# Core comparison
# ----------------------------------------------------------------------------------
def compare(A, B, ref_topology=None):
    """
    Compare two loaded structures. A = reference (model species), B = target.

    Sequence identity AND lDDT are computed from a SEQUENCE (homology) alignment. This is
    correct for orthologs and, crucially, robust for multi-domain proteins: a single rigid
    structural superposition (TM-align) misaligns domains that move relative to one another,
    which would corrupt per-region identity. TM-align is therefore used only for an optional
    supplementary global TM-score, never for the per-region numbers.
    """
    pairs, _ = correspondence_seq(A, B)
    if not pairs:
        raise ValueError("No aligned residues between the two structures.")

    iA = [p[0] for p in pairs]
    iB = [p[1] for p in pairs]
    g_rmsd, g_dist = superpose(B["ca"][iB], A["ca"][iA])
    g_ident = sum(A["seq"][a] == B["seq"][b] for a, b in pairs) / len(pairs)
    per_res_lddt = lddt_per_residue(A, B, pairs)
    g_lddt = 100 * float(np.nanmean(list(per_res_lddt.values()))) if per_res_lddt else float("nan")
    g_tm = global_tm_score(A, B) if HAVE_TMTOOLS else tm_score_from_dist(g_dist, len(A["seq"]))
    engine = "sequence alignment" + (" (+TM-align global TM-score)" if HAVE_TMTOOLS else "")

    result = {
        "engine": engine,
        "global": {
            "aligned_residues": len(pairs),
            "seq_identity_pct": round(100 * g_ident, 1),
            "lddt": round(g_lddt, 1),                 # superposition-free; robust per-region metric
            "tm_score": round(float(g_tm), 3),        # global, rigid; under-reads multi-domain
            "rmsd_A": round(g_rmsd, 2),               # global, rigid; unreliable for flexible domains
        },
        "regions": {},
    }

    if ref_topology:
        regions = [code_to_region(c) for c in ref_topology]
        dot_frac = ref_topology.count(".") / max(len(ref_topology), 1)
        if dot_frac > 0.3:
            sys.stderr.write(
                f"[warn] {dot_frac*100:.0f}% of the topology is '.', so TMbed did NOT label "
                "inside vs outside (default format). Re-run TMbed with `--out-format 4` to get "
                "the extracellular/intracellular split.\n")
        if len(regions) != len(A["seq"]):
            sys.stderr.write(
                f"[warn] topology length ({len(regions)}) != reference length "
                f"({len(A['seq'])}); aligning by index up to the shorter one. "
                "Make sure the TMbed .pred is for THIS reference sequence.\n")
        buckets = {}
        for a, b in pairs:
            if a < len(regions):
                buckets.setdefault(regions[a], []).append((a, b))
        for reg in sorted(buckets, key=lambda r: REGION_ORDER.index(r) if r in REGION_ORDER else 99):
            prs = buckets[reg]
            ia = [p[0] for p in prs]
            ib = [p[1] for p in prs]
            r_rmsd, r_dist = superpose(B["ca"][ib], A["ca"][ia])
            frac4 = float(np.nanmean(r_dist < 4.0)) if len(prs) >= 3 else float("nan")
            r_ident = sum(A["seq"][a] == B["seq"][b] for a, b in prs) / len(prs)
            lvals = [per_res_lddt[a] for a, _ in prs if per_res_lddt.get(a) == per_res_lddt.get(a)]
            r_lddt = 100 * float(np.mean(lvals)) if lvals else float("nan")
            result["regions"][reg] = {
                "aligned_residues": len(prs),
                "seq_identity_pct": round(100 * r_ident, 1),
                "lddt": (round(r_lddt, 1) if r_lddt == r_lddt else None),
                "rmsd_A": (round(r_rmsd, 2) if r_rmsd == r_rmsd else None),
                "frac_CA_within_4A": (round(frac4, 2) if frac4 == frac4 else None),
                "mean_plddt_ref": round(float(np.mean([A["plddt"][a] for a, _ in prs])), 1),
                "mean_plddt_target": round(float(np.mean([B["plddt"][b] for _, b in prs])), 1),
            }
    return result


# ----------------------------------------------------------------------------------
# Sequence-only comparison (no structures needed) -- AbScout's core output
# ----------------------------------------------------------------------------------
def compare_sequences(seqA, seqB, ref_topology=None, epitope=None):
    """
    Per-region % sequence identity between reference A and target B, from a sequence
    alignment alone -- NO structures required. This is AbScout's core, and it works for
    any UniProt accession (unlike AlphaFold DB, which doesn't cover every recent entry).
    Use compare() with structures if you additionally want lDDT/RMSD.

    epitope: optional (start, end) 1-based residue range on reference A (e.g. an antibody's
    immunogen region). If given, the result includes conservation of exactly that region in
    the target -- the direct "will this specific antibody cross-react" signal.
    """
    pairs, _ = correspondence_seq({"seq": seqA}, {"seq": seqB})
    if not pairs:
        raise ValueError("No aligned residues between the two sequences.")
    g_ident = sum(seqA[a] == seqB[b] for a, b in pairs) / len(pairs)
    # Identity is computed over ALIGNED positions only, so a sequence that aligns over just
    # part of its length can still score high identity. Coverage (aligned / reference length)
    # exposes that: low coverage means the identity describes only a fraction of the protein.
    coverage = len(pairs) / max(len(seqA), 1)
    g_score = score_substitutions([(seqA[a], seqB[b]) for a, b in pairs])
    result = {"engine": "sequence alignment",
              "global": {"aligned_residues": len(pairs),
                         "reference_residues": len(seqA),
                         "coverage_pct": round(100 * coverage, 1),
                         "seq_identity_pct": round(100 * g_ident, 1),
                         "similarity_pct": g_score["similarity_pct"],
                         "radical_subs": g_score["radical_subs"],
                         "mean_grantham_subs": g_score["mean_grantham_subs"]},
              "regions": {}}
    if ref_topology:
        regions = [code_to_region(c) for c in ref_topology]
        dot_frac = ref_topology.count(".") / max(len(ref_topology), 1)
        if dot_frac > 0.3:
            sys.stderr.write(f"[warn] {dot_frac*100:.0f}% of the topology is '.'; re-run TMbed "
                             "with `--out-format 4` for the extracellular/intracellular split.\n")
        if len(regions) != len(seqA):
            sys.stderr.write(f"[warn] topology length ({len(regions)}) != reference length "
                             f"({len(seqA)}); ensure the topology is for this reference.\n")
        buckets = {}
        for a, b in pairs:
            if a < len(regions):
                buckets.setdefault(regions[a], []).append((a, b))
        region_sizes = {}
        for r in regions:
            region_sizes[r] = region_sizes.get(r, 0) + 1
        for reg in sorted(buckets, key=lambda r: REGION_ORDER.index(r) if r in REGION_ORDER else 99):
            prs = buckets[reg]
            r_ident = sum(seqA[a] == seqB[b] for a, b in prs) / len(prs)
            r_size = region_sizes.get(reg, len(prs))
            r_score = score_substitutions([(seqA[a], seqB[b]) for a, b in prs])
            result["regions"][reg] = {"aligned_residues": len(prs),
                                      "reference_residues": r_size,
                                      "coverage_pct": round(100 * len(prs) / max(r_size, 1), 1),
                                      "seq_identity_pct": round(100 * r_ident, 1),
                                      "similarity_pct": r_score["similarity_pct"],
                                      "radical_subs": r_score["radical_subs"],
                                      "mean_grantham_subs": r_score["mean_grantham_subs"]}
    if epitope:
        es, ee = epitope                                   # 1-based inclusive, on reference A
        ep = [(a, b) for a, b in pairs if es - 1 <= a <= ee - 1]
        if ep:
            ep_ident = sum(seqA[a] == seqB[b] for a, b in ep) / len(ep)
            result["epitope"] = {"range": [es, ee], "reference_residues": ee - es + 1,
                                 "aligned_residues": len(ep),
                                 "seq_identity_pct": round(100 * ep_ident, 1)}
        else:
            result["epitope"] = {"range": [es, ee], "reference_residues": ee - es + 1,
                                 "aligned_residues": 0, "seq_identity_pct": None}
    return result


# ----------------------------------------------------------------------------------
# TMbed .pred parser
# ----------------------------------------------------------------------------------
def parse_tmbed_pred(path):
    """Parse a TMbed 3-line .pred file -> {id: label_string}."""
    lines = [ln.rstrip("\n") for ln in open(path)]
    out, i = {}, 0
    while i < len(lines):
        if lines[i].startswith(">"):
            hid = lines[i][1:].split()[0]
            out[hid] = lines[i + 2] if i + 2 < len(lines) else ""
            i += 3
        else:
            i += 1
    return out


# ----------------------------------------------------------------------------------
# Report + CLI
# ----------------------------------------------------------------------------------
def simple_report(reference, target, result):
    """Immunologist-facing summary: % sequence identity per membrane region + overall."""
    regs = result["regions"]
    L = [f"{target}  vs  {reference}   -- sequence identity by membrane region", ""]
    rows = [("extracellular", "extracellular (outside)"),
            ("TM", "transmembrane (in membrane)"),
            ("intracellular", "intracellular (inside)")]
    shown = False
    for key, label in rows:
        if key in regs:
            L.append(f"    {label:<30}{regs[key]['seq_identity_pct']:6.1f} %")
            shown = True
    if shown:
        L.append(f"    {'-' * 36}")
    L.append(f"    {'overall (whole protein)':<30}{result['global']['seq_identity_pct']:6.1f} %")
    L.append("")
    ec = regs.get("extracellular")
    ic = regs.get("intracellular")
    if ec:
        pct = ec["seq_identity_pct"]
        verdict = ("strong candidate" if pct >= 85 else
                   "borderline -- verify" if pct >= 70 else "unlikely to cross-react")
        L.append(f"  Surface / flow / IHC antibodies depend on the extracellular region: "
                 f"{pct:.0f}% identity -> {verdict}.")
    if ic:
        L.append(f"  Western / lysate antibodies depend on the intracellular region: "
                 f"{ic['seq_identity_pct']:.0f}% identity.")
    L.append("  (~85% identity is the usual antibody cross-reactivity rule of thumb -- calibrate it "
             "on your Phase 0 set. Structural detail is in the JSON; use --detailed to print it.)")
    return "\n".join(L)


def format_report(reference, target, result):
    L = []
    L.append(f"Structural comparison: target {target}  vs  reference {reference}")
    L.append(f"  alignment engine: {result['engine']}")
    g = result["global"]
    L.append(f"  GLOBAL  identity {g['seq_identity_pct']:5.1f}%   lDDT {g['lddt']:5.1f}   "
             f"TM-score {g['tm_score']:.3f}   RMSD {g['rmsd_A']:.2f} A   "
             f"({g['aligned_residues']} aligned)")
    if result["regions"]:
        L.append("  BY REGION:")
        L.append(f"    {'region':<15}{'ident%':>8}{'lDDT':>7}{'RMSD(A)':>9}"
                 f"{'pLDDT_ref':>11}{'pLDDT_tgt':>11}{'n':>6}")
        for reg, r in result["regions"].items():
            lddt = f"{r['lddt']:.1f}" if r["lddt"] is not None else "  -"
            rmsd = f"{r['rmsd_A']:.2f}" if r["rmsd_A"] is not None else "  -"
            L.append(f"    {reg:<15}{r['seq_identity_pct']:>8.1f}{lddt:>7}{rmsd:>9}"
                     f"{r['mean_plddt_ref']:>11.1f}{r['mean_plddt_target']:>11.1f}"
                     f"{r['aligned_residues']:>6}")
        ec = result["regions"].get("extracellular")
        if ec:
            lline = f"lDDT {ec['lddt']:.0f}" if ec["lddt"] is not None else ""
            L.append(f"\n  Antibody read: extracellular identity {ec['seq_identity_pct']:.1f}%, "
                     f"{lline} (mean pLDDT ref {ec['mean_plddt_ref']:.0f} / tgt {ec['mean_plddt_target']:.0f}). "
                     "This region drives surface/flow/IHC antibody cross-reactivity.")
        L.append("\n  lDDT (0-100, higher = more similar) is superposition-free and is the metric to "
                 "trust per region. RMSD is shown for reference but is inflated for large flexible/"
                 "multi-domain regions even when the fold is conserved. Low pLDDT in a region means "
                 "treat its structural numbers as uncertain rather than as real divergence.")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description="Per-region structural comparison via AlphaFold DB.")
    ap.add_argument("--reference", required=True, help="UniProt accession of the model-species protein")
    ap.add_argument("--target", required=True, help="UniProt accession of the target-species protein")
    ap.add_argument("--ref-topology", help="TMbed .pred file for the REFERENCE (enables per-region)")
    ap.add_argument("--outdir", default="afdb_cache", help="where to cache downloaded structures")
    ap.add_argument("--json", help="optional path to write the full result (all metrics) as JSON")
    ap.add_argument("--detailed", action="store_true",
                    help="print full structural metrics (lDDT/RMSD/pLDDT) instead of the simple % view")
    args = ap.parse_args()

    if HAVE_TMTOOLS:
        sys.stderr.write("[info] identity & lDDT come from the sequence alignment; "
                         "tmtools adds an exact global TM-score.\n")
    else:
        sys.stderr.write("[info] identity & lDDT come from the sequence alignment (correct for "
                         "orthologs). Global TM-score is approximated; `pip install tmtools` for "
                         "the exact value.\n")

    ref_path = fetch_afdb(args.reference, args.outdir)
    tgt_path = fetch_afdb(args.target, args.outdir)
    A = load_structure(ref_path)
    B = load_structure(tgt_path)

    topo = None
    if args.ref_topology:
        preds = parse_tmbed_pred(args.ref_topology)
        topo = next(iter(preds.values())) if preds else None

    result = compare(A, B, ref_topology=topo)
    print(format_report(args.reference, args.target, result) if args.detailed
          else simple_report(args.reference, args.target, result))
    if args.json:
        Path(args.json).write_text(json.dumps(result, indent=2))
        sys.stderr.write(f"\n[wrote] {args.json}\n")


if __name__ == "__main__":
    main()
