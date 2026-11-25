#!/usr/bin/env python3
"""
test_structure_compare.py — validates the alignment + per-region scoring logic of
structure_compare.py on synthetic CA-only structures with KNOWN per-region
differences. No network or real structures required.

Run:  python test_structure_compare.py
"""
import importlib.util
import tempfile
from pathlib import Path

import numpy as np

spec = importlib.util.spec_from_file_location("sc", "structure_compare.py")
sc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc)

AA = "ARNDCQEGHILKMFPSTWYV"


def helix_coords(n, seed=0):
    """A non-degenerate 3D backbone (alpha-helix-like), CA-CA ~3.8 A."""
    rng = np.random.default_rng(seed)
    t = np.arange(n)
    # alpha helix: ~100 deg/residue, 1.5 A rise; radius ~2.3 A -> CA-CA ~3.8 A
    coords = np.stack([2.3 * np.cos(np.deg2rad(100 * t)),
                       2.3 * np.sin(np.deg2rad(100 * t)),
                       1.5 * t], axis=1)
    return coords.astype(float)


def write_ca_pdb(path, coords, seq, plddt):
    aa3 = {v: k for k, v in sc.THREE_TO_ONE.items()}
    with open(path, "w") as fh:
        for i, (xyz, aa, b) in enumerate(zip(coords, seq, plddt), start=1):
            fh.write("ATOM  {:>5d}  CA  {:<3s} A{:>4d}    "
                     "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}           C\n".format(
                         i, aa3.get(aa, "GLY"), i, xyz[0], xyz[1], xyz[2], 1.00, b))
        fh.write("END\n")


def random_rigid_transform(coords, seed=1):
    """Apply a known rotation + translation (superposition must remove it)."""
    rng = np.random.default_rng(seed)
    a, b, c = rng.uniform(0, 2 * np.pi, 3)
    Rx = np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
    Ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
    Rz = np.array([[np.cos(c), -np.sin(c), 0], [np.sin(c), np.cos(c), 0], [0, 0, 1]])
    R = Rz @ Ry @ Rx
    return coords @ R.T + rng.uniform(-50, 50, 3)


def main():
    n = 40
    # regions: 1-15 extracellular ('o'), 16-25 TM ('H'), 26-40 intracellular ('i')
    topology = "o" * 15 + "H" * 10 + "i" * 15
    assert len(topology) == n

    ref_coords = helix_coords(n, seed=0)
    seq = "".join(AA[i % 20] for i in range(n))

    # target: same backbone + region-specific gaussian noise (EC small, TM large, IC medium),
    # then a rigid transform that superposition must undo.
    rng = np.random.default_rng(42)
    noise_std = np.where(np.arange(n) < 15, 0.2,
                         np.where(np.arange(n) < 25, 1.5, 0.5))
    tgt_coords = ref_coords + rng.normal(0, 1, (n, 3)) * noise_std[:, None]
    tgt_coords = random_rigid_transform(tgt_coords, seed=7)

    # pLDDT: reference high everywhere; target lower in the extracellular region
    ref_plddt = np.full(n, 90.0)
    tgt_plddt = np.where(np.arange(n) < 15, 60.0, 88.0)

    with tempfile.TemporaryDirectory() as d:
        ref_pdb = Path(d) / "ref.pdb"
        tgt_pdb = Path(d) / "tgt.pdb"
        write_ca_pdb(ref_pdb, ref_coords, seq, ref_plddt)
        write_ca_pdb(tgt_pdb, tgt_coords, seq, tgt_plddt)   # identical sequence -> 100% identity

        A = sc.load_structure(str(ref_pdb))
        B = sc.load_structure(str(tgt_pdb))
        assert A["seq"] == seq, "sequence round-trip through PDB failed"
        assert np.allclose(A["plddt"], 90.0), "pLDDT (B-factor) not read correctly"

        res = sc.compare(A, B, ref_topology=topology)

    print("engine:", res["engine"])
    print("global:", res["global"])
    for reg, r in res["regions"].items():
        print(f"  {reg:<15} {r}")

    g = res["global"]
    assert g["seq_identity_pct"] == 100.0, "identical sequences should be 100% identity"
    assert g["aligned_residues"] == n, "all residues should align for identical sequences"

    reg = res["regions"]
    for r in ("extracellular", "TM", "intracellular"):
        assert r in reg, f"missing region {r}"
    # KEY checks: per-region RMSD reflects the injected noise (TM most perturbed)
    assert reg["extracellular"]["rmsd_A"] < reg["TM"]["rmsd_A"], \
        "extracellular should superpose tighter than TM"
    assert reg["intracellular"]["rmsd_A"] < reg["TM"]["rmsd_A"], \
        "intracellular should superpose tighter than TM"
    # rigid transform was fully removed (EC RMSD should be small, ~0.2*sqrt(3))
    assert reg["extracellular"]["rmsd_A"] < 0.6, \
        f"rigid transform not removed? EC RMSD={reg['extracellular']['rmsd_A']}"
    # pLDDT correctly attributed per region
    assert abs(reg["extracellular"]["mean_plddt_target"] - 60.0) < 1e-6, "target EC pLDDT wrong"
    assert abs(reg["extracellular"]["mean_plddt_ref"] - 90.0) < 1e-6, "ref EC pLDDT wrong"

    # lDDT present and behaves: less-perturbed regions score higher; identical -> ~100
    assert reg["extracellular"]["lddt"] > reg["TM"]["lddt"], "EC should have higher lDDT than TM"
    res_self = sc.compare(A, A, ref_topology=topology)
    assert res_self["global"]["lddt"] > 99.0, \
        f"self-comparison lDDT should be ~100, got {res_self['global']['lddt']}"
    print("lDDT OK: EC %.1f > TM %.1f ; self-comparison %.1f" % (
        reg["extracellular"]["lddt"], reg["TM"]["lddt"], res_self["global"]["lddt"]))

    # --- also exercise the sequence-alignment fallback path ---
    saved = sc.HAVE_TMTOOLS
    try:
        sc.HAVE_TMTOOLS = False
        res2 = sc.compare(A, B, ref_topology=topology)
        assert res2["engine"].startswith("sequence"), "fallback should use sequence alignment"
        assert res2["global"]["seq_identity_pct"] == 100.0
        assert res2["regions"]["extracellular"]["rmsd_A"] < res2["regions"]["TM"]["rmsd_A"]
        print("fallback engine:", res2["engine"], "-> OK")
    finally:
        sc.HAVE_TMTOOLS = saved

    # --- parse_tmbed_pred ---
    with tempfile.NamedTemporaryFile("w", suffix=".pred", delete=False) as f:
        f.write(">sp|P00533|EGFR\n" + seq + "\n" + topology + "\n")
        pred_path = f.name
    parsed = sc.parse_tmbed_pred(pred_path)
    assert list(parsed.values())[0] == topology, "TMbed .pred parse mismatch"
    print("parse_tmbed_pred -> OK")

    # --- regression test for the multi-domain misalignment bug ---
    # Identity must come from sequence homology, NOT a rigid structural superposition:
    # rotating one domain in 3D must NOT change the reported % identity.
    m = 40
    mref = helix_coords(m, seed=3)
    mseq = "".join(AA[(i * 7) % 20] for i in range(m))
    mtgt_seq = list(mseq)
    for i in (3, 11, 27, 35):                        # 4 mutations -> 90% identity
        mtgt_seq[i] = AA[(AA.index(mseq[i]) + 5) % 20]
    mtgt_seq = "".join(mtgt_seq)
    mtgt = mref.copy()
    d2 = mtgt[20:]                                    # rigidly rotate the second "domain" 120 deg
    c2 = d2.mean(0)
    th = np.deg2rad(120)
    Rz = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    mtgt[20:] = (d2 - c2) @ Rz.T + c2
    with tempfile.TemporaryDirectory() as d:
        rp, tp = Path(d) / "mr.pdb", Path(d) / "mt.pdb"
        write_ca_pdb(rp, mref, mseq, np.full(m, 90.0))
        write_ca_pdb(tp, mtgt, mtgt_seq, np.full(m, 90.0))
        Am = sc.load_structure(str(rp))
        Bm = sc.load_structure(str(tp))
        rm = sc.compare(Am, Bm, ref_topology="o" * 20 + "i" * 20)
    assert 88.0 <= rm["global"]["seq_identity_pct"] <= 92.0, (
        "identity must reflect homology (~90%) regardless of domain motion, got "
        f"{rm['global']['seq_identity_pct']}% -- a structural correspondence would corrupt this")
    print("multi-domain identity OK: %.1f%% with a domain rotated 120deg (homology-based, unaffected)"
          % rm["global"]["seq_identity_pct"])

    # --- structure-free per-region identity (compare_sequences) ---
    order = "ACDEFGHIKLMNPQRSTVWY"
    sa = order                                   # first 10 = EC, last 10 (MNPQRSTVWY) = IC
    sb = order[:10] + "A" * 10                   # IC replaced with 'A' (absent from sa's IC) -> 0% there
    topo_s = "o" * 10 + "i" * 10
    rs = sc.compare_sequences(sa, sb, ref_topology=topo_s)
    assert rs["regions"]["extracellular"]["seq_identity_pct"] == 100.0, rs["regions"]
    assert rs["regions"]["intracellular"]["seq_identity_pct"] == 0.0, rs["regions"]
    assert rs["global"]["seq_identity_pct"] == 50.0, rs["global"]
    print("compare_sequences OK: EC 100%, IC 0%, global 50%% (no structures used)")

    # epitope-region conservation (antibody immunogen range on the reference)
    assert sc.compare_sequences(sa, sb, epitope=(1, 10))["epitope"]["seq_identity_pct"] == 100.0
    assert sc.compare_sequences(sa, sb, epitope=(11, 20))["epitope"]["seq_identity_pct"] == 0.0
    print("epitope conservation OK: residues 1-10 -> 100%, 11-20 -> 0%")

    print("\nALL STRUCTURE-COMPARE TESTS PASSED")




def test_conservation_scoring():
    """Grantham distances must reproduce published anchors; similarity must split
    conservative from radical substitutions."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("sc2", "structure_compare.py")
    sc2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(sc2)

    # Grantham (1974) published values, within rounding
    for a, b, expect in [("L", "I", 5), ("C", "W", 215), ("S", "R", 110),
                         ("G", "W", 184), ("R", "W", 101), ("D", "E", 45)]:
        d = sc2.grantham_distance(a, b)
        assert abs(d - expect) < 1.5, f"Grantham {a}->{b} = {d}, expected ~{expect}"
    assert sc2.substitution_class(sc2.grantham_distance("L", "I")) == "conservative"
    assert sc2.substitution_class(sc2.grantham_distance("C", "W")) == "radical"
    assert sc2.grantham_distance("X", "L") != sc2.grantham_distance("X", "L")   # nan

    # BLAST "positives" convention: strictly > 0, so Arg->Glu (score 0) is NOT conservative
    assert sc2.is_conservative("L", "I") and sc2.is_conservative("D", "E")
    assert not sc2.is_conservative("R", "E"), "charge reversal must not count as conservative"
    assert not sc2.is_conservative("C", "W")

    # identity is blind to severity; similarity is not
    cons = [("L", "I"), ("L", "I"), ("L", "I"), ("L", "L")]      # 25% identity, conservative
    rad = [("C", "W"), ("C", "W"), ("C", "W"), ("L", "L")]       # 25% identity, radical
    s_cons, s_rad = sc2.score_substitutions(cons), sc2.score_substitutions(rad)
    assert s_cons["seq_identity_pct"] == s_rad["seq_identity_pct"] == 25.0
    assert s_cons["similarity_pct"] == 100.0 and s_rad["similarity_pct"] == 25.0
    assert s_cons["radical_subs"] == 0 and s_rad["radical_subs"] == 3
    assert s_cons["mean_grantham_subs"] < 10 and s_rad["mean_grantham_subs"] > 200
    assert sc2.score_substitutions([("-", "L"), ("L", "-")])["compared"] == 0
    print("conservation scoring OK: Grantham anchors reproduce; identity 25%% for both, "
          "similarity 100%% vs 25%%")


def _protein(n, seed):
    import random
    r = random.Random(seed)
    return "".join(r.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(n))


def test_unmatched_stretches_and_whole_identity():
    """A target with a 20-residue N-terminal extension, a 6-residue insertion inside the
    extracellular region, and 10 residues deleted across the TM/intracellular boundary.
    Each stretch must be reported with its position and length, the insertion placed in
    the extracellular region, the terminal extension left unassigned, and the deletion split
    at the region boundary."""
    A = _protein(207, 7)
    topo = "S" * 22 + "o" * 104 + "H" * 26 + "i" * 55
    ext, ins = _protein(20, 8), "GGWGGW"
    B = ext + A[:100] + ins + A[100:150] + A[160:]
    res = sc.compare_sequences(A, B, ref_topology=topo)
    t_only = res["unmatched"]["target_only"]
    assert [(r["start"], r["end"], r["length"], r["where"], r["region"]) for r in t_only] == [
        (1, 20, 20, "N-terminal", None), (121, 126, 6, "internal", "extracellular")], t_only
    r_only = res["unmatched"]["reference_only"]
    assert [(r["start"], r["end"], r["region"]) for r in r_only] == [
        (151, 152, "TM"), (153, 160, "intracellular")], r_only
    g = res["global"]
    assert g["seq_identity_pct"] == 100.0, g                 # every aligned residue is identical
    assert g["target_only_residues"] == 26 and g["reference_only_residues"] == 10, g
    assert g["seq_identity_whole_pct"] == round(100 * 197 / (197 + 36), 1), g
    ec = res["regions"]["extracellular"]
    assert ec["seq_identity_pct"] == 100.0 and ec["target_only_residues"] == 6, ec
    assert ec["seq_identity_whole_pct"] == round(100 * 104 / 110, 1), ec
    tm = res["regions"]["TM"]
    assert tm["reference_only_residues"] == 2 and tm["seq_identity_whole_pct"] == round(
        100 * 24 / 26, 1), tm
    assert [r["placement"] for r in t_only] == ["terminal", "inside"], t_only
    assert all(r["placement"] == "inside" for r in r_only), r_only
    # a stretch the target lacks inside a part of the reference that is in no region (a GPI
    # propeptide, '.'), and the same comparison with no topology at all
    gpi = "S" * 19 + "o" * 111 + "." * 77
    res2 = sc.compare_sequences(A, A[:130], ref_topology=gpi)
    assert [r["placement"] for r in res2["unmatched"]["reference_only"]] == ["outside"], res2
    res3 = sc.compare_sequences(A, B)
    assert {r["placement"] for r in res3["unmatched"]["target_only"]} == {"no topology"}, res3
    print("unmatched OK: N-terminal 1-20 (20 aa, unassigned), 121-126 (6 aa, extracellular); "
          "deletion split TM 151-152 / intracellular 153-160; whole identity %.1f%%"
          % g["seq_identity_whole_pct"])


def test_project_topology():
    A = _protein(207, 7)
    topo = "S" * 22 + "o" * 104 + "H" * 26 + "i" * 55
    B = _protein(20, 8) + A[:100] + "GGWGGW" + A[100:]
    t = sc.project_topology(A, topo, B)
    assert len(t) == len(B)
    assert t[:20] == "." * 20, "a terminal extension has no side of the membrane to inherit"
    assert t[20:42] == "S" * 22 and t[120:126] == "oooooo", t[115:130]
    assert t.count("H") == 26 and t.endswith("i" * 55), t
    # an insertion exactly at a region boundary stays unassigned
    C = A[:126] + "PPPP" + A[126:]
    tc = sc.project_topology(A, topo, C)
    assert tc[126:130] == "....", tc[120:135]
    print("project_topology OK: extension '.', insertion inside a region inherits it, "
          "insertion at a boundary stays '.'")


def _diverge(seq, every, seed):
    """seq with one residue in `every` replaced by a different one (a diverged ortholog)."""
    import random
    r = random.Random(seed)
    out = list(seq)
    for i in range(0, len(seq), every):
        out[i] = r.choice([c for c in "ACDEFGHIKLMNPQRSTVWY" if c != seq[i]])
    return "".join(out)


def test_long_insertion_close_to_an_end():
    """Dog NCAM1 case: 254 residues more than human inside the cytoplasmic tail, 40 residues
    before the C-terminus. With free ends alone the aligner stopped early instead of opening
    the long gap and paired the last 40 human residues with the start of the insertion (5 of
    40 identical; their true partners give 39), so intracellular identity read 71.6% where it
    is 98.3%. The residues beyond the insertion must pair with their real partners, at either
    end of the protein, and the insertion must come out as one internal stretch."""
    A = _protein(400, 3)
    ins = _protein(270, 4)
    topo = "o" * 250 + "H" * 21 + "i" * 129
    for label, B, lo in (("identical", A[:365] + ins + A[365:], 100.0),
                         ("diverged", _diverge(A[:365], 7, 5) + ins + _diverge(A[365:], 7, 6),
                          80.0)):
        tail_true = [(i, i + 270) for i in range(365, 400)]
        pairs, info = sc.correspondence_seq({"seq": A}, {"seq": B})
        assert info["ends"] == "end to end", (label, info)
        # (a mismatched residue right at the insertion may sit on either side of the gap)
        kept = len(set(tail_true) & set(pairs))
        assert kept >= 33, (label, kept, "residues beyond the insertion mis-paired")
        # the free-ends alignment on its own is the failure being guarded against
        free = sc._align(A, B, free_ends=True)
        assert not set(tail_true) & set(free), (label, "free ends alone no longer fails?")
        res = sc.compare_sequences(A, B, ref_topology=topo)
        ic = res["regions"]["intracellular"]
        assert ic["seq_identity_pct"] >= lo and ic["aligned_residues"] == 129, (label, ic)
        assert res["global"]["alignment_ends"] == "end to end", res["global"]
        t_only = res["unmatched"]["target_only"]
        assert [(r["length"], r["where"], r["region"], r["placement"]) for r in t_only] == [
            (270, "internal", "intracellular", "inside")], (label, t_only)
        assert abs(t_only[0]["start"] - 366) <= 2, t_only
    # the same thing at the other end: the insertion 35 residues after the N-terminus
    B = A[:35] + ins + A[35:]
    pairs, info = sc.correspondence_seq({"seq": A}, {"seq": B})
    assert info["ends"] == "end to end", info
    assert sum(A[i] == B[j] for i, j in pairs) == 400 and len(pairs) == 400, len(pairs)
    print("long insertion near an end OK: residues beyond a 270-residue insertion keep their "
          "partners (C- and N-terminal side), one internal stretch, alignment end to end")


def test_end_extension_stays_one_stretch():
    """The other half of the rule. An extension at one end is not an insertion: it must stay
    one unmatched stretch and the rest must align as if it were not there. Ordinary pairs,
    where both candidate alignments agree, are reported as free-ends."""
    A = _protein(207, 7)
    for B, where, at in ((_protein(95, 9) + A, "N-terminal", (1, 95)),
                         (A + _protein(95, 9), "C-terminal", (208, 302))):
        pairs, info = sc.correspondence_seq({"seq": A}, {"seq": B})
        assert info["ends"] == "free", info
        res = sc.compare_sequences(A, B)
        t_only = res["unmatched"]["target_only"]
        assert [(r["start"], r["end"], r["where"]) for r in t_only] == [(at[0], at[1], where)], t_only
        assert res["global"]["seq_identity_pct"] == 100.0 and res["global"]["aligned_residues"] == 207
    same, info = sc.correspondence_seq({"seq": A}, {"seq": _diverge(A, 5, 2)})
    assert info["ends"] == "free" and len(same) == 207, info
    # the yardstick itself: pairs score BLOSUM62, a long internal gap costs the cap, ends are free
    s = sc._common_score("WWWW", "WWAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAWW", [(0, 0), (1, 1), (2, 42), (3, 43)])
    assert s == 4 * 11 - sc.LONG_GAP_CAP, s
    assert sc._common_score("WWWW", "AAAAAAAAAAWWWW", [(0, 10), (1, 11), (2, 12), (3, 13)]) == 44.0
    assert sc._common_score("WWAWW", "WWWW", [(0, 0), (1, 1), (3, 2), (4, 3)]) == 44.0 - sc.GAP_OPEN
    print("end extension OK: 95 residues at either end stay one stretch, free-ends alignment "
          "kept; yardstick charges a 40-residue internal gap %d, an overhang nothing"
          % sc.LONG_GAP_CAP)


# Real sequences, as Ensembl Compara holds them (November 2026): mouse THY1 and the naked
# mole-rat THY1 translation, whose first 95 residues have no counterpart in other species.
THY1_MOUSE = ("MNPAISVALLLSVLQVSRGQKVTSLTACLVNQNLRLDCRHENNTKDNSIQHEFSLTREKRKHVLSGTLGIPEHTYRSRVT"
              "LSNQPYIKVLTLANFTTKDEGDYFCELQVSGANPMSSNKSISVYRDKLVKCGGISLLVQNTSWMLLLLLSLSLLQALDFI"
              "SL")
THY1_NMR = ("CSLGCPTLIPGSQSLPPQARPAATGLPEVPLLQERATFLHPGSTATSSRAPRHPRPLPSASHWLSPGPSPSPSPPPLVKT"
            "AGSVQGAATGGGSIASRGSCSGGGRGPVLQVSRGQKVTSLTACLVNQTLRLDCHHENTTTLPVQYEFSLTREKKKHVLFG"
            "TMGVPEHAYRSRTNLTSNYNLKVLYLANFTNKDEGTYTCELRISGQQSIPNYKNVSVLRDKLVKCEGISLLAQNTSWLLL"
            "LLLSLSLLQATDFVSL")


def test_real_thy1_extension_is_not_split():
    """Naked mole-rat THY1 against mouse THY1: the case the free-ends alignment was introduced
    for. Aligned end to end, the 95-residue N-terminal extension is broken up, some of its
    residues paired with the start of the mouse protein. The two candidates differ here, so
    this is a real choice, and the free-ends alignment must be the one kept: the extension
    one stretch, 1-95, and the mouse signal peptide on the mole-rat's own."""
    free = sc._align(THY1_MOUSE, THY1_NMR, free_ends=True)
    whole = sc._align(THY1_MOUSE, THY1_NMR, free_ends=False)
    assert free != whole, "the two candidates agree: this test no longer exercises the choice"
    paired_in_extension = [j for _, j in whole if j < 95]
    assert paired_in_extension, "end to end no longer pulls extension residues onto mouse?"
    assert sc._common_score(THY1_MOUSE, THY1_NMR, free) > sc._common_score(
        THY1_MOUSE, THY1_NMR, whole)
    pairs, info = sc.correspondence_seq({"seq": THY1_MOUSE}, {"seq": THY1_NMR})
    assert info["ends"] == "free" and pairs == free, info
    res = sc.compare_sequences(THY1_MOUSE, THY1_NMR)
    assert [(r["start"], r["end"], r["where"]) for r in res["unmatched"]["target_only"]] == [
        (1, 95, "N-terminal")], res["unmatched"]["target_only"]
    assert res["global"]["aligned_residues"] == 161 and res["global"]["alignment_ends"] == "free"
    print("real THY1 OK: mole-rat extension 1-95 stays one stretch (end to end would pair %d of "
          "its residues with mouse); %d of 161 aligned residues identical"
          % (len(paired_in_extension), sum(THY1_MOUSE[i] == THY1_NMR[j] for i, j in pairs)))


if __name__ == "__main__":
    test_unmatched_stretches_and_whole_identity()
    test_project_topology()
    test_long_insertion_close_to_an_end()
    test_end_extension_stays_one_stretch()
    test_real_thy1_extension_is_not_split()
    main()
    test_conservation_scoring()
    print("\nALL STRUCTURE-COMPARE TESTS PASSED")
