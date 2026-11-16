#!/usr/bin/env python3
"""Test AbScout's report + epitope logic (no network/TMbed). Run: python test_abscout.py"""
import importlib.util

spec = importlib.util.spec_from_file_location("ab", "abscout.py")
ab = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ab)

TARGET = {"species": "danio_rerio", "taxid": 7955, "uniprot": "A0A8M1PEX5"}
MODELS = [
    {"species": "homo_sapiens", "common": "human", "uniprot": "P00533",
     "regions": {"extracellular": 78.2, "TM": 65.0, "intracellular": 70.1}, "overall": 73.0},
    {"species": "mus_musculus", "common": "mouse", "uniprot": "Q01279",
     "regions": {"extracellular": 71.5, "TM": 60.0, "intracellular": 68.0}, "overall": 69.0},
    {"species": "rattus_norvegicus", "common": "rat", "uniprot": "A0A8I5ZMJ2",
     "regions": {"extracellular": 72.0, "TM": 61.0, "intracellular": 69.0}, "overall": 70.0},
    {"species": "oryctolagus_cuniculus", "common": "rabbit", "uniprot": "G1SZE4",
     "regions": {"extracellular": 69.0, "TM": 58.0, "intracellular": 66.0}, "overall": 67.0},
]


def test_region_table_no_recommendation():
    report, result = ab.rank_and_report("EGFR", TARGET, [TARGET], MODELS)
    assert "recommendation" not in result, "recommendation section should be removed"
    assert "which antibody to try" not in report.lower()
    assert "strong cross-reactivity" not in report.lower(), "verdicts should be removed"
    assert "extracellular" in report.lower() and "78.2" in report
    print("region table (neutral, no recommendation/verdicts) OK")


def test_handles_missing_region_and_paralogs():
    models = [
        {"species": "homo_sapiens", "common": "human", "uniprot": "P1",
         "regions": {"extracellular": None, "TM": 40.0, "intracellular": 55.0}, "overall": 50.0},
        {"species": "mus_musculus", "common": "mouse", "uniprot": "P2",
         "regions": {"extracellular": 88.0, "TM": 70.0, "intracellular": 80.0}, "overall": 82.0},
    ]
    two_targets = [TARGET, {"species": "danio_rerio", "taxid": 7955, "uniprot": "OTHER"}]
    report, result = ab.rank_and_report("GENEX", TARGET, two_targets, models)
    assert result["n_target_paralogs"] == 2
    assert "2 target paralogs" in report
    # missing extracellular renders as a dash, doesn't crash
    assert "-" in report
    print("missing-region + paralog handling OK")


def test_uniprot_topology_parsing():
    feats = [
        {"type": "Signal", "location": {"start": {"value": 1}, "end": {"value": 24}}},
        {"type": "Topological domain", "description": "Extracellular",
         "location": {"start": {"value": 25}, "end": {"value": 645}}},
        {"type": "Transmembrane", "description": "Helical",
         "location": {"start": {"value": 646}, "end": {"value": 668}}},
        {"type": "Topological domain", "description": "Cytoplasmic",
         "location": {"start": {"value": 669}, "end": {"value": 1210}}},
    ]
    ts = ab._topology_string_from_features(feats, 1210)
    assert ts[0] == "S" and ts[24] == "o" and ts[645] == "H" and ts[668] == "i"
    soluble = [{"type": "Topological domain", "description": "Cytoplasmic",
                "location": {"start": {"value": 1}, "end": {"value": 100}}}]
    assert ab._topology_string_from_features(soluble, 100) is None
    tm_only = [{"type": "Transmembrane",
                "location": {"start": {"value": 50}, "end": {"value": 72}}}]
    assert ab._topology_string_from_features(tm_only, 200) is None, "TM-only must fall back to TMbed"
    print("uniprot topology parsing OK: full split mapped; TM-only and no-TM -> TMbed fallback")


def test_epitope_across_species():
    block = "ACDEFGHIKLMNPQRSTVWY"
    human = block * 3                                   # 60 aa
    target = human                                      # identical to human
    mouse = human[:20] + "A" * 20 + human[40:]          # residues 21-40 diverged
    ep = ab.epitope_across_species([("human", human), ("mouse", mouse)], target,
                                   "human", (21, 40))
    by = {r["common"]: r for r in ep["rows"]}
    assert by["human"]["identity_to_target"] == 100.0, by["human"]      # human == target
    assert by["human"]["length"] == 20 and by["human"]["coordinates"] == "21-40"
    assert by["mouse"]["identity_to_target"] < 50.0, by["mouse"]        # mouse diverged here
    assert "target" in by and by["target"]["identity_to_target"] is None
    aln = ep["alignment"]                                                # alignment for display
    assert len(aln["human"]) == 20 == len(aln["mouse"]) == len(aln["target"])
    assert aln["mouse"] == "A" * 20, aln["mouse"]                        # the diverged residues
    assert all(r["gaps_vs_target"] == 0 for r in ep["rows"] if r["common"] != "target")
    print("epitope_across_species OK: table + alignment, human 100%%, mouse <50%%")


def test_epitope_indel_is_reported():
    """A deletion in the target epitope must not hide behind a high identity."""
    block = "ACDEFGHIKLMNPQRSTVWY"
    human = block * 2                                   # 40 aa epitope
    target = block[:10] + block[15:] + block            # 5-residue deletion inside it
    ep = ab.epitope_across_species([("human", human)], target, "human", (1, 40))
    row = [r for r in ep["rows"] if r["common"] == "human"][0]
    assert row["identity_to_target"] == 100.0, row      # identity alone still looks perfect
    assert row["gaps_vs_target"] == 5, row              # ...but the indel is now surfaced
    assert row["epitope_coverage_pct"] < 90.0, row
    report, _ = ab.rank_and_report("GENEX", TARGET, [TARGET], MODELS,
                                   {"ref": "human", "range": [1, 40],
                                    "rows": ep["rows"], "alignment": ep["alignment"]})
    assert "insertion/deletion inside the epitope" in report, report
    print("epitope indel OK: 5-residue deletion reported as indels=5 despite 100%% identity")


def test_epitope_range_validation():
    """Out-of-range coordinates must raise, not silently yield an empty epitope table."""
    ref = "ACDEFGHIKLMNPQRSTVWY" * 2                  # 40 aa
    try:
        ab.epitope_across_species([("human", ref)], ref, "human", (1005, 1044))
    except ValueError as e:
        assert "outside" in str(e) and "signal peptide" in str(e), e
    else:
        raise AssertionError("coordinates past the sequence end should raise")
    # an end past the sequence is truncated rather than fatal
    ep = ab.epitope_across_species([("human", ref)], ref, "human", (35, 60))
    assert [r for r in ep["rows"] if r["common"] == "human"][0]["length"] == 6
    print("epitope range validation OK: out-of-range rejected, overhang truncated")


def test_peptide_match_floor():
    """An absent or mistyped peptide must be refused, not placed at its best-scoring window."""
    ref = "ACDEFGHIKLMNPQRSTVWY" * 5
    assert ab.locate_peptide("W" * 20, ref)[2] < 30.0   # the primitive still returns a window
    assert ab.detect_epitope_from_peptide([("human", ref)], "W" * 20) is None
    good = ref[40:60]                                    # a real slice of the reference
    hit = ab.detect_epitope_from_peptide([("human", ref)], good)
    assert hit and hit[0] == "human" and hit[3] == 100.0, hit
    print("peptide floor OK: 5%% match refused, exact match accepted")


def test_rank_report_epitope_table():
    epi = {"ref": "human", "range": [310, 360], "rows": [
        {"common": "human", "coordinates": "310-360", "length": 51, "identity_to_target": 94.0,
         "gaps_vs_target": 0, "epitope_coverage_pct": 100.0},
        {"common": "mouse", "coordinates": "308-358", "length": 51, "identity_to_target": 91.0,
         "gaps_vs_target": 0, "epitope_coverage_pct": 100.0},
        {"common": "target", "coordinates": "305-355", "length": 51, "identity_to_target": None,
         "gaps_vs_target": None, "epitope_coverage_pct": 100.0},
    ]}
    report, result = ab.rank_and_report("EGFR", TARGET, [TARGET], MODELS, epi)
    assert "epitope" in report.lower() and "310-360" in report and "94" in report
    assert "recommendation" not in result and result["epitope"]["ref"] == "human"
    print("epitope table render + no recommendation OK")


def test_locate_peptide():
    ref = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKL"
    assert ab.locate_peptide("GHIKLMNP", ref) == (6, 13, 100.0)          # exact
    assert ab.locate_peptide("g h i k l m n p", ref) == (6, 13, 100.0)   # whitespace tolerant
    loc = ab.locate_peptide("GHXKLMNP", ref)                             # one mismatch
    assert loc[0] == 6 and loc[1] == 13 and 80.0 <= loc[2] < 100.0, loc
    assert ab.locate_peptide("A" * 100, ref) is None                     # too long
    print("locate_peptide OK: exact, whitespace, fuzzy, too-long")


def test_detect_epitope_species():
    mouse = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKL"
    human = "ACDEFGHWWLMNPQRSTVWYACDEFGHIKL"          # differs from mouse at positions 8-9
    peptide = "GHIKLMNP"                                # exact slice of mouse; 2 mismatches in human
    best = ab.detect_epitope_from_peptide([("human", human), ("mouse", mouse)], peptide)
    assert best[0] == "mouse" and best[3] == 100.0, best  # picks mouse as the immunogen species
    print("detect_epitope_species OK: peptide -> mouse reference automatically")


def test_run_tmbed_rejects_short_sequence():
    # TMbed crashes (UnboundLocalError in make_batches) on sequences it cannot batch;
    # we must refuse before calling it.
    try:
        ab.run_tmbed("MKV", "/tmp", "frag")
    except ValueError as e:
        assert "too short" in str(e), e
        print("run_tmbed short-sequence guard OK:", e)
    else:
        raise AssertionError("run_tmbed should reject a 3 aa sequence")


def test_soluble_protein_report():
    # ZO-1 (TJP1) is cytoplasmic: no TM, so the region split must not be reported
    models = [
        {"species": "homo_sapiens", "common": "human", "uniprot": "Q07157", "soluble": True,
         "regions": {"extracellular": None, "TM": None, "intracellular": None}, "overall": 93.0},
        {"species": "mus_musculus", "common": "mouse", "uniprot": "P39447", "soluble": True,
         "regions": {"extracellular": None, "TM": None, "intracellular": None}, "overall": 92.3},
    ]
    report, result = ab.rank_and_report("TJP1", TARGET, [TARGET], models)
    assert "soluble" in report.lower() and "does not apply" in report.lower(), report
    assert "extracellular" not in report.split("Whole-protein")[1].lower()
    assert "93.0" in report
    print("soluble-protein report OK: no membrane split, whole-protein identity only")



def test_cache_roundtrip(tmpdir="/tmp/abscout_cache_test"):
    import shutil, pathlib
    shutil.rmtree(tmpdir, ignore_errors=True)
    old_dir, old_en = ab.CACHE_DIR, ab.CACHE_ENABLED
    ab.CACHE_DIR, ab.CACHE_ENABLED = pathlib.Path(tmpdir), True
    try:
        assert ab.cache_get("seq", "P00533") is None          # empty cache
        ab.cache_put("seq", "P00533", "MRPSGTAGAALLALLAALCPAS")
        assert ab.cache_get("seq", "P00533") == "MRPSGTAGAALLALLAALCPAS"
        # empty string distinguishes "checked, no annotation" from "not checked"
        ab.cache_put("topo_uniprot", "G1SZE4", "")
        assert ab.cache_get("topo_uniprot", "G1SZE4") == ""
        ab.CACHE_ENABLED = False                              # --no-cache must miss
        assert ab.cache_get("seq", "P00533") is None
    finally:
        ab.CACHE_DIR, ab.CACHE_ENABLED = old_dir, old_en
        shutil.rmtree(tmpdir, ignore_errors=True)
    print("cache roundtrip OK: hit, miss, empty-marker, and --no-cache honoured")


def test_coverage_reported():
    import structure_compare as sc
    full = "ACDEFGHIKLMNPQRSTVWY" * 3                   # 60 aa reference
    half = full[:30]                                     # target covers only half
    res = sc.compare_sequences(full, half)
    assert res["global"]["coverage_pct"] <= 55.0, res["global"]
    assert res["global"]["seq_identity_pct"] == 100.0    # perfect over what aligned
    print("coverage OK: identity 100%% but coverage %.0f%% -- exposes a partial sequence"
          % res["global"]["coverage_pct"])


def test_low_coverage_note():
    models = [{"species": "mus_musculus", "common": "mouse", "uniprot": "P1",
               "regions": {"extracellular": 90.0, "TM": 80.0, "intracellular": 88.0},
               "overall": 89.0, "coverage_pct": 55.0, "length": 400}]
    report, result = ab.rank_and_report("GENEX", TARGET, [TARGET], models)
    assert "low alignment coverage" in report, report
    assert result["abscout_version"] and result["generated"], "provenance missing"
    print("low-coverage note + provenance OK")


def test_batch_csv(tmp="/tmp/abscout_batch_test"):
    import shutil, pathlib, csv as _csv
    shutil.rmtree(tmp, ignore_errors=True)
    pathlib.Path(tmp).mkdir(parents=True)
    inp, outp = f"{tmp}/cases.csv", f"{tmp}/out.csv"
    pathlib.Path(inp).write_text(
        "case,gene,target_taxon\ncetuximab_mouse,EGFR,10090\nbad_row,,\n")
    real = ab.analyze                                     # stub: no network needed
    ab.analyze = lambda *a, **k: (
        None, {"species": "mus_musculus", "taxid": 10090, "uniprot": "Q01279"},
        [{"uniprot": "Q01279"}],
        [{"species": "homo_sapiens", "common": "human", "uniprot": "P00533",
          "regions": {"extracellular": 89.2, "TM": 73.9, "intracellular": 91.9},
          "overall": 89.7, "coverage_pct": 99.0, "length": 1210,
          "topology_source": "UniProt"}], None)
    try:
        ab.run_batch(inp, outp, 120)
    finally:
        ab.analyze = real
    rows = list(_csv.DictReader(open(outp)))
    assert len(rows) == 1, rows                           # blank row skipped
    assert rows[0]["case"] == "cetuximab_mouse"
    assert rows[0]["extracellular"] == "89.2" and rows[0]["model"] == "human"
    assert set(_csv.DictReader(open(outp)).fieldnames) == set(ab.BATCH_COLUMNS)
    shutil.rmtree(tmp, ignore_errors=True)
    print("batch mode OK: tidy CSV, one row per (case, model), bad rows skipped")




def test_glyco_sequon_detection():
    """N-X-[S/T] with X != P; overlapping sequons found; gains/losses reported in ref numbering."""
    assert ab.sequon_positions("AANASAA") == {2}            # N-A-S
    assert ab.sequon_positions("AANPSAA") == set()          # X = P is not a sequon
    assert ab.sequon_positions("AANATAA") == {2}            # N-A-T counts
    assert ab.sequon_positions("NNSSAA") == {0, 1}          # overlapping sequons
    print("sequon detection OK: N-X-S/T, proline excluded, overlaps found")


def test_glyco_gain_in_epitope():
    """A sequon gained in the target must be reported even at high identity."""
    block = "ACDEFGHIKLMNPQRSTVWY"
    ref = "AAAQATAAAA" + block          # no sequon at position 4 (Q-A-T)
    tgt = "AAANATAAAA" + block          # Q4N creates N-A-T -> sequon gained
    ep = ab.epitope_across_species([("human", ref)], tgt, "human", (1, 10))
    row = [r for r in ep["rows"] if r["common"] == "human"][0]
    assert row["identity_to_target"] == 90.0, row           # one substitution only
    assert ep["glyco"]["gained_in_target"] == [4], ep["glyco"]
    assert ep["glyco"]["lost_in_target"] == []
    report, _ = ab.rank_and_report("GENEX", TARGET, [TARGET], MODELS,
                                   {"ref": "human", "range": [1, 10], "rows": ep["rows"],
                                    "alignment": ep["alignment"],
                                    "substitutions": ep["substitutions"],
                                    "glyco": ep["glyco"]})
    assert "N-glycosylation sequon" in report and "gained in target at 4" in report, report
    print("glyco gain OK: Q4N sequon gained flagged despite 90%% identity")


def test_epitope_substitution_detail():
    """Per-residue substitutions carry Grantham distance and a severity class."""
    ref = "AAACAAAAAA"
    tgt = "AAAWAAAAAA"                  # C->W, the most radical pair (Grantham 215)
    ep = ab.epitope_across_species([("human", ref)], tgt, "human", (1, 10))
    subs = ep["substitutions"]
    assert len(subs) == 1, subs
    s = subs[0]
    assert s["position"] == 4 and s["from"] == "C" and s["to"] == "W"
    assert s["grantham"] >= 210 and s["class"] == "radical" and not s["conservative"]
    row = [r for r in ep["rows"] if r["common"] == "human"][0]
    assert row["radical_subs"] == 1 and row["similarity_to_target"] == 90.0
    print("substitution detail OK: C4W -> Grantham %d, radical" % s["grantham"])


def test_signal_peptide_numbering():
    """Leading 'S' in the topology gives the signal peptide length for the numbering note."""
    assert ab.signal_peptide_len("SSSSSoooooHHHHHiiiii") == 5
    assert ab.signal_peptide_len("ooooHHHHiiii") == 0
    assert ab.signal_peptide_len(None) == 0
    assert ab.signal_peptide_len("") == 0
    print("signal peptide length OK: 5 leading S, 0 when absent")




def _feat(t, s, e, d=""):
    return {"type": t, "description": d, "location": {"start": {"value": s}, "end": {"value": e}}}


def test_nterm_processing_real_feature_tables():
    """Mature-chain start from UniProt's own feature tables (as returned for these entries)."""
    cd3e = [_feat("Signal", 1, 22), _feat("Chain", 23, 207)]                      # P07766
    cdh1 = [_feat("Signal", 1, 22), _feat("Propeptide", 23, 154),                 # P12830
            _feat("Chain", 155, 882), _feat("Chain", 701, 882)]
    tfrc = [_feat("Chain", 1, 760), _feat("Topological domain", 1, 67, "Cytoplasmic")]   # P02786
    actb = [_feat("Chain", 1, 375), _feat("Initiator methionine", 1, 1, "Removed")]      # P60709
    n = ab._nterm_processing_from_features
    assert n(cd3e) == {"signal_len": 22, "mature_start": 23,
                       "parts": ["22-residue signal peptide"]}
    got = n(cdh1)
    assert got["mature_start"] == 155 and got["signal_len"] == 22, got     # not 23: propeptide
    assert "propeptide" in got["parts"]
    assert n(tfrc)["mature_start"] == 1 and n(tfrc)["parts"] == []         # type II, no signal
    assert n(actb)["mature_start"] == 2                                    # initiator Met removed
    print("N-terminal processing OK: CD3E 23, CDH1 155 (signal+propeptide), TFRC 1, ACTB 2")


def test_epitope_ref_must_be_a_compared_species():
    """An unavailable reference must raise, never silently anchor on another species."""
    seq = "ACDEFGHIKLMNPQRSTVWY" * 2
    try:
        ab.epitope_across_species([("human", seq)], seq, "rat", (1, 10))
    except ValueError as e:
        assert "'rat'" in str(e) and "human" in str(e), e
    else:
        raise AssertionError("a reference species that was not compared must raise")
    print("epitope reference validation OK: unavailable species rejected, available listed")


def test_range_used_is_reported():
    seq = "ACDEFGHIKLMNPQRSTVWY" * 2                    # 40 aa
    ep = ab.epitope_across_species([("human", seq)], seq, "human", (35, 60))
    assert ep["range_used"] == [35, 40], ep["range_used"]
    print("range_used OK: 35-60 on a 40 aa protein reported as 35-40")


def test_batch_survives_failure_and_carries_notes(tmp="/tmp/abscout_batch_test2"):
    """One failing case (even SystemExit) must not sink the batch; notes reach the CSV."""
    import shutil, pathlib, csv as _csv
    shutil.rmtree(tmp, ignore_errors=True)
    pathlib.Path(tmp).mkdir(parents=True)
    inp, outp = f"{tmp}/cases.csv", f"{tmp}/out.csv"
    pathlib.Path(inp).write_text(
        "case,gene,target_taxon\nfirst,EGFR,10090\nmissing,NOPE,10090\n"
        "exits,QUIT,10090\nlast,EGFR,10090\n")
    good = ({"notes": ["rat (X1) skipped: partial sequence, 100 aa vs target 500 aa"]},
            {"species": "mus_musculus", "taxid": 10090, "uniprot": "Q01279"},
            [{"uniprot": "Q01279"}],
            [{"species": "homo_sapiens", "common": "human", "uniprot": "P00533",
              "regions": {"extracellular": 89.2, "TM": 73.9, "intracellular": 91.9},
              "overall": 89.7, "coverage_pct": 99.0, "length": 1210,
              "topology_source": "UniProt"}], None)

    def fake(gene, *a, **k):
        if gene == "NOPE":
            raise ab.AbScoutError("No target-species ortholog")
        if gene == "QUIT":
            raise SystemExit("legacy exit")
        return good
    real = ab.analyze
    ab.analyze = fake
    try:
        ab.run_batch(inp, outp, 120)
    finally:
        ab.analyze = real
    rows = {r["case"]: r for r in _csv.DictReader(open(outp))}
    assert set(rows) == {"first", "missing", "exits", "last"}, set(rows)
    assert rows["missing"]["error"].startswith("No target-species ortholog")
    assert rows["exits"]["error"] == "legacy exit"
    assert rows["last"]["extracellular"] == "89.2" and not rows["last"]["error"]
    assert "partial sequence" in rows["first"]["note"]
    shutil.rmtree(tmp, ignore_errors=True)
    print("batch resilience OK: failures recorded per case, later cases still run, notes kept")




def test_duplicate_species_rejected():
    """Two orthologs of one species must not silently overwrite each other's epitope row."""
    block = "ACDEFGHIKLMNPQRSTVWY"
    same, different = block * 2, "W" * 20 + block
    try:
        ab.epitope_across_species([("human", same), ("rabbit", same), ("rabbit", different)],
                                  same, "human", (1, 20))
    except ValueError as e:
        assert "rabbit" in str(e) and "one ortholog per species" in str(e), e
    else:
        raise AssertionError("duplicate species names must be refused, not silently collapsed")
    print("duplicate-species guard OK: repeated species name refused")


def test_one_ortholog_chosen_per_species(monkey=None):
    """analyze() collapses several orthologs of a model species to the best one."""
    import resolve_orthologs as ro
    block = "ACDEFGHIKLMNPQRSTVWY"
    full, short = block * 10, block * 2                  # 200 aa vs 40 aa
    seqs = {"HUM": full, "RAB_GOOD": full, "RAB_SHORT": short, "TGT": full}

    def o(common, acc, status, otype, tax, is_t=False):
        return {"species": common, "common": common, "taxid": tax, "ensembl_gene": "G",
                "uniprot": acc, "uniprot_status": status, "note": "",
                "orthology_type": otype, "is_target": is_t}

    real_resolve, real_fetch, real_topo, real_tmbed, real_avail = (
        ro.resolve, ab.fetch_uniprot_seq, ab.uniprot_topology, ab.run_tmbed,
        ab.tmbed_available)
    ab.tmbed_available = lambda: True                 # this case is about TMbed being run once
    tmbed_calls = []
    ro.resolve = lambda g, t: {"gene": g, "target_taxid": t, "orthologs": [
        o("human", "HUM", "reviewed", "ortholog_one2one", 9606),
        o("rabbit", "RAB_SHORT", "unreviewed", "ortholog_one2many", 9986),
        o("rabbit", "RAB_GOOD", "reviewed", "ortholog_one2one", 9986),
        o("target", "TGT", "reviewed", "ortholog_one2one", 10181, True)]}
    ab.fetch_uniprot_seq = lambda acc: seqs[acc]
    ab.uniprot_topology = lambda acc, n: (None, "")
    def fake_tmbed(seq, work, tag, acc=None):
        tmbed_calls.append(acc)
        return "o" * 20 + "H" * 20 + "i" * (len(seq) - 40)
    ab.run_tmbed = fake_tmbed
    try:
        resolved, target, targets, models, _ = ab.analyze("GENEX", 10181, "/tmp")
    finally:
        (ro.resolve, ab.fetch_uniprot_seq, ab.uniprot_topology,
         ab.run_tmbed, ab.tmbed_available) = (real_resolve, real_fetch, real_topo, real_tmbed,
                                              real_avail)

    commons = [m["common"] for m in models]
    assert "RAB_GOOD" in tmbed_calls, tmbed_calls
    assert commons.count("rabbit") == 1, commons
    rab = [m for m in models if m["common"] == "rabbit"][0]
    assert rab["uniprot"] == "RAB_GOOD", rab          # reviewed, one2one, full length
    assert tmbed_calls.count("RAB_SHORT") == 0, "TMbed must not run on the discarded ortholog"
    note = " ".join(resolved["notes"])
    assert "rabbit: 2 orthologs" in note and "RAB_SHORT" in note, note
    print("one-per-species OK: RAB_GOOD kept over RAB_SHORT, discard noted, TMbed run once")




def test_short_reviewed_model_is_rescued_not_skipped():
    """A short reviewed model (human CDH1 496 aa vs 882) must be upgraded, not dropped --
    dropping it also took out the epitope reference and failed the whole case."""
    import resolve_orthologs as ro
    block = "ACDEFGHIKLMNPQRSTVWY"
    full, short = block * 44, block * 24                 # 880 aa vs 480 aa
    seqs = {"SHORT_REV": short, "FULL_REV": full, "TGT": full}
    upgrade_args = []

    def o(common, acc, status, tax, is_t=False):
        return {"species": common, "common": common, "taxid": tax, "ensembl_gene": "G",
                "uniprot": acc, "uniprot_status": status, "note": "",
                "orthology_type": "ortholog_one2one", "is_target": is_t}

    saved = (ro.resolve, ab.fetch_uniprot_seq, ab.uniprot_topology,
             ab.run_tmbed, ro.upgrade_partial)
    ro.resolve = lambda g, t: {"gene": g, "target_taxid": t, "orthologs": [
        o("human", "SHORT_REV", "reviewed", 9606),
        o("target", "TGT", "reviewed", 10181, True)]}
    ab.fetch_uniprot_seq = lambda acc: seqs[acc]
    ab.uniprot_topology = lambda acc, n: ("o" * 20 + "H" * 20 + "i" * (n - 40), "tm")
    ab.run_tmbed = lambda *a, **k: (_ for _ in ()).throw(AssertionError("TMbed should not run"))

    def fake_upgrade(acc, status, gene, taxid, ratio=1.3):
        upgrade_args.append((acc, status))
        return ("FULL_REV", "reviewed", "partial 480aa -> complete 880aa")
    ro.upgrade_partial = fake_upgrade
    try:
        resolved, target, targets, models, epi = ab.analyze(
            "CDH1", 10181, "/tmp", epitope=(1, 20), epitope_ref="human")
    finally:
        (ro.resolve, ab.fetch_uniprot_seq, ab.uniprot_topology,
         ab.run_tmbed, ro.upgrade_partial) = saved

    assert upgrade_args == [("SHORT_REV", "reviewed")], upgrade_args
    assert len(models) == 1 and models[0]["uniprot"] == "FULL_REV", models
    assert epi is not None and epi["ref"] == "human", "the epitope reference must survive"
    assert any("was partial" in n for n in resolved["notes"]), resolved["notes"]
    print("short-reviewed rescue OK: 480aa SHORT_REV -> 880aa FULL_REV, epitope still reported")


def test_full_length_model_costs_no_extra_lookup():
    """The rescue must not fire for a normal model -- it would add UniProt calls per run."""
    import resolve_orthologs as ro
    block = "ACDEFGHIKLMNPQRSTVWY"
    full = block * 44
    saved = (ro.resolve, ab.fetch_uniprot_seq, ab.uniprot_topology, ro.upgrade_partial)
    ro.resolve = lambda g, t: {"gene": g, "target_taxid": t, "orthologs": [
        {"species": "human", "common": "human", "taxid": 9606, "ensembl_gene": "G",
         "uniprot": "HUM", "uniprot_status": "reviewed", "note": "",
         "orthology_type": "ortholog_one2one", "is_target": False},
        {"species": "t", "common": "t", "taxid": 10181, "ensembl_gene": "G",
         "uniprot": "TGT", "uniprot_status": "reviewed", "note": "",
         "orthology_type": "ortholog_one2one", "is_target": True}]}
    ab.fetch_uniprot_seq = lambda acc: full
    ab.uniprot_topology = lambda acc, n: ("o" * 20 + "H" * 20 + "i" * (n - 40), "tm")
    ro.upgrade_partial = lambda *a, **k: (_ for _ in ()).throw(
        AssertionError("upgrade_partial must not be called for a full-length model"))
    try:
        ab.analyze("GENEX", 10181, "/tmp")
    finally:
        (ro.resolve, ab.fetch_uniprot_seq, ab.uniprot_topology,
         ro.upgrade_partial) = saved
    print("no-extra-lookup OK: full-length model triggers no upgrade call")




def _stub_run(ro_mod, orths, seqs, curated=None, tmbed=None, upgrade=None, topology=None,
              reviewed=None, tmbed_ok=True):
    """Run analyze() with the network stubbed out. Returns (resolved, target, models, tmbed_calls)."""
    calls = []
    saved = (ro_mod.resolve, ab.fetch_uniprot_seq, ab.uniprot_topology,
             ab.run_tmbed, ro_mod.upgrade_partial, ab.reviewed_alternatives, ab.tmbed_available)
    ab.tmbed_available = lambda: tmbed_ok
    ab.reviewed_alternatives = reviewed or (lambda gene, taxid: [])
    ro_mod.resolve = lambda g, t: {"gene": g, "target_taxid": t, "orthologs": orths}
    ab.fetch_uniprot_seq = lambda acc: seqs[acc]
    _cur = curated or (lambda acc, n: None)
    ab.uniprot_topology = topology or (
        lambda acc, n: ((_cur(acc, n)), "tm" if _cur(acc, n) else ""))

    def tm(seq, work, tag, acc=None):
        calls.append(acc)
        return (tmbed or (lambda s: "i" * len(s)))(seq)
    ab.run_tmbed = tm
    ro_mod.upgrade_partial = upgrade or (lambda *a, **k: (None, None, ""))
    try:
        resolved, target, targets, models, epi = ab.analyze("GENEX", 10090, "/tmp")
    finally:
        (ro_mod.resolve, ab.fetch_uniprot_seq, ab.uniprot_topology,
         ab.run_tmbed, ro_mod.upgrade_partial, ab.reviewed_alternatives,
         ab.tmbed_available) = saved
    return resolved, target, models, calls


def _orth(common, acc, status, tax, is_t=False, otype="ortholog_one2one"):
    return {"species": common, "common": common, "taxid": tax, "ensembl_gene": "G",
            "uniprot": acc, "uniprot_status": status, "note": "",
            "orthology_type": otype, "is_target": is_t}


def _seq(n, seed=11):
    import random
    r = random.Random(seed)
    return "".join(r.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(n))


def test_target_prefers_complete_over_reviewed_fragment():
    """Mouse CDH1: Q9R1K5 is reviewed at 493 aa, A0A0R4IZW5 unreviewed at 884. Every identity
    is measured against the target, so completeness must outrank curation."""
    import resolve_orthologs as ro
    full = _seq(884)
    seqs = {"FULL_TREMBL": full, "SHORT_REV": full[:493], "HUM": full}
    orths = [_orth("human", "HUM", "reviewed", 9606),
             _orth("mouse", "SHORT_REV", "reviewed", 10090, True),
             _orth("mouse", "FULL_TREMBL", "unreviewed", 10090, True)]
    _, target, _, _ = _stub_run(ro, orths, seqs)
    assert target["uniprot"] == "FULL_TREMBL", target
    # order of the candidate list must not matter
    _, target2, _, _ = _stub_run(ro, list(reversed(orths)), seqs)
    assert target2["uniprot"] == "FULL_TREMBL", target2
    print("target selection OK: 884 aa unreviewed chosen over 493 aa reviewed, order-independent")


def test_fragment_cannot_declare_the_protein_soluble():
    """TMbed's 'no TM' must come from the longest model: a fragment that excludes the helix
    previously marked a single-pass membrane protein as soluble."""
    import resolve_orthologs as ro
    full = _seq(880, seed=3)
    seqs = {"TGT": full, "FRAG": full[:300], "LONG": full}
    orths = [_orth("human", "FRAG", "reviewed", 9606),          # listed first, has no TM
             _orth("rat", "LONG", "unreviewed", 10116),
             _orth("mouse", "TGT", "reviewed", 10090, True)]
    # TMbed finds a helix only in the full-length sequence
    tmbed = lambda s: ("o" * 600 + "H" * 21 + "i" * (len(s) - 621)) if len(s) > 700 else "i" * len(s)
    _, _, models, calls = _stub_run(ro, orths, seqs, tmbed=tmbed)
    assert calls and calls[0] == "LONG", f"TMbed must run on the longest model first, got {calls}"
    assert not any(m.get("soluble") for m in models), "a fragment must not decide solubility"
    assert any(m["regions"]["TM"] is not None for m in models), models
    print("soluble short-circuit OK: longest model decides, membrane protein kept as membrane")


def test_curated_annotation_outranks_prediction_for_solubility():
    """If any model carries curated TM annotation, TMbed is not consulted on the question."""
    import resolve_orthologs as ro
    full = _seq(880, seed=5)
    seqs = {"TGT": full, "HUM": full, "RAB": full}
    orths = [_orth("human", "HUM", "reviewed", 9606),
             _orth("rabbit", "RAB", "unreviewed", 9986),
             _orth("mouse", "TGT", "reviewed", 10090, True)]
    curated = lambda acc, n: ("o" * 600 + "H" * 21 + "i" * (n - 621)) if acc == "HUM" else None
    _, _, models, calls = _stub_run(ro, orths, seqs, curated=curated)
    assert not any(m.get("soluble") for m in models)
    assert "RAB" in calls, "the unannotated model still needs its own topology"
    hum = [m for m in models if m["common"] == "human"][0]
    assert hum["topology_source"] == "UniProt" and hum["regions"]["TM"] is not None
    print("curated-outranks-prediction OK: annotated TM settles solubility")


def test_partial_target_is_rescued_against_ortholog_lengths():
    """A short target is detected using the orthologs already fetched, with no extra lookup
    in the normal case."""
    import resolve_orthologs as ro
    full = _seq(880, seed=9)
    seqs = {"SHORT": full[:400], "FULL": full, "HUM": full, "RAT": full}
    orths = [_orth("human", "HUM", "reviewed", 9606),
             _orth("rat", "RAT", "reviewed", 10116),
             _orth("mouse", "SHORT", "reviewed", 10090, True)]
    seen = []

    def upgrade(acc, status, gene, taxid, ratio=1.3):
        seen.append(acc)
        return ("FULL", "unreviewed", "partial") if acc == "SHORT" else (None, None, "")
    resolved, target, models, _ = _stub_run(ro, orths, seqs, upgrade=upgrade)
    assert seen == ["SHORT"], seen
    assert target["uniprot"] == "FULL", target
    assert any("was partial" in n for n in resolved["notes"]), resolved["notes"]

    # and a normal-length target triggers no lookup at all
    seen.clear()
    orths_ok = [_orth("human", "HUM", "reviewed", 9606),
                _orth("rat", "RAT", "reviewed", 10116),
                _orth("mouse", "FULL", "reviewed", 10090, True)]
    _stub_run(ro, orths_ok, seqs, upgrade=upgrade)
    assert seen == [], f"no upgrade lookup expected for a full-length target, got {seen}"
    print("target rescue OK: 400 aa target upgraded to 880 aa, no lookup when length is normal")




def test_gpi_anchored_protein_is_extracellular_not_soluble():
    """THY1/CD90 is GPI-anchored: a real surface protein with no transmembrane helix.
    Reporting it as 'soluble' is wrong for exactly the antibodies people check before
    buying. Feature table is UniProt P04216 verbatim."""
    feats = [_feat("Signal", 1, 19),
             _feat("Chain", 20, 130, "Thy-1 membrane glycoprotein"),
             _feat("Propeptide", 131, 161, "Removed in mature form"),
             _feat("Lipidation", 130, 130, "GPI-anchor amidated cysteine; alternate")]
    topo = ab._gpi_topology_from_features(feats, 161)
    assert topo is not None and len(topo) == 161
    assert topo[:19] == "S" * 19, topo[:25]                 # signal peptide
    assert set(topo[19:130]) == {"o"}, "the mature chain must be entirely extracellular"
    assert set(topo[130:]) == {"."}, "the cleaved GPI signal is not part of the protein"
    assert "H" not in topo and "i" not in topo

    # a protein with no GPI annotation, and one with a helix as well, both decline this path
    assert ab._gpi_topology_from_features([_feat("Chain", 1, 100)], 100) is None
    both = feats + [_feat("Transmembrane", 60, 80, "Helical")]
    assert ab._gpi_topology_from_features(both, 161) is None, "a helix is described by the TM path"
    print("GPI topology OK: THY1 mature chain 20-130 all extracellular, signal and "
          "propeptide excluded")


def test_gpi_run_reports_extracellular_and_not_soluble():
    """End to end: a GPI-anchored protein must come out with an extracellular identity."""
    import resolve_orthologs as ro
    full = _seq(161, seed=21)
    seqs = {"HUM": full, "TGT": full[:120] + _seq(41, seed=22)}
    orths = [_orth("human", "HUM", "reviewed", 9606),
             _orth("mouse", "TGT", "reviewed", 10090, True)]
    gpi = "S" * 19 + "o" * 111 + "." * 31
    _, _, models, calls = _stub_run(
        ro, orths, seqs,
        topology=lambda acc, n: ((gpi, "gpi") if acc == "HUM" else (None, "")))
    hum = [m for m in models if m["common"] == "human"][0]
    assert not hum.get("soluble"), "a GPI-anchored protein is not soluble"
    assert hum["regions"]["extracellular"] is not None, hum["regions"]
    assert hum["regions"]["TM"] is None, "there is no helix to report"
    assert hum["topology_source"] == "UniProt (GPI-anchored)", hum["topology_source"]
    assert calls == [], "a curated GPI anchor settles topology; TMbed must not run"
    print("GPI end-to-end OK: extracellular %.1f%%, no TM row, TMbed skipped"
          % hum["regions"]["extracellular"])


def test_no_compara_ortholog_message():
    """The message must say what actually happened rather than imply a retry will fix it."""
    import resolve_orthologs as ro
    saved = ro.resolve
    ro.resolve = lambda g, t: {"gene": g, "target_taxid": t, "orthologs": [
        _orth("human", "HUM", "reviewed", 9606)]}          # no target ortholog at all
    try:
        ab.analyze("ACTB", 9031, "/tmp")
    except ab.AbScoutError as e:
        msg = str(e)
        assert "No Compara ortholog" in msg and "check it manually" in msg, msg
        assert "does not fall back to matching on gene name" in msg, msg
    else:
        raise AssertionError("a missing target ortholog must raise")
    finally:
        ro.resolve = saved
    print("no-ortholog message OK: states the cause, declines name-based guessing")




# Real UniProt sequences (fetched this session) for the cases found in the panel run.
REAL = {'P07766': 'MQSGTHWRVLGLCLLSVGVWGQDGNEEMGGITQTPYKVSISGTTVILTCPQYPGSEILWQHNDKNIGGDEDDKNIGSDEDHLSLKEFSELEQSGYYVCYPRGSKPEDANFYLYLRARVCENCMEMDVMSVATIVIVDICITGGLLLLVYYWSKNRKAKAKPVTRGAGAGGRQRGQNKERPPPVPNPDYEPIRKGQRDLYSGLNQRRI', 'P22646': 'MRWNTFWGILCLSLLAVGTCQDDAENIEYKVSISGTSVELTCPLDSDENLKWEKNGQELPQKHDKHLVLQDFSEVEDSGYYVCYTPASNKNTYLYLKARVCEYCVEVDLTAVAIIIIVDICITLGLLMVIYYWSKNRKAKAKPVTRGTGAGSRPRGQNKERPPPVPNPDYEPIRKGQRDLYSGLNQRAV', 'P27597': 'MQSRNLWRILGLCLLSVGAWGQDEDFKASDDLTSISPEKRFKVSISGTEVVVTCPDVFGYDNIKWEKNDNLVEGASNRELSQKEFSEVDDSGYYACYADSIKEKSYLYLRARVCANCIEVNLMAVVTIIVADICLTLGLLLMVYYWSKTRKANAKPVMRGTGAGSRPRGQNKEKPPPVPNPDYEPIRKGQQDLYSGLNQRGI', 'A0A8I3MGF0': 'CLVPLENVYLQFKLKPRLLSYGVFLFSCSDDLTSISPEKRFKVSISGTEVVVTCPDVFGYDNIKWEKNDNLVEGASNRELSQKEFSEVDDSGYYACYADSVKEKSYLYLRARVCANCIEVNLMAVVTIIVADICLTLGLLLMVYYWSKTRKANAKPVMRGTGAGSRPRGKTMELSWKSLTNGDVQPGPGWVASPQLGQDNFYRTFTVWPNSGQTVSHLRLLLSFLP', 'Q9TUF9': 'MRAGTLWRVLALWLLSVAAWGQEDDDHADDYTQKLFTVSISGTRVVLTCPVEAEGGDIHWERDEKSLPNTKKELDLTDFSEMEHSGYYSCYVGTKNKENEHILYLKARVCEACMEVDLTTVASIVVADVCVTLGLLLLVYYWSKNRKAKCKPVTRGAGAGGRPRGQNKERPPPVPNPDYEPIRKGQRDLYSGLNQRGI', 'A0A5F9CJ85': 'APPGSLVTFITGLVLLSLACLFWNVPGAPHHHHHQFPCPLSHCPVFAVSISGNHVVLTCPVEADAGDIHWGKDEKSLPDTDKELDLENFSEMEHSGYYSCYVGTKNKENQHILYLKARVCEACMEVDLTTVASIVVADVCVTLGLLLLVYYWSKNRKAKCKPVMRGAGAGGRPRGEPQPLSIPPQPIRKGQRDLYSGLNQRGI', 'P12830': 'MGPWSRSLSALLLLLQVSSWLCQEPEPCHPGFDAESYTFTVPRRHLERGRVLGRVNFEDCTGRQRTAYFSLDTRFKVGTDGVITVKRPLRFHNPQIHFLVYAWDSTYRKFSTKVTLNTVGHHHRPPPHQASVSGIQAELLTFPNSSPGLRRQKRDWVIPPISCPENEKGPFPKNLVQIKSNKDKEGKVFYSITGQGADTPPVGVFIIERETGWLKVTEPLDRERIATYTLFSHAVSSNGNAVEDPMEILITVTDQNDNKPEFTQEVFKGSVMEGALPGTSVMEVTATDADDDVNTYNAAIAYTILSQDPELPDKNMFTINRNTGVISVVTTGLDRESFPTYTLVVQAADLQGEGLSTTATAVITVTDTNDNPPIFNPTTYKGQVPENEANVVITTLKVTDADAPNTPAWEAVYTILNDDGGQFVVTTNPVNNDGILKTAKGLDFEAKQQYILHVAVTNVVPFEVSLTTSTATVTVDVLDVNEAPIFVPPEKRVEVSEDFGVGQEITSYTAQEPDTFMEQKITYRIWRDTANWLEINPDTGAISTRAELDREDFEHVKNSTYTALIIATDNGSPVATGTGTLLLILSDVNDNAPIPEPRTIFFCERNPKPQVINIIDADLPPNTSPFTAELTHGASANWTIQYNDPTQESIILKPKMALEVGDYKINLKLMDNQNKDQVTTLEVSVCDCEGAAGVCRKAQPVEAGLQIPAILGILGGILALLILILLLLLFLRRRAVVKEPLLPPEDDTRDNVYYYDEEGGGEEDQDFDLSQLHRGLDARPEVTRNDVAPTLMSVPRYLPRPANPDEIGNFIDENLKAADTDPTAPPYDSLLVFDYEGSGSEAASLSSLNSSESDKDQDYDYLNEWGNRFKKLADMYGGGEDD', 'Q9UM11': 'MDQDYERRLLRQIVIQNENTMPRVTEMRRTLTPASSPVSSPSKHGDRFIPSRAGANWSVNFHRINENEKSPSQNRKAKDATSDNGKDGLAYSALLKNELLGAGIEKVQDPQTEDRRLQPSTPEKKGLFTYSLSTKRSSPDDGNDVSPYSLSPVSNKSQKLLRSPRKPTRKISKIPFKVLDAPELQDDFYLNLVDWSSLNVLSVGLGTCVYLWSACTSQVTRLCDLSVEGDSVTSVGWSERGNLVAVGTHKGFVQIWDAAAGKKLSMLEGHTARVGALAWNAEQLSSGSRDRMILQRDIRTPPLQSERRLQGHRQEVCGLKWSTDHQLLASGGNDNKLLVWNHSSLSPVQQYTEHLAAVKAIAWSPHQHGLLASGGGTADRCIRFWNTLTGQPLQCIDTGSQVCNLAWSKHANELVSTHGYSQNQILVWKYPSLTQVAKLTGHSYRVLYLAMSPDGEAIVTGAGDETLRFWNVFSKTRSTKVKWESVSVLNLFTRIR'}


def test_reviewed_preferred_for_dog_cd3e_target():
    """Real case: Compara's xref gave dog CD3E as TrEMBL A0A8I3MGF0, whose mispredicted
    C-terminus drops the conserved epitope; the reviewed P27597 has it."""
    import resolve_orthologs as ro
    orths = [_orth("human", "P07766", "reviewed", 9606),
             _orth("mouse", "P22646", "reviewed", 10090),
             _orth("dog", "A0A8I3MGF0", "unreviewed", 9615, True)]
    resolved, target, _, _ = _stub_run(
        ro, orths, REAL, reviewed=lambda g, t: ["P27597"] if t == 9615 else [])
    assert target["uniprot"] == "P27597", target
    assert any("used reviewed P27597" in n for n in resolved["notes"]), resolved["notes"]
    print("reviewed preference OK: dog target A0A8I3MGF0 -> P27597")


def test_reviewed_preferred_for_rabbit_cd3e_model():
    """Real case: rabbit CD3E A0A5F9CJ85 replaces the epitope with EPQPLSIPPQP."""
    import resolve_orthologs as ro
    orths = [_orth("human", "P07766", "reviewed", 9606),
             _orth("rabbit", "A0A5F9CJ85", "unreviewed", 9986),
             _orth("mouse", "P22646", "reviewed", 10090, True)]
    _, _, models, _ = _stub_run(
        ro, orths, REAL, reviewed=lambda g, t: ["Q9TUF9"] if t == 9986 else [])
    rab = [m for m in models if m["common"] == "rabbit"][0]
    assert rab["uniprot"] == "Q9TUF9", rab
    print("reviewed preference OK: rabbit model A0A5F9CJ85 -> Q9TUF9")


def test_synonym_collision_is_refused():
    """Real case: the CDH1 symbol also matches FZR1 (Q9UM11). Offered as the only reviewed
    alternative, it must be refused; offered alongside E-cadherin, E-cadherin must win."""
    import resolve_orthologs as ro
    ecad = REAL["P12830"]
    seqs = dict(REAL, RAB_ECAD=ecad[:400] + ("A" if ecad[400] != "A" else "G") + ecad[401:])
    orths = [_orth("human", "P12830", "reviewed", 9606),
             _orth("rabbit", "RAB_ECAD", "unreviewed", 9986),
             _orth("mouse", "P12830", "reviewed", 10090, True)]
    resolved, _, models, _ = _stub_run(ro, orths, seqs, reviewed=lambda g, t: ["Q9UM11"])
    rab = [m for m in models if m["common"] == "rabbit"][0]
    assert rab["uniprot"] == "RAB_ECAD", "FZR1 must never replace an E-cadherin ortholog"
    assert any("Q9UM11 exists but resembles the other species less" in n
               for n in resolved["notes"]), resolved["notes"]

    _, _, models2, _ = _stub_run(ro, orths, seqs, reviewed=lambda g, t: ["Q9UM11", "P12830"])
    rab2 = [m for m in models2 if m["common"] == "rabbit"][0]
    assert rab2["uniprot"] == "P12830", rab2
    print("synonym guard OK: FZR1 refused alone, E-cadherin chosen when both are offered")


def test_partial_rescue_cannot_swap_in_a_different_protein():
    """The length rescue searches by gene symbol too, so it carries the same guard."""
    import resolve_orthologs as ro
    ecad = REAL["P12830"]
    seqs = dict(REAL, FRAG=ecad[:300], TGT=ecad)
    orths = [_orth("human", "FRAG", "reviewed", 9606),
             _orth("rat", "P12830", "reviewed", 10116),
             _orth("mouse", "TGT", "reviewed", 10090, True)]
    # the symbol search hands back FZR1 as "the longer entry"
    resolved, _, models, _ = _stub_run(
        ro, orths, seqs, upgrade=lambda acc, st, g, t, ratio=1.3: ("Q9UM11", "reviewed", "x"))
    assert not any(m["uniprot"] == "Q9UM11" for m in models), models
    assert any("does not resemble the other species" in n for n in resolved["notes"])
    # and the genuine full-length entry is accepted
    _, _, models2, _ = _stub_run(
        ro, orths, seqs, upgrade=lambda acc, st, g, t, ratio=1.3: ("P12830", "reviewed", "x"))
    hum = [m for m in models2 if m["common"] == "human"][0]
    assert hum["uniprot"] == "P12830", hum
    print("rescue guard OK: FZR1 refused as a 'longer entry', real E-cadherin accepted")


def test_non_ortholog_is_flagged_for_manual_check():
    """Whatever path brought it in, a sequence resembling no other species is flagged."""
    import resolve_orthologs as ro
    orths = [_orth("human", "P12830", "reviewed", 9606),
             _orth("rabbit", "Q9UM11", "reviewed", 9986),      # FZR1 masquerading as CDH1
             _orth("mouse", "P12830", "reviewed", 10090, True)]
    resolved, _, _, _ = _stub_run(ro, orths, REAL)
    flagged = [n for n in resolved["notes"] if "may not be the ortholog" in n]
    assert flagged and "rabbit" in flagged[0] and "check it manually" in flagged[0], flagged
    print("non-ortholog flag OK:", flagged[0][:70] + "...")


def test_absent_epitope_and_notes_reach_both_reports():
    """A truncated ortholog shows 'absent' (never a blank that reads as 0%), and the run's
    notes appear in the terminal report and the HTML page."""
    full = _seq(120, seed=31)
    ep = ab.epitope_across_species([("human", full), ("rabbit", full[30:])], full,
                                   "human", (2, 16))
    rab = [r for r in ep["rows"] if r["common"] == "rabbit"][0]
    assert rab["absent"] and rab["length"] == 0, rab
    epi = {"ref": "human", "range": [2, 16], "rows": ep["rows"], "alignment": ep["alignment"],
           "substitutions": ep["substitutions"], "glyco": ep["glyco"]}
    notes = ["rabbit: used reviewed X instead of unreviewed Y"]
    report, _ = ab.rank_and_report("ACTB", TARGET, [TARGET], MODELS, epi, notes=notes)
    assert "absent" in report and "this is not 0%" in report, report
    assert "Notes on this run" in report and notes[0] in report
    page = ab.render_html("ACTB", TARGET, MODELS, epi, notes=notes)
    assert "Notes on this run" in page and "used reviewed X" in page
    assert ">absent</td>" in page
    print("reports OK: 'absent' shown instead of a blank; notes in terminal and HTML")


def test_topology_segments():
    """Contiguous runs per region in 1-based coordinates; signal and '.' are not regions."""
    t = "SSSooooHHHiiiHHHoo..."
    seg = ab.topology_segments(t)
    assert seg["extracellular"] == [(4, 7), (17, 18)], seg
    assert seg["TM"] == [(8, 10), (14, 16)], seg
    assert seg["intracellular"] == [(11, 13)], seg
    assert ab.topology_segments("") == {"extracellular": [], "TM": [], "intracellular": []}
    print("topology segments OK: multi-pass runs split, signal/unassigned left out")


def test_region_alignments_real_cd3e():
    """Real human/mouse/dog CD3E on the curated human topology: three regions, target
    coordinates reported, and the conserved ERPPPVPNPDYEP motif sits in the intracellular
    alignment of every species."""
    hum = REAL["P07766"]
    topo = "S" * 22 + "o" * 104 + "H" * 26 + "i" * 55
    models = [("human", hum), ("mouse", REAL["P22646"])]
    r = ab.region_alignments(models, REAL["P27597"], "human", topo)
    assert r["ref"] == "human"
    labels = [g["label"] for g in r["regions"]]
    assert labels == ["Extracellular", "Transmembrane", "Intracellular"], labels
    tm = r["regions"][1]["segments"][0]
    assert tm["range"] == [127, 152], tm["range"]
    assert all(len(v) == 26 for v in tm["alignment"].values()), tm["alignment"]
    intra = r["regions"][2]["segments"][0]
    tgt_row = [x for x in intra["rows"] if x["common"] == "target"][0]
    assert tgt_row["coordinates"] and tgt_row["length"] > 40, tgt_row
    for sp in ("human", "mouse"):
        assert "PPPVPNPDYEP" in intra["alignment"][sp].replace("-", ""), sp
    html = ab.render_html("CD3E", {"common": "dog", "species": "canis_lupus_familiaris", "uniprot": "P27597", "uniprot_status": "reviewed", "taxid": 9615, "length": len(REAL["P27597"])},
                          [], None, regions=r)
    for h in ("<h3>Extracellular</h3>", "<h3>Transmembrane</h3>", "<h3>Intracellular</h3>",
              "human 127&ndash;152"):
        assert h in html, h
    print("region alignments OK: Extracellular / Transmembrane / Intracellular sections "
          "rendered, target coordinates %s" % tgt_row["coordinates"])


def test_derived_cache_is_per_version():
    """Real case: human THY1 kept a "no topology" answer cached before GPI support existed,
    so every later version sent it to TMbed. Parsed answers must not cross versions; raw
    sequences may."""
    import pathlib, tempfile
    old_dir, old_en = ab.CACHE_DIR, ab.CACHE_ENABLED
    with tempfile.TemporaryDirectory() as tmpdir:
        ab.CACHE_DIR, ab.CACHE_ENABLED = pathlib.Path(tmpdir), True
        try:
            for stale in ("topo_uniprot", "topo_uniprot@0.4.3"):      # what older versions wrote
                d = pathlib.Path(tmpdir) / stale
                d.mkdir(parents=True, exist_ok=True)
                (d / "P04216.txt").write_text("")
            assert ab.cache_get("topo_uniprot", "P04216") is None, "stale parse was served"
            ab.cache_put("topo_uniprot", "P04216", "SSo")
            assert ab.cache_get("topo_uniprot", "P04216") == "SSo"
            ab.cache_put("seq", "P04216", "MNLAIS")
            assert (pathlib.Path(tmpdir) / "seq" / "P04216.txt").exists(), "raw data is unversioned"
        finally:
            ab.CACHE_DIR, ab.CACHE_ENABLED = old_dir, old_en
    print("cache OK: parsed UniProt answers are per version; sequences are shared")


def test_gpi_alternate_description_real_thy1():
    """The real human THY1 entry says "GPI-anchor amidated cysteine; alternate"."""
    f = [_feat("Signal", 1, 19), _feat("Chain", 20, 130, "Thy-1 membrane glycoprotein"),
         _feat("Propeptide", 131, 161, "Removed in mature form"),
         _feat("Lipidation", 130, 130, "GPI-anchor amidated cysteine; alternate")]
    topo = ab._gpi_topology_from_features(f, 161)
    assert topo and topo.count("o") == 111 and topo[:19] == "S" * 19 and topo.endswith("." * 31), topo
    print("GPI OK: P04216 as UniProt annotates it -> 111 extracellular residues")


def test_deleted_uniprot_entry_is_reported_not_used():
    """Real case: naked mole-rat THY1 resolves to A0A0P6J0R8, deleted from UniProtKB. The empty
    FASTA used to become a 0 aa target and every species failed with "zero length"."""
    import resolve_orthologs as ro

    class R:
        def __init__(self, text="", js=None):
            self.text, self._js, self.status_code = text, js, 200
        def json(self):
            return self._js
    inactive = {"entryType": "Inactive", "inactiveReason": {
        "inactiveReasonType": "DELETED", "deletedReason": "Not part of a reference proteome"}}
    saved, saved_en = ab._uniprot_get, ab.CACHE_ENABLED
    ab._uniprot_get = lambda url, tries=4: R(js=inactive) if url.endswith(".json") else R("")
    ab.CACHE_ENABLED = False
    try:
        ab.fetch_uniprot_seq("A0A0P6J0R8")
    except ab.UniProtUnavailable as e:
        assert str(e) == "deleted from UniProtKB (Not part of a reference proteome)", str(e)
    else:
        raise AssertionError("an empty FASTA must not become a 0 aa sequence")
    finally:
        ab._uniprot_get, ab.CACHE_ENABLED = saved, saved_en

    class Seqs(dict):
        def __getitem__(self, acc):
            if acc == "A0A0P6J0R8":
                raise ab.UniProtUnavailable("deleted from UniProtKB (Not part of a reference "
                                            "proteome)")
            return dict.__getitem__(self, acc)
    orths = [_orth("human", "P04216", "reviewed", 9606),
             _orth("heterocephalus_glaber_female", "A0A0P6J0R8", "unreviewed", 10181, True)]
    try:
        _stub_run(ro, orths, Seqs({"P04216": _seq(161)}))
    except ab.AbScoutError as e:
        msg = str(e)
        assert "A0A0P6J0R8" in msg and "deleted from UniProtKB" in msg, msg
        assert "check it manually" in msg, msg
    else:
        raise AssertionError("a run with no usable target sequence must stop")
    print("deleted entry OK: reported as deleted from UniProtKB, not run as a 0 aa protein")


def test_without_tmbed_uncurated_species_keeps_whole_protein_numbers():
    """The web version has no TMbed. A species with no curated topology must still be
    compared (whole protein + epitope), not skipped, and TMbed must never be called."""
    import resolve_orthologs as ro
    hum, mou, tgt = REAL["P07766"], REAL["P22646"], REAL["P27597"]
    topo_h = "S" * 22 + "o" * 104 + "H" * 26 + "i" * 55
    orths = [_orth("human", "P07766", "reviewed", 9606),
             _orth("mouse", "P22646", "reviewed", 10090),
             _orth("dog", "P27597", "reviewed", 9615, True)]
    resolved, _, models, calls = _stub_run(
        ro, orths, REAL, tmbed_ok=False,
        topology=lambda acc, n: ((topo_h, "tm") if acc == "P07766" else (None, "")))
    assert calls == [], "TMbed must not be called when it is not available"
    mou_r = [m for m in models if m["common"] == "mouse"][0]
    assert mou_r["regions"] == {"extracellular": None, "TM": None, "intracellular": None}
    assert mou_r["overall"] and mou_r["topology_source"] == "not annotated in UniProt", mou_r
    assert not mou_r.get("soluble"), "unknown topology is not soluble"
    hum_r = [m for m in models if m["common"] == "human"][0]
    assert hum_r["regions"]["TM"] is not None, hum_r["regions"]
    assert resolved["region_alignments"], "human's curated topology still anchors the regions"
    note = " ".join(resolved["notes"])
    assert "no curated membrane topology in UniProt for mouse (P22646)" in note, note
    assert "region-by-region" in note, note
    print("no-TMbed OK: mouse kept with whole-protein %.1f%%, human regions intact, notes say why"
          % mou_r["overall"])


def test_without_tmbed_nothing_curated_is_unknown_not_soluble():
    import resolve_orthologs as ro
    orths = [_orth("human", "P07766", "reviewed", 9606),
             _orth("dog", "P27597", "reviewed", 9615, True)]
    resolved, _, models, calls = _stub_run(ro, orths, REAL, tmbed_ok=False,
                                           topology=lambda acc, n: (None, ""))
    assert calls == [] and models and not any(m.get("soluble") for m in models), models
    note = " ".join(resolved["notes"])
    assert "membrane topology unknown" in note and "not assumed to be soluble" in note, note
    assert not resolved["region_alignments"], "no topology, no region sections"
    print("no-TMbed OK: nothing curated -> 'topology unknown', never 'soluble'")


def test_no_target_message_says_why():
    """Real case: TP53 in the Damara mole-rat (885580). Ensembl has that genome but not in
    Compara, so no gene can resolve there; the message should say so rather than suggest a
    per-gene manual check, and the web page must not suggest a command-line option."""
    import resolve_orthologs as ro
    human_only = {"orthologs": [_orth("human", "P04637", "reviewed", 9606)],
                  "homology_query_ok": True}
    saved = (ro.compara_genomes, ab.IN_BROWSER)
    try:
        ro.compara_genomes = lambda t: [("Damara mole rat", "fukomys_damarensis", False)]
        m = ab.no_target_message("TP53", 885580, human_only)
        assert "Damara mole rat (fukomys_damarensis)" in m and "Compara's gene trees" in m, m
        assert "--timeout" not in m, "a timeout cannot explain a genome outside Compara"

        ro.compara_genomes = lambda t: [("Naked mole-rat female", "heterocephalus_glaber_female",
                                         True)]
        m = ab.no_target_message("XYZ1", 10181, human_only)
        assert m.startswith("No Compara ortholog for XYZ1 in taxid 10181") and "--timeout" in m, m
        ab.IN_BROWSER = True
        m = ab.no_target_message("XYZ1", 10181, human_only)
        assert "--timeout" not in m and "run it again" in m, m
        ab.IN_BROWSER = False

        ro.compara_genomes = lambda t: []
        m = ab.no_target_message("TP53", 123456789, human_only)
        assert "no genome for taxid 123456789" in m, m

        m = ab.no_target_message("TP53", 9615, {"orthologs": [], "homology_query_ok": False})
        assert "Ensembl returned no orthologs for TP53" in m, m

        no_xref = {"orthologs": [dict(_orth("dog", None, None, 9615, True),
                                      ensembl_gene="ENSCAFG00845012345")],
                   "homology_query_ok": True}
        m = ab.no_target_message("TP53", 9615, no_xref)
        assert "Compara has a TP53 ortholog" in m and "no UniProt entry" in m, m

        # and analyze() really stops with it
        ro.compara_genomes = lambda t: [("Damara mole rat", "fukomys_damarensis", False)]
        try:
            _stub_run(ro, [_orth("human", "P07766", "reviewed", 9606)], REAL)
        except ab.AbScoutError as e:
            assert "Compara's gene trees" in str(e), str(e)
        else:
            raise AssertionError("no target must stop the run")
    finally:
        ro.compara_genomes, ab.IN_BROWSER = saved
    print("no-target OK: species outside Compara, unknown taxon, failed query and missing "
          "UniProt entry each get their own message")


class _DeletedSeqs(dict):
    """UniProt sequences where some accessions behave like deleted entries (empty FASTA)."""
    def __init__(self, base, deleted):
        super().__init__(base)
        self.deleted = set(deleted)

    def __getitem__(self, acc):
        if acc in self.deleted:
            raise ab.UniProtUnavailable("deleted from UniProtKB (Not part of a reference "
                                        "proteome)")
        return dict.__getitem__(self, acc)


def _cd3e_topology(asked):
    hum = "S" * 22 + "o" * 104 + "H" * 26 + "i" * 55
    def topo(acc, n):
        asked.append(acc)
        return (hum, "tm") if acc == "P07766" else (None, "")
    return topo


def test_deleted_target_uses_ensembl_then_swissprot():
    """Real case: dog TP53's Compara link (A0A8C0RLV2) is deleted from UniProt. Shown with real
    CD3E sequences: Ensembl's sequence of the same Compara ortholog takes the deleted entry's
    place, and the reviewed Swiss-Prot entry (P27597) then replaces it as the curated one.
    Nothing is matched on name alone: the Swiss-Prot entry must resemble the other species
    at least as well as Compara's own sequence does."""
    import resolve_orthologs as ro
    dog = dict(_orth("dog", "A0A_DELETED", "unreviewed", 9615, True),
               ensembl_protein="ENSCAFP_T")
    orths = [_orth("human", "P07766", "reviewed", 9606),
             _orth("mouse", "P22646", "reviewed", 10090), dog]
    asked, saved = [], ab.fetch_ensembl_seq
    ab.fetch_ensembl_seq = lambda pid: REAL["A0A8I3MGF0"]     # the gene model Ensembl holds
    try:
        resolved, target, models, _ = _stub_run(
            ro, orths, _DeletedSeqs(REAL, ["A0A_DELETED"]), tmbed_ok=False,
            topology=_cd3e_topology(asked),
            reviewed=lambda g, t: ["P27597"] if t == 9615 else [])
    finally:
        ab.fetch_ensembl_seq = saved
    assert target["uniprot"] == "P27597" and target["uniprot_status"] == "reviewed", target
    notes = " ".join(resolved["notes"])
    assert "UniProt A0A_DELETED is deleted from UniProtKB" in notes, notes
    assert "Ensembl's sequence of the Compara ortholog, ENSCAFP_T" in notes, notes
    assert "used reviewed P27597 instead of Ensembl's sequence Ensembl:ENSCAFP_T" in notes, notes
    assert not any(str(a).startswith("Ensembl:") for a in asked), "UniProt asked about an Ensembl ID"
    print("deleted target OK: A0A_DELETED -> Ensembl ENSCAFP_T -> reviewed P27597")


def test_deleted_target_without_swissprot_keeps_ensembl_sequence():
    """Naked mole-rat THY1 case: no reviewed entry exists, so the run goes ahead on Ensembl's
    sequence and every report labels it as such."""
    import resolve_orthologs as ro
    dog = dict(_orth("dog", "A0A_DELETED", "unreviewed", 9615, True),
               ensembl_protein="ENSCAFP_T")
    orths = [_orth("human", "P07766", "reviewed", 9606),
             _orth("mouse", "P22646", "reviewed", 10090), dog]
    asked, saved = [], ab.fetch_ensembl_seq
    ab.fetch_ensembl_seq = lambda pid: REAL["P27597"]
    try:
        resolved, target, models, _ = _stub_run(
            ro, orths, _DeletedSeqs(REAL, ["A0A_DELETED"]), tmbed_ok=False,
            topology=_cd3e_topology(asked))
    finally:
        ab.fetch_ensembl_seq = saved
    assert target["uniprot"] == "Ensembl:ENSCAFP_T" and target["uniprot_status"] == "ensembl"
    assert target["replaced_uniprot"] == "A0A_DELETED"
    assert {m["common"] for m in models} == {"human", "mouse"}, models
    html = ab.render_html("CD3E", target, models, None, notes=resolved["notes"])
    assert "Ensembl:ENSCAFP_T" in html, "the report must say where the target came from"
    print("deleted target OK: no Swiss-Prot entry -> run on Ensembl:ENSCAFP_T, labelled")


def test_deleted_model_uses_ensembl_sequence_without_uniprot_lookups():
    import resolve_orthologs as ro
    mouse = dict(_orth("mouse", "A0A_GONE", "unreviewed", 10090), ensembl_protein="ENSMUSP_T")
    orths = [_orth("human", "P07766", "reviewed", 9606), mouse,
             _orth("dog", "P27597", "reviewed", 9615, True)]
    asked, saved = [], ab.fetch_ensembl_seq
    ab.fetch_ensembl_seq = lambda pid: REAL["P22646"]
    try:
        resolved, target, models, _ = _stub_run(
            ro, orths, _DeletedSeqs(REAL, ["A0A_GONE"]), tmbed_ok=False,
            topology=_cd3e_topology(asked))
    finally:
        ab.fetch_ensembl_seq = saved
    mou = [m for m in models if m["common"] == "mouse"][0]
    assert mou["uniprot"] == "Ensembl:ENSMUSP_T", mou
    assert mou["topology_source"] == "no UniProt entry (Ensembl sequence)", mou
    assert mou["overall"] and not any(str(a).startswith("Ensembl:") for a in asked), asked
    print("deleted model OK: mouse compared on Ensembl's sequence, UniProt not asked about it")


def test_no_sequence_anywhere_is_reported():
    import resolve_orthologs as ro
    dog = dict(_orth("dog", "A0A_DELETED", "unreviewed", 9615, True),
               ensembl_protein="ENSCAFP_T")
    saved = ab.fetch_ensembl_seq

    def none(pid):
        raise ab.UniProtUnavailable(f"Ensembl returned no sequence for {pid}")
    ab.fetch_ensembl_seq = none
    try:
        _stub_run(ro, [_orth("human", "P07766", "reviewed", 9606), dog],
                  _DeletedSeqs(REAL, ["A0A_DELETED"]), tmbed_ok=False)
    except ab.AbScoutError as e:
        assert "neither UniProt nor Ensembl gave a sequence" in str(e), str(e)
    else:
        raise AssertionError("must stop when no sequence can be had")
    finally:
        ab.fetch_ensembl_seq = saved
    print("no-sequence OK: stops and says neither source had one")


if __name__ == "__main__":
    test_region_table_no_recommendation()
    test_handles_missing_region_and_paralogs()
    test_uniprot_topology_parsing()
    test_epitope_across_species()
    test_epitope_indel_is_reported()
    test_peptide_match_floor()
    test_epitope_range_validation()
    test_glyco_sequon_detection()
    test_glyco_gain_in_epitope()
    test_epitope_substitution_detail()
    test_signal_peptide_numbering()
    test_rank_report_epitope_table()
    test_locate_peptide()
    test_detect_epitope_species()
    test_run_tmbed_rejects_short_sequence()
    test_soluble_protein_report()
    print("-" * 60)
    test_cache_roundtrip()
    test_coverage_reported()
    test_low_coverage_note()
    test_batch_csv()
    test_nterm_processing_real_feature_tables()
    test_epitope_ref_must_be_a_compared_species()
    test_range_used_is_reported()
    test_batch_survives_failure_and_carries_notes()
    test_duplicate_species_rejected()
    test_one_ortholog_chosen_per_species()
    test_short_reviewed_model_is_rescued_not_skipped()
    test_full_length_model_costs_no_extra_lookup()
    test_target_prefers_complete_over_reviewed_fragment()
    test_fragment_cannot_declare_the_protein_soluble()
    test_curated_annotation_outranks_prediction_for_solubility()
    test_partial_target_is_rescued_against_ortholog_lengths()
    test_gpi_anchored_protein_is_extracellular_not_soluble()
    test_gpi_run_reports_extracellular_and_not_soluble()
    test_no_compara_ortholog_message()
    test_reviewed_preferred_for_dog_cd3e_target()
    test_reviewed_preferred_for_rabbit_cd3e_model()
    test_synonym_collision_is_refused()
    test_partial_rescue_cannot_swap_in_a_different_protein()
    test_non_ortholog_is_flagged_for_manual_check()
    test_absent_epitope_and_notes_reach_both_reports()
    test_topology_segments()
    test_region_alignments_real_cd3e()
    test_derived_cache_is_per_version()
    test_gpi_alternate_description_real_thy1()
    test_deleted_uniprot_entry_is_reported_not_used()
    test_without_tmbed_uncurated_species_keeps_whole_protein_numbers()
    test_without_tmbed_nothing_curated_is_unknown_not_soluble()
    test_no_target_message_says_why()
    test_deleted_target_uses_ensembl_then_swissprot()
    test_deleted_target_without_swissprot_keeps_ensembl_sequence()
    test_deleted_model_uses_ensembl_sequence_without_uniprot_lookups()
    test_no_sequence_anywhere_is_reported()
    print("\nALL ABSCOUT TESTS PASSED")
