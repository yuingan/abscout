#!/usr/bin/env python3
"""Mock tests for resolve_orthologs.py — no network. Run: python test_resolve_orthologs.py"""
import importlib.util

spec = importlib.util.spec_from_file_location("ro", "resolve_orthologs.py")
ro = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ro)

# --- mock Ensembl homology response (condensed), incl. one-to-many fish orthologs ---
CONDENSED = {"data": [{"id": "ENSG00000146648", "homologies": [
    {"species": "mus_musculus", "id": "ENSMUSG00000020122", "protein_id": "ENSMUSP_MOUSE",
     "type": "ortholog_one2one"},
    {"species": "rattus_norvegicus", "id": "ENSRNOG00000004332", "protein_id": "ENSRNOP_RAT",
     "type": "ortholog_one2one"},
    {"species": "danio_rerio", "id": "ENSDARG_A", "protein_id": "ENSDARP_A",
     "type": "ortholog_one2many"},
    {"species": "danio_rerio", "id": "ENSDARG_B", "protein_id": "ENSDARP_B",
     "type": "ortholog_one2many"},
]}]}


def test_parse_homologies_condensed():
    src, homs = ro.parse_homologies(CONDENSED)
    assert src == "ENSG00000146648", src
    assert len(homs) == 4
    assert {h["species"] for h in homs} == {"mus_musculus", "rattus_norvegicus", "danio_rerio"}
    assert sum(1 for h in homs if h["species"] == "danio_rerio") == 2, "one-to-many fish lost"
    assert homs[0]["gene_id"] == "ENSMUSG00000020122"
    assert homs[0]["protein_id"] == "ENSMUSP_MOUSE"
    print("parse_homologies (condensed) OK")


def test_parse_homologies_nested():
    nested = {"data": [{"id": "ENSG1", "homologies": [
        {"species": "mus_musculus", "type": "ortholog_one2one",
         "target": {"id": "ENSMUSG1", "protein_id": "ENSMUSP1", "species": "mus_musculus"}},
    ]}]}
    src, homs = ro.parse_homologies(nested)
    assert homs[0]["gene_id"] == "ENSMUSG1" and homs[0]["protein_id"] == "ENSMUSP1"
    print("parse_homologies (nested/full) OK")


def test_parse_xref():
    xref = [{"primary_id": "P00533", "dbname": "Uniprot/SWISSPROT", "display_id": "EGFR_HUMAN"}]
    assert ro.parse_xref_uniprot(xref) == "P00533"
    assert ro.parse_xref_uniprot([]) is None
    assert ro.parse_xref_uniprot(None) is None
    print("parse_xref_uniprot OK")


def test_resolve_end_to_end():
    XREF = {"ENSG00000146648": "P00533", "ENSMUSP_MOUSE": "Q01279",
            "ENSRNOP_RAT": "Q9QX70", "ENSDARP_A": "DARA1", "ENSDARP_B": "DARB1"}

    def fake_get(path, params=None):
        if path.startswith("/homology/symbol/human/EGFR"):
            return CONDENSED
        if path.startswith("/xrefs/id/"):
            eid = path.rsplit("/", 1)[-1]
            return [{"primary_id": XREF[eid]}] if eid in XREF else []
        return None

    ro.ensembl_get = fake_get
    ro.uniprot_by_gene_organism = lambda g, t: (None, None)   # force the xref path
    res = ro.resolve("EGFR", 7955)                            # 7955 = zebrafish

    by_common = {o["common"]: o for o in res["orthologs"]}
    assert by_common["human"]["uniprot"] == "P00533"
    assert by_common["mouse"]["uniprot"] == "Q01279"
    assert by_common["rat"]["uniprot"] == "Q9QX70"
    zf = [o for o in res["orthologs"] if o["taxid"] == 7955]
    assert len(zf) == 2, "should keep both fish orthologs"
    assert all(o["is_target"] for o in zf), "fish should be flagged as target"
    assert {o["uniprot"] for o in zf} == {"DARA1", "DARB1"}
    # human flagged not-target, mouse not-target
    assert by_common["human"]["is_target"] is False
    print("resolve() end-to-end OK:",
          [(o["common"], o["uniprot"], o["orthology_type"]) for o in res["orthologs"]])


def test_dedup_same_accession():
    dup = {"data": [{"id": "ENSG1", "homologies": [
        {"species": "danio_rerio", "id": "ENSDARG_A", "protein_id": "P_A", "type": "ortholog_one2many"},
        {"species": "danio_rerio", "id": "ENSDARG_B", "protein_id": "P_B", "type": "ortholog_one2many"},
    ]}]}

    def fake_get(path, params=None):
        if path.startswith("/homology/symbol/human/GENEX"):
            return dup
        if path.startswith("/xrefs/id/"):
            return [{"primary_id": "SAME"}]          # both genes map to the same accession
        return None

    ro.ensembl_get = fake_get
    ro.uniprot_by_gene_organism = lambda g, t: (None, None)
    res = ro.resolve("GENEX", 7955)
    zf = [o for o in res["orthologs"] if o["species"] == "danio_rerio"]
    assert len(zf) == 1, f"duplicate accession should collapse to one, got {len(zf)}"
    print("dedup OK: two genes -> same accession collapsed to one")


def test_partial_upgrade():
    XREF_TREMBL = {"P_RAT": "SHORT1", "ENSG1": "FULL_H"}

    def fake_get(path, params=None):
        if path.startswith("/homology/symbol/human/GENEY"):
            return {"data": [{"id": "ENSG1", "homologies": [
                {"species": "rattus_norvegicus", "id": "ENSRNOG1", "protein_id": "P_RAT",
                 "type": "ortholog_one2one"}]}]}
        if path.startswith("/xrefs/id/"):
            eid = path.rsplit("/", 1)[-1]
            db = dict(params or []).get("external_db")
            if db == "Uniprot/SPTREMBL":                 # all unreviewed
                return [{"primary_id": XREF_TREMBL.get(eid)}] if eid in XREF_TREMBL else []
            return []
        return None

    ro.ensembl_get = fake_get
    ro.uniprot_by_gene_organism = lambda g, t: (None, None)
    # SHORT1 is 300 aa, and a 1200 aa complete rat entry exists; human's is already full-length
    ro.uniprot_length = lambda acc: {"SHORT1": 300, "FULL_H": 1200, "RAT_FULL": 1200}.get(acc, 1200)
    ro.longest_uniprot_by_gene = lambda g, t, reviewed_only=False: (
        ("RAT_FULL", "unreviewed", 1200) if t == 10116 else ("FULL_H", "unreviewed", 1200))

    res = ro.resolve("GENEY", 10116)
    rat = [o for o in res["orthologs"] if o["species"] == "rattus_norvegicus"][0]
    assert rat["uniprot"] == "RAT_FULL", rat
    assert "partial" in rat["note"]
    human = [o for o in res["orthologs"] if o["species"] == "homo_sapiens"][0]
    assert human["uniprot"] == "FULL_H" and human["note"] == "", "full-length must be untouched"
    print("partial upgrade OK: rat SHORT1(300aa) -> RAT_FULL(1200aa); full human left alone")




def test_xref_selection_is_deterministic():
    """Ensembl does not guarantee xref order; the same response must give the same accession."""
    resp = [{"primary_id": "Q9UM11"}, {"primary_id": "P12830"}]
    assert ro.parse_xref_uniprot(resp) == ro.parse_xref_uniprot(list(reversed(resp)))
    assert ro.parse_xref_uniprot([{"primary_id": "B"}, {"primary_id": "A"}]) == "A"
    assert ro.parse_xref_uniprot([]) is None and ro.parse_xref_uniprot(None) is None
    print("xref determinism OK: order of the response no longer changes the accession")


def test_reviewed_partial_is_upgraded_within_swissprot():
    """A short REVIEWED entry must be upgraded, but only to another reviewed one."""
    calls = []
    real_len, real_longest = ro.uniprot_length, ro.longest_uniprot_by_gene
    ro.uniprot_length = lambda acc: {"SHORT_REV": 496, "FULL_REV": 882}.get(acc, 882)

    def fake_longest(gene, taxid, reviewed_only=False):
        calls.append(reviewed_only)
        return ("FULL_REV", "reviewed", 882) if reviewed_only else ("FULL_TREMBL", "unreviewed", 900)
    ro.longest_uniprot_by_gene = fake_longest
    try:
        acc, status, note = ro.upgrade_partial("SHORT_REV", "reviewed", "CDH1", 9606)
        assert acc == "FULL_REV" and status == "reviewed", (acc, status)
        assert calls == [True], "a reviewed entry must only be compared against reviewed ones"
        assert "496aa -> complete 882aa" in note, note
        # a full-length entry is left alone
        assert ro.upgrade_partial("FULL_REV", "reviewed", "CDH1", 9606)[0] == "FULL_REV"
    finally:
        ro.uniprot_length, ro.longest_uniprot_by_gene = real_len, real_longest
    print("reviewed-partial upgrade OK: 496aa -> 882aa, stays within Swiss-Prot")


if __name__ == "__main__":
    test_parse_homologies_condensed()
    test_parse_homologies_nested()
    test_parse_xref()
    test_resolve_end_to_end()
    test_dedup_same_accession()
    test_partial_upgrade()
    test_xref_selection_is_deterministic()
    test_reviewed_partial_is_upgraded_within_swissprot()
    print("\nALL RESOLVE-ORTHOLOGS TESTS PASSED")
