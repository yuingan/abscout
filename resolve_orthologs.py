#!/usr/bin/env python3
"""
resolve_orthologs.py
====================
Given a human gene symbol and a TARGET species, find its orthologs in the
antibody-rich model species (human, mouse, rat, rabbit) AND the target species
via Ensembl Compara, and resolve each to a UniProt accession -- so the results
feed straight into structure_compare.py / TMbed.

Per run (Ensembl REST + UniProt REST, all lightweight API calls):
  1. /homology/symbol/human/<GENE>   Compara orthologues (condensed)  -> ortholog gene IDs per species
  2. /xrefs/id/<ensembl_id>          UniProt cross-reference          -> UniProt accession per ortholog
  3. UniProt gene+organism search    fallback if no xref              -> UniProt accession

Handles one-to-many orthology (e.g. duplicated genes in fish) by returning ALL
orthologs found per species, not just one.

Deps: requests  (already in your env)

USAGE
  python resolve_orthologs.py --gene EGFR --target-taxon 7955 --json egfr_orthologs.json
  # taxids: 7955 zebrafish | 9031 chicken | 9823 pig | 9913 cow | 9615 dog | 8364 X.tropicalis

NOTE: built without live access to the Ensembl/UniProt APIs, so the parsing is
unit-tested on mock responses (see test_resolve_orthologs.py) but the live calls
are untested -- run it once on EGFR and eyeball the accessions before a batch.
"""
import argparse
import json
import sys
import time

import requests

ENSEMBL = "https://rest.ensembl.org"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"

# Ensembl production name -> (common name, NCBI taxid)
MODEL_SPECIES = {
    "homo_sapiens": ("human", 9606),
    "mus_musculus": ("mouse", 10090),
    "rattus_norvegicus": ("rat", 10116),
    "oryctolagus_cuniculus": ("rabbit", 9986),
}
NONHUMAN_MODEL_TAXIDS = [10090, 10116, 9986]
COMMON_NAME = {k: v[0] for k, v in MODEL_SPECIES.items()}


# ----------------------------------------------------------------------------------
# HTTP
# ----------------------------------------------------------------------------------
ENSEMBL_TIMEOUT = 120          # seconds; raise via abscout --timeout for slow gene families


def ensembl_get(path, params=None):
    # JSON is requested with Ensembl's content-type URL parameter rather than a Content-Type
    # header. Ensembl documents both; the parameter keeps this a plain GET, which a web
    # browser can send to another site without a CORS pre-flight (the web version of
    # AbScout runs this same function inside the browser).
    if isinstance(params, dict):
        params = list(params.items())
    params = list(params or []) + [("content-type", "application/json")]
    headers = {"Accept": "application/json"}
    for attempt in range(1, 4):
        try:
            # escalate the timeout each attempt -- big gene families (TJP1, OCLN) are slow
            r = requests.get(ENSEMBL + path, params=params, headers=headers,
                             timeout=ENSEMBL_TIMEOUT * attempt)
            if r.status_code == 200:
                return r.json()
            if r.status_code in (400, 404):
                return None
            time.sleep(1.5 * attempt)          # 429/5xx -> back off
        except requests.RequestException as exc:
            sys.stderr.write(f"  [warn] Ensembl request error ({attempt}/3): {exc}\n")
            time.sleep(1.5 * attempt)
    return None


def ensembl_protein_seq(protein_id):
    """The protein sequence Ensembl holds for a translation (ENS...P...): the sequence Compara
    itself used when it called the ortholog. None if Ensembl gave nothing back."""
    resp = ensembl_get(f"/sequence/id/{protein_id}", [("type", "protein")])
    seq = resp.get("seq", "") if isinstance(resp, dict) else ""
    return seq.replace("*", "").strip() or None


def compara_genomes(taxid):
    """
    The Ensembl genomes for a taxon as [(display_name, name, in_compara)], [] when Ensembl
    has none, or None when Ensembl could not be asked. A genome that is in Ensembl but not in
    Compara's gene trees (has_peptide_compara = 0) has no Compara orthologs for any gene --
    the Damara mole-rat, for one, and the paternal naked mole-rat assembly.
    """
    resp = ensembl_get(f"/info/genomes/taxonomy/{taxid}")
    if resp is None:
        return None
    return parse_genomes(resp)


def parse_genomes(resp):
    out = []
    for g in resp if isinstance(resp, list) else []:
        flag = str(g.get("has_peptide_compara", "0")).strip().lower()
        out.append((g.get("display_name") or g.get("name", "?"), g.get("name", "?"),
                    flag not in ("0", "", "none", "false")))
    return out


# ----------------------------------------------------------------------------------
# Parsing (pure functions -- unit-tested)
# ----------------------------------------------------------------------------------
def parse_homologies(resp):
    """
    Parse an Ensembl /homology condensed (or full) response.
    Returns (source_gene_id, [ {species, gene_id, protein_id, type}, ... ]).
    Tolerant of both the flat 'condensed' layout and the nested 'target' layout.
    """
    src_gene = None
    out = []
    for entry in (resp or {}).get("data", []):
        src_gene = src_gene or entry.get("id")
        for h in entry.get("homologies", []):
            tgt = h.get("target") if isinstance(h.get("target"), dict) else h
            out.append({
                "species": h.get("species") or tgt.get("species"),
                "gene_id": tgt.get("id") or h.get("id"),
                "protein_id": tgt.get("protein_id") or h.get("protein_id"),
                "type": h.get("type"),
            })
    return src_gene, out


def parse_xref_uniprot(resp):
    """
    Return one UniProt accession from an Ensembl /xrefs/id response list.

    Ensembl can list several accessions for one gene and does not guarantee their order,
    so taking whichever came first made the same query resolve differently between runs
    (human CDH1 gave both the full-length entry and a 496 aa one). Sort first so a run is
    reproducible; `upgrade_partial` then fixes the choice if a longer entry exists.
    """
    accs = sorted({x["primary_id"] for x in (resp or []) if x.get("primary_id")})
    return accs[0] if accs else None


# ----------------------------------------------------------------------------------
# Ensembl id -> UniProt accession (xref, then UniProt gene+organism fallback)
# ----------------------------------------------------------------------------------
def ensembl_to_uniprot(ensembl_id):
    if not ensembl_id:
        return None, None
    for db, status in (("Uniprot/SWISSPROT", "reviewed"), ("Uniprot/SPTREMBL", "unreviewed")):
        resp = ensembl_get(f"/xrefs/id/{ensembl_id}",
                           [("external_db", db), ("all_levels", "1")])
        acc = parse_xref_uniprot(resp)
        if acc:
            return acc, status
    return None, None


def _uniprot_search(params, tries=4):
    """UniProt search with polite backoff; returns the results list (possibly empty)."""
    for i in range(1, tries + 1):
        try:
            r = requests.get(UNIPROT_SEARCH, params=params, timeout=30)
        except requests.RequestException:
            time.sleep(1.5 * i)
            continue
        if r.status_code == 200:
            return r.json().get("results", [])
        if r.status_code in (429, 500, 502, 503, 504) and i < tries:
            time.sleep(min(float(r.headers.get("Retry-After", 1.5 * i)), 20))
            continue
        return []
    return []


def uniprot_by_gene_organism(gene, taxid):
    if not gene or not taxid:
        return None, None
    results = _uniprot_search({"query": f"(gene_exact:{gene}) AND (organism_id:{taxid})",
                               "fields": "accession,reviewed", "format": "json", "size": "5"})
    results.sort(key=lambda e: 0 if "reviewed" in e.get("entryType", "").lower() else 1)
    if results:
        e = results[0]
        status = "reviewed" if "reviewed" in e.get("entryType", "").lower() else "unreviewed"
        return e.get("primaryAccession"), status
    return None, None


def uniprot_length(acc):
    """Sequence length for an accession (or None)."""
    results = _uniprot_search({"query": f"accession:{acc}", "fields": "length",
                               "format": "json", "size": "1"})
    if results:
        e = results[0]
        return (e.get("sequence") or {}).get("length") or e.get("length")
    return None


def longest_uniprot_by_gene(gene, taxid, reviewed_only=False):
    """
    Longest entry for gene+organism. Returns (acc, status, length).
    With reviewed_only, considers Swiss-Prot entries alone, so a curated accession is never
    traded for an unreviewed one. Ties break on the accession so the result is reproducible.
    """
    q = f"(gene_exact:{gene}) AND (organism_id:{taxid})"
    if reviewed_only:
        q += " AND (reviewed:true)"
    results = _uniprot_search({"query": q, "fields": "accession,length,reviewed",
                               "format": "json", "size": "25"})
    best = None
    for e in results:
        ln = (e.get("sequence") or {}).get("length") or e.get("length") or 0
        rev = "reviewed" in e.get("entryType", "").lower()
        acc = e.get("primaryAccession") or ""
        key = (ln, 1 if rev else 0, acc)
        if best is None or key > best[0]:
            best = (key, acc, "reviewed" if rev else "unreviewed")
    return (best[1], best[2], best[0][0]) if best else (None, None, None)


def upgrade_partial(acc, status, gene, taxid, ratio=1.3):
    """
    If a substantially longer entry exists for gene+organism, prefer it -- this catches
    partial sequences that UniProt does not formally flag as fragments, which otherwise
    drop whole regions and deflate identity.

    Applies to reviewed entries too: Swiss-Prot carries short entries (human CDH1 resolved
    to a 496 aa reviewed accession against an 882 aa protein). A reviewed accession is only
    ever replaced by another reviewed one, so curation is never traded away for length.
    """
    cur = uniprot_length(acc)
    if not cur:
        return acc, status, ""
    best_acc, best_status, best_len = longest_uniprot_by_gene(
        gene, taxid, reviewed_only=(status == "reviewed"))
    if best_acc and best_acc != acc and best_len and best_len >= cur * ratio:
        return best_acc, best_status, f"partial {cur}aa -> complete {best_len}aa"
    return acc, status, ""


def taxid_for(species_name, target_taxid):
    """Model species -> known taxid; any non-model species in the response is the target."""
    info = MODEL_SPECIES.get(species_name)
    return info[1] if info else target_taxid


# ----------------------------------------------------------------------------------
# Orchestration
# ----------------------------------------------------------------------------------
def resolve(gene, target_taxid):
    target_taxa = sorted(set(NONHUMAN_MODEL_TAXIDS) | {target_taxid})
    params = [("type", "orthologues"), ("format", "condensed"), ("sequence", "none")]
    params += [("target_taxon", str(t)) for t in target_taxa]
    resp = ensembl_get(f"/homology/symbol/human/{gene}", params)
    src_gene, homs = parse_homologies(resp)

    result = {"gene": gene, "query_species": "homo_sapiens",
              "human_gene_id": src_gene, "target_taxid": target_taxid, "orthologs": [],
              # False when Ensembl gave no answer at all (request failed, or the symbol is
              # not one it recognises) -- not the same as "no ortholog in that species"
              "homology_query_ok": resp is not None}

    def add(species, taxid, ensembl_gene, protein_id, otype):
        eid = protein_id or ensembl_gene
        acc, status = ensembl_to_uniprot(eid)
        if not acc and ensembl_gene and ensembl_gene != eid:
            acc, status = ensembl_to_uniprot(ensembl_gene)
        if not acc:
            acc, status = uniprot_by_gene_organism(gene, taxid)
        note = ""
        # Partial sequences (often not flagged as 'fragment') silently drop whole regions.
        # Unreviewed entries are the common offenders and are checked here; a reviewed entry
        # can also be short, but verifying every one would cost two extra UniProt calls per
        # ortholog, so that case is caught in abscout.analyze(), which already holds the
        # sequence and the target length and can act only when the length actually looks wrong.
        if acc and status == "unreviewed":
            acc, status, note = upgrade_partial(acc, status, gene, taxid)
            if note:
                sys.stderr.write(f"  [{species}] {note}\n")
        result["orthologs"].append({
            "species": species, "common": COMMON_NAME.get(species, species),
            "taxid": taxid, "ensembl_gene": ensembl_gene, "ensembl_protein": protein_id,
            "uniprot": acc, "uniprot_status": status, "note": note,
            "orthology_type": otype, "is_target": (taxid == target_taxid),
        })

    if src_gene:
        add("homo_sapiens", 9606, src_gene, None, "self")
    for h in homs:
        add(h["species"], taxid_for(h["species"], target_taxid),
            h["gene_id"], h["protein_id"], h["type"])

    # collapse exact (species, UniProt) duplicates -- e.g. paralogs mapping to one accession
    seen, deduped = set(), []
    for o in result["orthologs"]:
        key = (o["species"], o["uniprot"])
        if o["uniprot"] and key in seen:
            continue
        seen.add(key)
        deduped.append(o)
    result["orthologs"] = deduped
    return result


# ----------------------------------------------------------------------------------
# Report + CLI
# ----------------------------------------------------------------------------------
def format_table(result):
    L = [f"Orthologs of {result['gene']} (target taxid {result['target_taxid']}):", ""]
    L.append(f"  {'common':<14}{'ensembl_gene':<24}{'UniProt':<13}{'status':<12}"
             f"{'orthology':<20}{'tgt'}")
    for o in result["orthologs"]:
        acc = o["uniprot"] or "(none)"
        L.append(f"  {(o['common'] or '?'):<14}{(o['ensembl_gene'] or '-'):<24}{acc:<13}"
                 f"{(o['uniprot_status'] or '-'):<12}{(o['orthology_type'] or '-'):<20}"
                 f"{'*' if o['is_target'] else ''}")
    missing = [o for o in result["orthologs"] if not o["uniprot"]]
    if missing:
        L.append(f"\n  [warn] {len(missing)} ortholog(s) had no UniProt accession -- "
                 "no AlphaFold structure will be fetchable for those.")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description="Resolve gene orthologs to UniProt accessions via Ensembl Compara.")
    ap.add_argument("--gene", required=True, help="human gene symbol, e.g. EGFR")
    ap.add_argument("--target-taxon", required=True, type=int,
                    help="NCBI taxid of the target species (e.g. 7955 zebrafish)")
    ap.add_argument("--json", help="optional path to write the result as JSON")
    args = ap.parse_args()

    result = resolve(args.gene, args.target_taxon)
    print(format_table(result))
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(result, fh, indent=2)
        sys.stderr.write(f"\n[wrote] {args.json}\n")


if __name__ == "__main__":
    main()
