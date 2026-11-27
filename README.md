# AbTaxa

**Predict cross-species antibody cross-reactivity for membrane proteins, by membrane region.**

AbTaxa helps you decide whether a commercial antibody — raised against human, mouse, rat, or
rabbit — is likely to cross-react with the same protein in *your* species. It compares
per-region sequence identity (extracellular / transmembrane / intracellular) between your target
species and each antibody-host species, because whether an antibody binds depends on how
conserved the specific region it targets is.

Give it a gene and a target species; AbTaxa finds the orthologs, works out the membrane
topology, and reports how conserved each region is — and, if you have a specific antibody, how
conserved that antibody's exact epitope is.

> **AbTaxa is a prioritisation / screening tool, not a guarantee of binding.** It reports
> sequence identities for you to interpret. Conformational epitopes, glycosylation, and single
> hot-spot residues can change reactivity independently of overall identity.

**No installation needed for single cases:** the [web version](https://yuingan.github.io/abtaxa/)
runs this same code in your browser (see [Web version](#web-version)).

---

## What it does

- Resolves orthologs of your gene across human, mouse, rat, rabbit **and your target species**
  (Ensembl Compara), preferring complete over partial sequences.
- Determines each model protein's membrane topology, best evidence first: its own curated
  UniProt annotation; otherwise the curated topology of the first species that has one, in the
  order human → mouse → rat → rabbit, carried across through the alignment; "soluble" when no
  species' entry annotates a membrane segment and at least one of them is a reviewed
  (Swiss-Prot) entry; TMbed prediction only when none of these applies.
- Reports **per-region % identity** (extracellular / transmembrane / intracellular / overall)
  between the target and each model species.
- Optionally, for a specific antibody: locates its immunogen (paste the peptide **or** give
  coordinates) and reports the epitope's conservation across species, with a difference-highlighted
  alignment, per-residue substitution severity, and N-glycosylation changes.
- Outputs a terminal table, a JSON file, and a self-contained HTML report.

## Web version

<https://yuingan.github.io/abtaxa/> runs AbTaxa's own Python files (`abtaxa.py`,
`resolve_orthologs.py`, `structure_compare.py`, copied unchanged) inside the browser with
[Pyodide](https://pyodide.org). Nothing is installed and no AbTaxa server is involved:
queries go from the browser straight to Ensembl and UniProt, as they do from the command line.
The first visit downloads about 15 MB (Python, numpy, Biopython), which the browser then caches.

What it leaves out:

- **TMbed.** It is only needed when no species in the run has a curated topology and no
  reviewed entry rules out a membrane segment; a species without its own curated topology
  borrows one exactly as on the command line. In that remaining case the web version reports
  the regions as unknown (the protein is not assumed to be soluble); use the command line.
- **Batch mode.** One case at a time.

Everything else (ortholog choice, the curated-entry preference and its checks, region and
epitope identities, substitutions, glycosylation sequons, the numbering check, the
region-by-region alignments) is the same code as the command line and gives the same numbers.
The site is rebuilt from this repository on every push to `main`, after the tests pass
(`.github/workflows/pages.yml`).

## Installation

AbTaxa depends on a few bioinformatics packages — notably **TMbed** (a ProtT5-based topology
predictor). The reliable route is conda:

```bash
# 1. create the environment
conda env create -f environment.yml
conda activate abtaxa

# 2. install TMbed separately, with --no-deps so it does not override the pinned transformers
pip install --no-deps git+https://github.com/BernhoferM/TMbed.git

# 3. install AbTaxa itself (provides the `abtaxa` command)
pip install .
```

**Notes on TMbed**

- It downloads the ProtT5 model (~2–3 GB) on first run, so the first prediction of a session is
  slow (~20–40 s on CPU); later ones are faster.
- It requires `transformers==4.30.2` (pinned in `environment.yml`); newer versions break it. The
  `--no-deps` install in step 2 stops TMbed pulling an incompatible version.
- On **Windows, use WSL2** — the bioinformatics dependencies are Linux/macOS only — or the
  web version.
- **TMbed is optional.** Without it (`pip install .` alone), AbTaxa still runs and gives the
  same numbers whenever any species has a curated topology to borrow, or a reviewed entry shows
  the protein is soluble. Only when neither applies are the regions reported as unknown. This
  is what the web version does.

## Quick start

```bash
# CD3E in dog (NCBI taxid 9615), with the immunogen peptide of the CD3-12 antibody
abtaxa --gene CD3E --target-taxon 9615 --epitope-seq ERPPPVPNPDYEP --html dog_cd3e.html
```

The region table from that run (version 0.7.0, November 2026):

```
Sequence identity to target, by membrane region (%):

  model       extracellular  transmembrane  intracellular   overall
  human                48.4           65.4           83.6      63.8
  mouse                48.8           65.4           87.3      61.9
  rat                  59.5           66.7           82.4      62.2
  rabbit               45.4           71.4           86.0      61.7
```

The epitope itself (`ERPPPVPNPDYEP`, intracellular) is 92.3% identical to dog in all four
species, with one conservative change (R→K).

Extracellular identity is what matters for surface / flow-cytometry / IHC antibodies;
intracellular for western-blot / lysate antibodies.

## Checking a specific antibody's epitope

Find the **Immunogen** field on the antibody's datasheet and give it to AbTaxa. Two ways:

**Paste the immunogen peptide** (easiest — AbTaxa finds where it is *and which species it came
from*, so you don't have to specify):

```bash
abtaxa --gene EGFR --target-taxon 10181 --epitope-seq "CQRLQKDVEKMGVGWKHLGSSKNLKELPMNNL"
```

**Or give coordinates** (a residue range on a stated species — coordinates only mean something
relative to a specific sequence):

```bash
abtaxa --gene EGFR --target-taxon 10181 --epitope-coor 1005-1036 --epitope-ref human
```

The epitope section reports each species' coordinates, length, indel count, and identity to the
target over that region, and shows an alignment with the residues that differ from the target
highlighted — so you can see exactly which residues might weaken binding.

Two safeguards worth knowing about. A pasted peptide must match a model species at **≥70%** or
AbTaxa refuses to place it and reports no epitope — without that floor, a mistyped or
wrong-species peptide would still land at its best-scoring window and look authoritative. And
`--epitope-coor` coordinates are validated against the reference length: residue numbers only
mean something relative to one sequence, and datasheets sometimes number the **mature** protein
(signal peptide removed) while UniProt numbers from Met1. If the two conventions differ your
range is silently shifted, so prefer `--epitope-seq` when the datasheet gives a peptide.

## Command-line options

| Option | Description |
|---|---|
| `--gene` | Human gene symbol, e.g. `EGFR`. **Required.** |
| `--target-taxon` | NCBI taxid of your target species, e.g. `7955`. **Required.** |
| `--epitope-seq` | Immunogen peptide sequence from the datasheet; auto-located. |
| `--epitope-coor` | Immunogen residue range, e.g. `300-350`, on the reference species. |
| `--epitope-ref` | Species the `--epitope-coor` numbers refer to (default: `human`). |
| `--html` | Write a self-contained HTML report to this path. |
| `--json` | Write the full result as JSON to this path. |
| `--batch` | CSV of cases to run in one go (see **Batch mode**). |
| `--batch-out` | Where to write batch results (default `abtaxa_batch.csv`). |
| `--timeout` | Ensembl timeout in seconds; raise for large gene families (default 120). |
| `--no-cache` | Ignore the on-disk cache of sequences and topology. |
| `--version` | Print the AbTaxa version. |

Common target taxids: `7955` zebrafish · `10181` naked mole-rat · `9031` chicken · `9823` pig ·
`9913` cow · `9615` dog · `8364` *Xenopus tropicalis*.

## Batch mode

To run many antibodies/species at once — for a validation set, or just a panel of targets —
put the cases in a CSV and use `--batch`:

```csv
case,gene,target_taxon,epitope_coor,epitope_ref,known_outcome
cetuximab_cyno,EGFR,9541,313-482,human,cross-reactive
cetuximab_mouse,EGFR,10090,313-482,human,not-cross-reactive
zo1_1a12_nmr,TJP1,10181,,human,cross-reactive
```

```bash
abtaxa --batch cases.csv --batch-out results.csv
```

Required columns are `gene` and `target_taxon`; `case`, `epitope_coor`, `epitope_seq` and
`epitope_ref` are optional. **Any other columns you add (such as `known_outcome`) are copied
through to the results**, so the output CSV can be analysed directly against ground truth
without a join. Results are tidy — one row per (case, model species) — with the per-region
identities, the epitope measures, and coverage.

Three columns carry what the terminal report would otherwise tell you:

- `note` — anything that qualifies the numbers: a model skipped as a partial sequence, an
  immunogen peptide that matched only approximately or was refused, a range cut to the
  protein's length, low alignment coverage.
- `epitope_residues` and `epitope_alt_numbering` — the residues found at the range you gave,
  and the range and residues under the *other* numbering convention (see **Coordinate
  numbering**). Compare both with the datasheet before trusting an epitope row.
- `error` — why a case produced no result. One failing case (an unresolved ortholog, bad
  coordinates) is recorded and the batch carries on; rows are written as each case finishes,
  so an interrupted run keeps what it has done.

Two worked examples are included: therapeutic antibodies with structurally mapped epitopes in
[`validation/validation_cases.csv`](validation/validation_cases.csv), and everyday research
antibodies with datasheet immunogens in
[`validation/common_antibodies.csv`](validation/common_antibodies.csv).

## Caching

Sequences and topology never change for a given accession, so AbTaxa caches them under
`~/.abtaxa_cache` (override with `ABTAXA_CACHE`, disable with `--no-cache`). The first run of
a protein pays the TMbed cost (~20–40 s); repeat and batch runs are near-instant, which also
keeps load off the public Ensembl and UniProt APIs. What AbTaxa *derives* from a UniProt entry
(topology, how the protein is anchored, signal peptide) is cached per AbTaxa version, so an
upgrade that reads entries better re-reads them once; sequences and TMbed predictions are kept.

## Beyond percent identity

Percent identity treats every mismatch alike, but an antibody does not — Leu→Ile is nearly
invisible to it, while Arg→Glu reverses the charge. Two epitopes at the same identity can
differ in whether binding survives, so AbTaxa reports three further measures alongside it.

**Similarity** counts positions that are identical *or* a conservative substitution by
BLOSUM62, using BLAST's "positives" convention of a strictly positive score. (Score 0 is
excluded deliberately: Arg→Glu scores 0 yet reverses the charge.)

**Substitution severity** uses the Grantham (1974) physicochemical distance over composition,
polarity and molecular volume. Values are computed from Grantham's own formula rather than
transcribed from a table, and the unit tests check them against his published anchors
(Leu→Ile ≈ 5, Cys→Trp ≈ 215). Substitutions are banded as *conservative* (≤ 50),
*moderate* (≤ 100) and *radical* (> 100), and each one in the epitope is listed with its
position, change and class — so a charge reversal at a contact residue is visible rather than
averaged away.

**N-glycosylation sequons** (N-X-[S/T], X ≠ P) are compared between the immunogen species and
the target across the epitope. A glycan is bulky and can physically block antibody access, so
a sequon *gained* in the target can abolish binding even at 95% identity — a failure mode
sequence identity cannot see at all.

None of these makes the output a binding prediction; they make the reasons behind a percentage
visible.

## Coordinate numbering

Residue numbers only mean something relative to one sequence. UniProt counts from Met1, while
datasheets sometimes number the **mature** protein — after the signal peptide, and after the
propeptide where there is one. The shift is 24 residues for EGFR and 154 for E-cadherin, enough
to point at the wrong region entirely.

A real example: a widely used anti-CD3 immunogen is printed on datasheets as "aa 156-168
(ERPPPVPNPDYEP) from human CD3 epsilon (UniProt P07766)". In the UniProt sequence that peptide
is at **178-190**; residues 156-168 are `KAKAKPVTRGAGA`. The datasheet is counting from the
mature protein (22-residue signal peptide removed). Taken at face value the range gives 92%
identity to mouse for the wrong thirteen residues, instead of 100% for the right ones.

When you pass `--epitope-coor`, AbTaxa therefore:

- validates the range against the reference length (an out-of-range start is an error, an
  overhanging end is cut and noted);
- prints the residues actually found at that range;
- reads the mature-chain start from UniProt (signal peptide, propeptide, or a removed initiator
  methionine) and prints the range *and residues* under the mature-numbering reading too — or
  says that reading is impossible when it would run past the end of the protein.

Compare both sets of residues with the datasheet and re-run with the second range if that is
the one that matches. Passing the immunogen peptide with `--epitope-seq` avoids the ambiguity
entirely and is the better route whenever the datasheet gives a sequence — including when the
peptide carries a non-native linker residue such as a terminal cysteine, which simply shows up
as a match slightly under 100%.

The epitope reference must be one of the species compared in the run. If it is not (it is the
target itself, or its ortholog was skipped as a partial sequence), AbTaxa stops and lists the
species that are available rather than reading your coordinates against a different sequence.

## How to read the output

- **Region table** — % sequence identity between each antibody-host species and your target,
  split by membrane region. Higher means more conserved, so an antibody targeting that region is
  more likely to cross-react. Each cell is identity over aligned residues; where some residues
  are present in only one of the two proteins, a second figure ("counting unmatched") counts
  them as differences, so a protein of very different length shows up as different. The
  `topology` line under each species says where its regions came from (its own UniProt entry,
  or borrowed from another species).
- **Target sequence by region** — the whole target, in its own numbering, mapped onto the
  reference topology through the alignment: which stretch is extracellular, transmembrane,
  intracellular or signal peptide, and which have no counterpart in the reference.
- **Residues present in only one protein** — every stretch with no counterpart in the other
  protein: who has it, its residues, its length, whether it is N-terminal, C-terminal or
  internal, and its region. A stretch inside one region is placed in it; one at a terminus or a
  region boundary is not placed. Example: Ensembl's naked mole-rat THY1 is 256 aa against
  161 aa in the other species, and its first 95 residues have no counterpart in any of them.
  That may be real biology or a gene-model error; AbTaxa shows it rather than deciding.
- **Epitope table** (if an epitope is given) — the same, restricted to the antibody's actual
  immunogen region, with per-species coordinates, length, and an **indel count**. Identity is
  computed over columns present in both species, so an insertion or deletion is excluded from
  the percentage rather than counted against it — a deletion inside an epitope can therefore
  sit behind a high identity. Read the indel count alongside the percentage: an indel reshapes
  the local backbone and is usually more disruptive to binding than a substitution.
- **Alignment** — the epitope residues across species, with differences from the target
  highlighted.
- **Membrane regions, residue by residue** (HTML and JSON) — separate Extracellular,
  Transmembrane and Intracellular sections. Each segment is defined on the human topology (or
  the first model species with a topology) and aligned across every species, like the epitope:
  per-species coordinates, length, indel count and identity to the target, then the residues
  with differences from the target highlighted. A multi-pass protein gets one block per loop
  or helix. The reference is the species whose topology is borrowed (the first curated one in
  the order human → mouse → rat → rabbit). Identities here are measured on the
  reference-defined segments, so they can differ slightly from the region table, where a
  species with its own curated topology uses that one. GPI-anchored proteins show only an
  extracellular section; soluble proteins have none.
- **Similarity / radical substitutions** — see **Beyond percent identity** above. A large gap
  between identity and similarity means the differences are mostly conservative; a non-zero
  radical count points at the specific changes most likely to matter.
- **Coverage** — the fraction of the reference protein that actually aligned. Identity is
  computed over aligned positions only, so a partial or divergent sequence can show high
  identity over a small part of the protein; AbTaxa flags coverage below 80%.

There is deliberately **no verdict** — the percentages are yours to interpret. A common vendor
rule of thumb is that ~85% identity in the relevant region suggests likely cross-reactivity, but
this threshold is not calibrated here.

## Limitations

- Sequence identity is a *proxy* for epitope conservation. It cannot see conformational epitopes
  (residues far apart in sequence but adjacent in 3D), differential glycosylation, or single
  hot-spot residues that dominate binding.
- The ~85% rule of thumb is borrowed, not validated in this tool.
- Orthologs come from Ensembl Compara; species outside Ensembl's vertebrate set will not resolve.
- Requires internet access (queries Ensembl and UniProt).
- Epitope identity is reference-anchored: each species is projected onto the reference's epitope
  columns rather than aligned directly to the target, which is what makes the columns comparable
  when species carry different indels.
- For soluble (non-membrane) proteins the extracellular/transmembrane/intracellular split does
  not apply; AbTaxa detects this and reports whole-protein identity instead.
- **GPI-anchored proteins** (CD90/THY1, CD59, CD73, PrP and other common surface targets) are
  tethered to the outer leaflet by a lipid rather than a transmembrane helix, so a predictor
  finds no helix and would call them soluble. Where UniProt annotates the GPI anchor, AbTaxa
  reports the mature chain as entirely extracellular and labels the source
  `UniProt (GPI-anchored)`. The cleaved C-terminal GPI attachment signal is excluded.
- **Curated entries are preferred.** Compara's cross-reference sometimes points at an
  unreviewed TrEMBL model even where a curated Swiss-Prot entry exists, and those models can
  be wrong exactly where it matters: dog CD3E `A0A8I3MGF0` has a mispredicted C-terminus that
  drops the conserved `ERPPPVPNPDYEP` epitope (30.8% identity) where the reviewed `P27597`
  retains it (92.3%). AbTaxa prefers a reviewed entry, but only after checking it against
  the other species — see below.
- **Every sequence found by gene symbol is checked.** UniProt's symbol search matches
  synonyms: asking for human `CDH1` returns `Q9UM11`, which is FZR1, not E-cadherin. Any step
  that could pull in a same-named protein (the curated-entry preference and the partial-length
  rescue) first checks that the candidate resembles the other species in the run, and refuses
  it otherwise. A sequence resembling none of the others is flagged in the notes for manual
  checking rather than silently used.
- **No ortholog found**: if Ensembl Compara returns no ortholog for the target species,
  AbTaxa stops and says so. It deliberately does not fall back to matching on gene symbol,
  because a gene with the same name is not necessarily the ortholog and a wrong target
  sequence would silently distort every identity in the run. Check such cases by hand.
- **Genomes outside Compara.** Some genomes are in Ensembl but not in Compara's gene trees, so
  no gene has a Compara ortholog there. In November 2026 these included the Damara mole-rat
  (*Fukomys damarensis*, 885580) and the paternal naked mole-rat assembly (the maternal one,
  which AbTaxa uses for 10181, is included). AbTaxa says so instead of giving a per-gene
  message. To check a species: `https://rest.ensembl.org/info/genomes/taxonomy/<taxid>?content-type=application/json`
  and look at `has_peptide_compara`.

- **Withdrawn UniProt entries.** UniProt deletes unreviewed entries that are not part of a
  reference proteome, and Ensembl's cross-references still point at some of them (dog TP53 ->
  `A0A8C0RLV2`, naked mole-rat THY1 -> `A0A0P6J0R8`). AbTaxa then takes the sequence Ensembl
  holds for that same Compara ortholog (labelled `Ensembl:<protein ID>` in every report), and
  the usual Swiss-Prot step still replaces it with a reviewed entry when one exists and
  resembles the other species at least as well (dog TP53 -> `Q29537`). Nothing is matched on
  gene name alone. Ensembl sequences have no UniProt topology, so a model species taken from
  Ensembl borrows one like any other species without its own.
- **Borrowed topology** assumes orthologs share their membrane topology, which is the norm. The
  borrowed boundaries land where the alignment puts them (mouse CD3E's borrowed helix falls at
  109-134, exactly where mouse's own UniProt entry puts it), but a species-specific change in
  topology would not be seen.
- **Each pair is aligned twice and the better alignment is kept.** One alignment leaves
  overhanging ends free, the other runs end to end, and neither is right every time:
  - *Free ends* suit a protein that simply carries an extension. Naked mole-rat THY1 from
    Ensembl starts with ~95 proline/serine-rich residues that no other species has (and no
    initiator Met); with free ends they stay one N-terminal stretch. Aligned end to end, the
    stretch is broken up and 12 of its residues are paired with the mouse signal peptide.
  - *End to end* suits an insertion close to an end. Dog NCAM1 has 254 more residues than human
    NCAM1 inside the cytoplasmic tail, 40 residues before the C-terminus. With free ends alone,
    stopping early cost nothing while opening that gap cost more than the last 40 residues
    were worth, so they were paired with the start of the insertion and human-dog
    intracellular identity read 71.6% where it is 98.3% (versions 0.9.0-0.9.2).

  Both candidates are scored by one yardstick: BLOSUM62 over the aligned pairs, overhanging
  ends free, and an internal gap charged the usual penalty (11 to open, 1 to extend) but
  never more than 30, so a long insertion counts as one event however long it is. The higher
  score is kept, a tie keeps the free-ends alignment, and the JSON records the choice per
  comparison (`alignment_ends`: `free` or `end to end`). Most pairs give the same alignment
  either way.

  Checked against Ensembl Compara's own alignments of the same sequences, human against eight
  mammals for 36 antibody targets (278 comparisons, November 2026): the two pair 98.9% of
  residues the same way in the median comparison, and where the human protein has curated
  regions (249 comparisons) region identity is within 5 percentage points in 245 of them for
  the extracellular region, 227 for the transmembrane region and 221 for the intracellular
  region. Compara aligns the whole gene family at once; AbTaxa aligns two proteins, because
  it often compares a curated UniProt entry rather than Ensembl's own translation and has to
  place immunogen peptides. Short, poorly conserved stretches (a transmembrane helix, a
  cytoplasmic tail) can be lined up in several equally good ways, and there the two disagree
  most; neither figure is exact there.

  What remains: an insertion just inside a terminus can still be reported as a terminal
  overhang when the residues beyond it are few or poorly conserved. Rat CD3E A0A8I6AAG4 has 74
  extra residues after residue 16 of its signal peptide; they are reported as rat 1-74,
  N-terminal, and the dog signal peptide is aligned to the inserted stretch. The mature
  regions (extracellular, transmembrane, intracellular) align the same either way; only the
  signal peptide, and so the overall figure, moves.
- **A gene symbol that is also another gene's synonym.** Ensembl's symbol lookup matches
  synonyms, so one symbol can come back as two human genes: `CDH1` returns CDH1 (E-cadherin)
  and FZR1, which lists CDH1 as a synonym. AbTaxa keeps the gene whose official symbol is
  the one you typed, sets the other aside and says so in the notes. If a symbol is no gene's
  official symbol and a synonym of several, AbTaxa stops and lists them: use the official
  symbol of the one you mean. (Up to 0.9.2 the orthologs of both genes were pooled, and the
  right protein was reached only because FZR1 is much shorter than E-cadherin.)

## Components

- `abtaxa.py` — the driver (ortholog resolution → topology → per-region comparison → report).
- `resolve_orthologs.py` — Ensembl Compara ortholog resolution to UniProt accessions.
- `structure_compare.py` — sequence-based per-region comparison; also a standalone
  structure-based comparison (lDDT/RMSD) for a pair that both have AlphaFold structures
  (needs the optional `tmtools`).
- `web/` — the web version: `index.html` (the page), `worker.js` (loads Pyodide and runs the
  Python files above), `build_site.py` (assembles the site, copying the Python files as-is).

## Development

```bash
python test_abtaxa.py
python test_resolve_orthologs.py
python test_structure_compare.py
```

The tests stub every Ensembl/UniProt call, so they need no network and no TMbed; GitHub Actions
runs them on every push (`.github/workflows/tests.yml`).

To try the web version locally before pushing:

```bash
python web/build_site.py               # writes _site/
python -m http.server -d _site 8000    # then open http://localhost:8000
```

## Citing

If AbTaxa is useful in your work, please cite the repository
(https://github.com/yuingan/abtaxa).

## License

MIT — see [LICENSE](LICENSE).
