# CLAUDE.md — MMRRC

This project is a hackathon involving the MMRRC (Mutant Mouse Resource &
Research Centers, https://www.mmrrc.org/). We'll flesh out these notes as we go.

## The MMRRC catalog dataset

`data/mmrrc_catalog_data.csv` (gitignored, ~147 MB) was downloaded from
https://www.mmrrc.org/methods/data_download.php. 588,980 rows × 18 columns.

Header note: `RESEARCH_AREAS ` has a trailing space in the file. Strip it on
load with `pl.read_csv(path).rename(str.strip)`.

### Grain: long format, two disjoint kinds of row — *not* one row per strain

`STRAIN/STOCK_ID` is **not** unique: 69,388 distinct strains over 588,980 rows
(~8.5 rows/strain, no nulls). Strain-level fields are repeated verbatim on every
row of a strain.

Every row is one of two disjoint kinds — **zero rows carry both**:

| Kind | Rows | `GENE_SYMBOL` / `MGI_GENE_ACCESSION_ID` | `ALLELE_SYMBOL` / `MGI_ALLELE_ACCESSION_ID` |
|---|---:|---|---|
| gene row | 534,057 | set | null |
| allele row | 44,591 | null | set |
| neither | 10,332 | null | null |

So **a gene is never joined to an allele within a row** — they connect only
through `STRAIN/STOCK_ID`. `MMRRC:000001-UNC` shows the shape exactly: three
gene rows (`Fgg`, `Fgb`, `Fga`) plus one allele row (`Tg(Fga,Fgb,Fgg)1Unc`).

Per strain: mean 7.7 genes (max 206 — `MMRRC:042022-MU`, a lesion spanning 206
genes across 20 chromosomes), and 9,178 strains have no gene row at all.
Alleles per strain: 37,225 have 0, 19,877 have 1, 12,162 have 2, 113 have 3,
10 have 4, one has 6.

**Consequence: any per-strain count taken off the raw frame is inflated.** Use
`catalog.unique(subset=["STRAIN/STOCK_ID"])` for strain-level questions, and
treat the raw frame as the strain↔gene edge list. Beware `n_unique()` on
nullable columns — it counts null as a value, which silently turns "no allele"
into "1 allele". Use `.drop_nulls().n_unique()`.

`(STRAIN/STOCK_ID, ALLELE_SYMBOL, GENE_SYMBOL)` is *nearly* a key — 97 violating
groups covering 229 rows. All of them are "neither" rows; within those groups
only `MUTATION_TYPE` and `CHROMOSOME` ever differ, and 82 rows are byte-identical
duplicates.

### The ENU trap — read before ranking anything

**8,229 chemically-induced strains (`MUTATION_TYPE == "CI"`) carry a mean of 58
genes each and account for 478,364 of the 588,980 rows (81%).** Every other
mutation type averages one to three genes per strain.

An ENU strain's gene list is a set of *candidate variants found by sequencing*,
not deliberate manipulations. Longer genes collect more random hits, so ranking
genes by strain count with `CI` included returns `Ttn`, `Obscn`, `Neb`, `Rsf1`,
`Dst` — a gene-length ranking, not a research-interest one. Excluding `CI` gives
`Hprt1`, `EGFP`, `cre`, `lacZ`, `Ctbp2`. **8,167 of 24,593 genes (33%) appear
only in CI strains** and have never been deliberately targeted.

Default to excluding `CI` for any "which genes matter" question, and say so.

### Reading the nomenclature

Allele and strain symbols are structured MGI/ILAR nomenclature, not free text.
28,406 distinct allele symbols break down as:

| Shape | Count | Meaning | Example |
|---|---:|---|---|
| `Gene<tmN…>` | 20,857 | targeted mutation | `Diaph3<tm1a(KOMP)Mbp>` |
| `Gene<emN…>` | 4,647 | endonuclease-mediated (CRISPR) | `Prim1<em1(IMPC)J>` |
| `Tg(…)` | 1,867 | transgene, free-standing | `Tg(Rarb-EGFP)IT82Gsat` |
| `Gene<other>` | 590 | named/spontaneous allele | `Nr2e1<frc>` |
| `Gene<GtN…>` | 348 | gene trap | `Sigmar1<Gt(OST422756)Lex>` |
| bare | 75 | unstructured | `aspb`, `Et(icre)21733Rdav` |
| `Del(…)` | 22 | deletion spanning a region | `Del(7Cyp2s1-Cyp2f2)2Ding` |

**Angle brackets are this file's rendering of a superscript.** `Diaph3<tm1a(KOMP)Mbp>`
is printed elsewhere as Diaph3^tm1a(KOMP)Mbp. Parse on `<`/`>`, and expect the
gene symbol before the bracket to duplicate `GENE_SYMBOL` — except that allele
rows have `GENE_SYMBOL` null (see grain above), so the symbol prefix is often the
*only* place the gene name appears on that row.

Worked example — `MMRRC:000001-UNC`, designation
`C57BL/6-Tg(Fga,Fgb,Fgg)1Unc/Mmnc`:

| Token | Meaning |
|---|---|
| `C57BL/6` | background inbred strain |
| `-` | joins background to the allele it carries |
| `Tg` | allele class: transgene, inserted at a random site |
| `(Fga,Fgb,Fgg)` | **contents of the insert** — here three fibrinogen genes |
| `1` | serial number for that lab (can be alphanumeric, e.g. `IT82`) |
| `Unc` | ILAR lab registration code (University of North Carolina) |
| `/Mmnc` | holding facility — MMRRC at UNC |

The parenthesised part of a `Tg(…)` is the **payload**, not a locus. Parenthesised
`(KOMP)`, `(IMPC)`, `(EUCOMM)` in a `tm` symbol are the originating project, not a
gene. On KOMP/EUCOMM alleles the `a`/`b`/`c`/`d` suffix on `tm1` marks the
conditional-ready cassette series (`tm1a` knockout-first → `tm1b`/`tm1c`/`tm1d`
after Cre and/or Flp), so `tm1a` and `tm1b` on the same gene are the same project,
not independent alleles.

### Gotchas in the gene and chromosome columns

- `CHROMOSOME` is free text: 50 distinct values for what should be 22, including
  `unknown`, `UN`, `unk`, `N/A`, `Chr 1`, `Chr11:4938754-4948064 bp`, `8q21.13`
  (a *human* cytoband) and `919`. Normalise with uppercase → strip a `CHR`
  prefix → keep `1`–`19`, `X`, `Y`, `MT`, bucket the rest as unmapped.
- Rows on "chromosome" 20, 21 and 22 are **human transgenes** — mice have 19
  autosomes plus X/Y. `SOD1` and `APP` on chr21 are human coordinates.
- **For transgenics, `CHROMOSOME` on a gene row is the gene's normal address,
  not the insertion site.** A `Tg` insert lands at random, so the gene rows
  report where the *endogenous* copy lives while the allele row admits the truth:
  of 2,110 `TG` allele rows, 1,735 say `unknown` and 82 are null — yet 2,192 of
  the 4,001 `TG` gene rows (2,048 strains) assert a numeric chromosome.
  `MMRRC:000001-UNC` is the pattern: three gene rows on chromosome 3 (the mouse
  fibrinogen cluster) plus one allele row, `Tg(Fga,Fgb,Fgg)1Unc`, on `unknown`.
  Any "which chromosome is most manipulated" count inherits this — it measures
  where the copied genes normally sit, not where anything was engineered.
- A transgene's genes are not necessarily from another species. `Tg(Fga,Fgb,Fgg)1Unc`
  is a ~100 kb P1 clone of the mouse's *own* fibrinogen genes, used to build a
  hyperfibrinogenemia model (PMID 11521996) — hence mouse MGI ids on those rows.
  Check the case of the symbol before assuming a humanised model.
- 723 gene symbols are ALL-CAPS (592 with no MGI id): non-mouse transgenes.
- `EGFP`, `cre`, `lacZ`, `tTA` are cassettes, not loci — no MGI id, no
  chromosome. They rank near the top of any non-CI gene ranking.
- 23 gene symbols map to more than one MGI id; 19 genes appear on more than one
  chromosome value.

### Phenotypes: `MPT_IDS` and the Mammalian Phenotype Ontology

`MPT_IDS` holds pipe-separated `label [MP:id]` pairs. 4,604 strains are annotated
(**6.6% of 69,388** — absence means absence of characterisation, not absence of an
effect), 46,860 entries, 5,463 distinct MP terms once repaired. Annotated strains
skew to targeted mutations (2,666 `TM` vs 374 `CI`).

**The labels in this column are damaged; the MP ids are intact. Read the ids and
throw the labels away.** Filed upstream as
[issue #1](https://github.com/gaurav/mmrrc/issues/1).

MMRRC's exporter joins the phenotype names into one comma-separated string,
re-splits that string on `", "`, and zips the pieces against the id list,
truncating to its length:

```python
labels = ", ".join(names).split(", ")[:len(ids)]   # what MMRRC appears to do
```

MP names legitimately contain commas (`decreased CD4-positive, alpha-beta T cell
number`), so each one is torn in two, everything after slides by one, and the
names at the end fall off. **1,450 of 4,604 annotated strains are affected and
2,148 phenotype names are dropped from the file.** Reconstructing the column from
the ids by exactly that procedure reproduces 43,900 of 48,069 label slots (91.3%)
character for character — that is what identifies the mechanism.

All 48,069 ids are valid MP terms, so recovery is a one-liner:

```python
pl.col("MPT_IDS").str.extract_all(r"MP:\d+")   # then explode + join to mp.owl
```

**Do not try to recover terms by matching label text.** That was tried first and
is wrong: it discards correct ids, invents terms from torn fragments, cannot
recover names that were never published, and lost ~1,250 real annotations.

The stored labels are a stale snapshot too (`hypoactivity` → `decreased locomotor
activity`, `retinal degeneration` → `retina degeneration`), which is a second
reason to render every label from `mp.owl`.

Watch for `MP:0002169` "no abnormal phenotype detected" — 442 strains, the second
most common entry. It is a negative result, not a phenotype. The ontology files
it under *normal phenotype*, which is the cheapest way to spot it.

### Parsing `data/mp.owl`

101 MB RDF/XML (gitignored), 15,288 MP terms (457 obsolete) plus merged
PATO/UBERON/GO/CHEBI/CL imports — 126,454 `owl:Class` elements in total.

**No ontology library needed.** A single `xml.etree.ElementTree.iterparse` pass
pulling `owl:Class` elements whose `rdf:about` starts with `…/obo/MP_` takes
~1.5 s — fast enough to run live in a cell. `rdflib`/`pronto`/`owlready2` would
all be slower and add a dependency.

Two traps:

- **Only `el.clear()` the `owl:Class` elements.** Clearing every element wipes
  child text before the parent's `end` event can read it, and `findtext` silently
  returns `''` instead of the label. Cost a debugging round.
- Take **named parents only** from `rdfs:subClassOf` — read `rdf:resource` and
  skip the anonymous `owl:Restriction` children, which carry no resource.

The 28 children of `MP:0000001` "mammalian phenotype" are the body-system
categories (adipose tissue, behavior/neurological, … vision/eye) — the natural
grouping level. **MP is a DAG, not a tree**: 4,477 terms have more than one
parent, a mean of 1.49 categories per term and up to 5. Category counts therefore
overlap and never sum to the total; say so wherever they are displayed.

### Code legends

Verified against https://www.mmrrc.org/methods/data_download.php.

`MUTATION_TYPE`: SM spontaneous · TM targeted · TG transgenic · GT gene trap ·
CI chemically induced · RAD radiation induced · CH chromosomal aberration ·
RB Robertsonian translocation · TL reciprocal translocation · TP transposition ·
INV inversion · INS insertion · DEL deletion · DP duplication · OTH other.

`STRAIN_TYPE`: IS inbred · UN unclassified · MSR mutant strain · MSK mutant
stock · COI coisogenic · SEG segregating inbred · NON noninbred · WDS
wild-derived inbred · CSS consomic/chromosome substitution · RI recombinant
inbred · RC recombinant congenic. (`CON`, 916 strains, appears in the data but
not in the published legend — presumably congenic; confirm before relying on it.)

`STATE` (comma-separated, multi-valued): LM live mouse · CA cryo-archived ·
EM cryopreserved embryos · SP cryopreserved/freeze-dried sperm · ES ES cell
lines · CU currently unavailable.

### Column sparsity (nulls out of 588,980)

| Column | Nulls |
|---|---|
| `MPT_IDS` | 578,098 |
| `MGI_ALLELE_ACCESSION_ID` | 544,820 |
| `PUBMED_IDS` | 490,125 |
| `MGI_GENE_ACCESSION_ID` | 60,220 |
| `GENE_SYMBOL` | 54,923 |
| `OTHER_NAMES` | 0 |

Cross-reference columns of interest: `MGI_GENE_ACCESSION_ID`,
`MGI_ALLELE_ACCESSION_ID`, `PUBMED_IDS`, `MPT_IDS`, `SDS_URL`, and `OTHER_NAMES`
(which carries the `RRID:MMRRC_*` identifier).

---

Below: notes on working with marimo, carried over from earlier pairing sessions.

## Editing the live notebook

**Never edit the notebook's `.py` file directly while a session is running.**
The kernel will overwrite it on save. Always use `marimo._code_mode` (`cm`) from
the `execute-code.sh` scratchpad.

```bash
bash /Users/gaurav/.claude/skills/marimo-pair/scripts/execute-code.sh \
  --url 'http://localhost:2718/' \
  --token "$(cat '/var/folders/.../caa62f-token.txt')" <<'PY'
import marimo._code_mode as cm

async with cm.get_context() as ctx:
    cid = ctx.create_cell("x = 1", hide_code=False)
    ctx.run_cell(cid)
PY
```

Read cell code before editing to avoid staleness errors (`StaleCellError`).
If you're sure the edit is safe to overwrite, pass `skip_staleness_check=True`
to `cm.get_context()`.

## Cell output: always end with an expression

Marimo displays **the last expression** in a cell. `if/else` is a statement, not
an expression, so values created inside branches are silently discarded:

```python
# WRONG — mo.vstack result is discarded; cell shows nothing
if condition:
    mo.vstack([...])   # dropped
else:
    mo.md("empty")     # dropped

# RIGHT — assign in branches, return at end
if condition:
    _output = mo.vstack([...])
else:
    _output = mo.md("empty")

_output   # ← this is the cell's displayed output
```

This applies to any conditional: `if/else`, `match`, etc.

## Private vs. public names

Names prefixed with `_` are private to their cell and invisible to other cells
and to the `execute-code.sh` scratchpad. Use `_` for intermediates that no
downstream cell should read; use public names for anything the reactive graph
needs to track.

## Paths: always relative to the notebook file

Use paths relative to the notebook's own directory:

```python
DATA_CSV = Path("../data/some-file.csv")
```

Never use absolute paths — they break portability.

The marimo **kernel's working directory is the notebook file's own directory**
(where `marimo edit` was launched), so relative paths resolve from there. For a
notebook one level down (e.g. `notebooks/`), a repo-root `data/` dir is
`../data/`. Confirm with `Path(...).exists()` in the kernel if unsure.

## Editing cells programmatically with `marimo._code_mode`

Beyond running code, `cm.get_context()` can durably restructure the notebook:

- `ctx.edit_cell(target, code=..., hide_code=...)` — replace an existing cell.
- `ctx.create_cell(code, before=<id> | after=<id>, hide_code=...)` → new cell id.
- `ctx.delete_cell(target)`, `ctx.move_cell(...)`.

Cell IDs (e.g. `MJUe`) are stable across reorders. Gotchas learned the hard way:

- **An exception anywhere in the `async with` block rolls back every queued
  edit/create in that block.** Don't inspect a just-created cell's
  `ctx.cells[new_id].status` in the same context — the lookup can raise before
  commit and discard your edits. Create/edit in one context, then run and verify
  in a *separate* context.
- `ctx.run_cell()` only runs the named cell (see "Running cells in dependency
  order"); after edits, run upstream→downstream yourself.

## `.removesuffix` vs `.stem` for compound extensions

`Path("data.tar.gz").stem` → `"data.tar"` (wrong — only strips `.gz`).

Use `path.name.removesuffix(".tar.gz")` to strip the full compound suffix.

## Polars and complex types

Polars DataFrames cannot hold arbitrary Python dicts or heterogeneous objects
natively. Serialize them to JSON strings for storage and parse on demand:

```python
"study_mappings": json.dumps(obj.get("metadata", {}).get("study_mappings", {}))
```

For large files, prefer `pl.scan_csv()` (lazy frame) over `pl.read_csv()` so
filters push down into the scan without loading everything into memory.

Use `infer_schema_length=N` (e.g. 6000) when building DataFrames from lists of
dicts where optional fields appear sparsely.

## mo.ui.table selection

`mo.ui.table(df, selection="single")` returns a **polars DataFrame** of the
selected rows via `.value`. When nothing is selected, `.value` is an empty
DataFrame (`len == 0`).

```python
_sel = my_table.value
if not len(_sel):
    _output = mo.md("Nothing selected.")
else:
    _row = _sel["column"][0]
    ...
```

The table is fully reactive: cells that reference `my_table` re-run whenever
the selection changes.

## Running cells in dependency order

`ctx.run_cell()` does **not** automatically resolve stale dependencies — it only
runs the cell you name. If an upstream cell went stale (e.g. because you edited
its dependency), its downstream public names become undefined, causing `NameError`
when the dependent cell runs.

Always run in dependency order manually:

```python
async with cm.get_context() as ctx:
    ctx.run_cell("upstream_cell_id")   # defines `my_var`
    ctx.run_cell("downstream_cell_id") # reads `my_var`
```

To find which cells are stale, check `ctx.cells["id"].status == "stale"` or
inspect the graph with `ctx.graph.ancestors("id")`.

## Debugging a cell that produces no output

If a cell runs without errors but shows nothing in the UI:

1. Check that the last line is a bare expression (not inside `if/else` — see above).
2. Add a debug line that always produces output regardless of branching:
   ```python
   _debug = mo.md(f"sel={some_widget.value!r}")
   # ... rest of cell ...
   _output = mo.vstack([_debug, ...])
   ```
3. Use `ctx.cells["id"].status` and `.errors` to confirm the cell is idle with no errors.

## Polars list aggregation in group_by

Within an eager `group_by().agg()`, you **cannot** chain list operations directly
on an aggregation expression — `.head(3).list.join(", ")` fails with
`InvalidOperationError: attempted list to_struct on non-list dtype`.

The fix: aggregate to a list in `agg()`, then apply list ops in a separate
`with_columns()`:

```python
# WRONG
df.group_by("key").agg(
    pl.col("label").head(3).list.join(", ").alias("examples")  # fails
)

# RIGHT
df.group_by("key")
.agg(pl.col("label").alias("_labels"))          # produces List[String]
.with_columns(
    pl.col("_labels").list.head(3).list.join(", ").alias("examples")
)
.drop("_labels")
```

## Side-by-side widget + detail layout with mo.hstack

A widget defined in one cell doesn't need to be displayed there. Define it without
a trailing display expression, then include it inside `mo.hstack` in the downstream
reactive cell. This achieves a side-by-side layout while respecting marimo's
no-cycles constraint (a cell can't depend on a widget it defines):

```python
# Cell A — defines widget, does NOT display it
variable_usage_data = mo.ui.table(df, selection="single")
# (no bare `variable_usage_data` here)

# Cell B — reactive to widget, displays both side by side
_sel = variable_usage_data.value
_detail = mo.md("select a row") if not len(_sel) else mo.vstack([...])
mo.hstack([variable_usage_data, _detail], widths=[1, 2], align="start")
```

`mo.hstack` accepts `widths` as relative proportions (e.g. `[1, 2]`) and
`align="start"` to top-align panels of unequal height.

## mo.ui.table hidden_columns and visible_columns

`hidden_columns=["col1", "col2"]` hides those columns from the default table view.
They are **not** user-toggleable from the UI — to see them, change to
`visible_columns=[...]` in code. Use this for "examples" columns that clutter the
default view but are useful for deeper inspection.

`visible_columns=["col1", "col2"]` is the inverse: only those columns appear. Use
this to create a slim selector view from a richer DataFrame without defining a
separate derived frame:

```python
# Full table for exploration (all columns):
variable_usage_data = mo.ui.table(full_df, hidden_columns=["examples"])

# Slim selector for the master-detail hstack (3 columns only):
variable_selector = mo.ui.table(
    variable_usage_data._data.drop("_marimo_row_id"),
    visible_columns=["name", "study_usage_count", "concept_count"],
    selection="single",
    page_size=15,
)
```

`widget._data` exposes the underlying DataFrame of any `mo.ui.table`. Marimo
automatically adds a `_marimo_row_id` column — drop it before passing to another
widget to avoid confusing the UI.

## Always define public widgets in all branches

If a downstream cell references a public name (e.g. `concept_table`), that name
must be defined in **every** branch of the defining cell's `if/else`. Otherwise
the downstream cell gets a `NameError` when the "missing" branch executes.

Use an empty typed DataFrame for the "nothing selected" case:

```python
_empty = pl.DataFrame({
    "col_a": pl.Series([], dtype=pl.String),
    "col_b": pl.Series([], dtype=pl.String),
})

if not len(_sel):
    concept_table = mo.ui.table(_empty, selection="single")
    _detail = mo.md("_Select a row._")
else:
    concept_table = mo.ui.table(_real_df, selection="single")
    _detail = mo.vstack([..., concept_table, ...])

_output = ...
_output
```

## Chained reactive drill-downs

You can chain reactive cells to build multi-level drill-downs by exposing each
intermediate widget as a public name:

```
uouT  defines variable_selector  (slim table, not displayed)
SRkg  reads variable_selector.value → shows detail for first selected variable
Cqiu  reads variable_selector.value → shows combined usage for all selected variables
```

Each downstream cell only needs to reference the widget name — marimo handles
the reactive re-run chain automatically.

## Multi-select tables

`mo.ui.table(df, selection="multi")` lets the user select multiple rows.
`.value` still returns a polars DataFrame — just with multiple rows. Check
`len(_sel)` to detect any selection; iterate or call `.to_list()` on a column
to get all selected values.

When multiple rows are selected and you want combined aggregated results (e.g.
union of studies/concepts across all selected variables), read all row IDs with:

```python
_ids = _sel["variable_id"].to_list()
_hits = df.filter(pl.col("variable_id").is_in(_ids))
```

The same cell can serve as both the primary detail view (first-selected item)
and a combined summary for all selected items — just separate the two outputs
into different UI sections.

## Public widgets inside mo.ui.tabs() still power downstream cells

A public widget embedded inside a `mo.ui.tabs()` dict is still a public name
in the cell — marimo tracks it regardless of how deeply it's nested in the
output hierarchy. The downstream cell reacts to `.value` changes just the same:

```python
# Cell A — concept_table is public even though it lives inside tabs
concept_table = mo.ui.table(_concept_hits, selection="multi")
_output = mo.vstack([
    ...,
    mo.ui.tabs({
        "Concepts": concept_table,   # ← widget inside tabs, still public
        "Other": mo.ui.table(_other_df),
    }),
])
_output

# Cell B — reacts whenever the user selects rows in concept_table
_csel = concept_table.value
```

The key is that `concept_table` is assigned at the top level of the cell body
(not as a `_private` name). Its position in the rendered output tree doesn't
affect the dataflow graph.

## Loading small reference files inside a cell

Small reference CSVs (a few KB) can be loaded directly inside a reactive cell
with `pl.read_csv(Path("../path/to/file.csv"))`. There's no need to create a
dedicated upstream cell for data that's only used in one place:

```python
_labels = (
    pl.read_csv(Path("../data/labels.csv"))
    .select(["id", "title"])
    .unique(subset=["id"])
)
```

## Python dict accumulation for complex JSON-column aggregations

When building a table from per-row JSON (e.g. `study_variable_mappings`) where
the aggregation logic is complex or nested, it's cleaner to iterate rows in
Python and accumulate into a dict, then convert to a DataFrame, than to force
everything through Polars expressions.

Pattern: deduplicate by key first, iterate with `iter_rows(named=True)`, parse
JSON, accumulate with `setdefault`, then convert:

```python
_study_map = {}  # study_id -> {selected_var_id -> [local_names]}
for _vid in _var_ids:
    _rows = df.filter(pl.col("variable_id").eq(_vid)).unique(subset=["variable_id"])
    for _row in _rows.iter_rows(named=True):
        _svm = _row["study_variable_mappings"]
        if _svm:
            for _sid, _local_names in json.loads(_svm).items():
                _study_map.setdefault(_sid, {})[_vid] = _local_names

_rows = [
    {"study_id": sid, "mappings": "; ".join(f"{', '.join(lnames)} → {vid}"
                                             for vid, lnames in sorted(vmap.items()))}
    for sid, vmap in sorted(_study_map.items())
]
_df = pl.DataFrame(_rows) if _rows else pl.DataFrame({"study_id": pl.Series([], dtype=pl.String), ...})
```

Use `variables_df.unique(subset=["variable_id"], keep="first")` before iterating
to avoid processing the same variable_id multiple times when it appears in
multiple CRFs.

## Session state

`__marimo__/session/` is ephemeral — it is gitignored. Do not commit it.

## Publishing: `html-wasm` on GitHub Pages

`.github/workflows/pages.yml` exports the notebook with

```bash
marimo export html-wasm gene-mapper-notebook.py -o ../_site --mode run
```

and deploys it to www.ggvaidya.com/mmrrc. The export ships the notebook *source*
plus a Pyodide runtime, so the browser runs every cell and the tables stay fully
interactive — search, sort, paginate, select. Nothing runs at build time and no
data is bundled; the export takes seconds and needs no `data/`.

Two dead ends before this, don't retry them:

- **`marimo export ipynb --include-outputs`**. `mo.ui.table` / `mo.ui.multiselect`
  serialise as `<marimo-table>` / `<marimo-multiselect>` custom elements with
  their data in `data-*` attributes and **no fallback content** — GitHub strips
  them, so the tables vanish rather than degrading to static ones.
- **`marimo export html`**. Tables render but are baked static: no sorting,
  filtering or selection. `marimo export md` is worse still — no
  `--include-outputs` at all, so the narrative appears as source.

### Running the same notebook in the browser

The notebook has to work under both CPython and Pyodide. Three things that
matter, all of them load-bearing:

- **Data comes over the network in WASM.** `data_bytes(name)` returns the local
  `../data/<name>` if it exists and otherwise fetches
  `https://cdn.jsdelivr.net/gh/gaurav/mmrrc@main/data/<name>` — via
  `pyodide.http.pyfetch` under `sys.platform == "emscripten"` (Pyodide's `urllib`
  can't reach the network), else `urllib.request`. `pyfetch` returns the error
  page's body on a 404 instead of raising, so the WASM branch calls
  `raise_for_status()`; `urlopen` already raises.
  jsDelivr over raw.githubusercontent.com for the cache: both send
  `access-control-allow-origin: *`, but raw sends `max-age=300`, so every visit
  re-downloaded all 21 MB, while jsDelivr sends a week. **Two consequences of
  pinning `@main`: a data change has to land on `main` first, and then waits out
  jsDelivr's 12h edge cache.** jsDelivr caps `/gh/` files at 20 MB and
  `mmrrc_catalog_data.csv.gz` is at 16 MB — if it outgrows that, raw still works.
- **`pl.read_csv` goes through pyarrow in WASM.** polars' own CSV reader isn't
  built for emscripten; marimo silently falls back to `pyarrow.csv`, which is why
  the script header carries `pyarrow; sys_platform == 'emscripten'`. That fallback
  reads a buffer, not a path, so it can't sniff gzip — hence the explicit
  `gzip.decompress` before the read.
  Then the trap: **`str.extract_all` on any frame derived from that pyarrow-read
  catalog panics with `PanicException: capacity overflow`** — including a 4,604-row
  filtered slice, and `rechunk()` does not help. Rebuilding the column through
  Python is fine, so the `extract_all()` helper next to `data_bytes` does the
  regex with `re.findall`. Only `extract_all` is affected; `str.extract`,
  `str.contains`, `str.join`, `group_by`, `unique` and friends are all fine.
- **polars in Pyodide is older than the local one.** This surfaced as
  `DataFrame.explode` having no `empty_as_null` there — the notebook no longer
  calls `explode` at all (the phenotype pairs are built in the same Python loop
  that works around `extract_all`), but expect more of these. The
  pyodide-distributed version is whatever Pyodide shipped, and pinning it higher
  in the script header can't change that. Keep to the older API where the two
  spellings mean the same thing.

### Checking it actually runs

The export succeeding proves nothing: every failure above happened at page load,
in the browser, with a green build. Serve the export and drive it headlessly:

```bash
python -m http.server --directory <export-dir>
uvx --from playwright playwright install chromium   # once
```

then a Playwright script that loads the page, waits for a marker string
(`588,980`), and prints console messages. Cell exceptions arrive as console logs
tagged `[STDERR] … (cellId)` with a full Python traceback, which is the only
place the real error appears — the page itself just shows an empty cell. Boot to
catalog-loaded is ~10 s. Cells behind a table selection (the gene detail panel)
never run until something is selected, so click a row before believing them.
