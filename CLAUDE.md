# CLAUDE.md — MMRRC

This project is a hackathon involving the MMRRC (Mutant Mouse Resource &
Research Centers, https://www.mmrrc.org/). We'll flesh out these notes as we go.

## The MMRRC catalog dataset

`data/mmrrc_catalog_data.csv` (gitignored, ~147 MB) was downloaded from
https://www.mmrrc.org/methods/data_download.php. 588,980 rows × 18 columns.

Header note: `RESEARCH_AREAS ` has a trailing space in the file. Strip it on
load with `pl.read_csv(path).rename(str.strip)`.

### Grain: one row per (strain, allele, gene) — *not* per strain

`STRAIN/STOCK_ID` is **not** unique: 69,388 distinct strains over 588,980 rows
(~8.5 rows/strain, no nulls). Strain-level fields are repeated verbatim on every
row of a strain; only the gene columns vary.

For example `MMRRC:042022-MU` has 206 rows in which `STRAIN/STOCK_DESIGNATION`,
`OTHER_NAMES`, `MGI_ALLELE_ACCESSION_ID`, `ALLELE_SYMBOL`, `MPT_IDS`,
`PUBMED_IDS` and `SDS_URL` are all constant, while `GENE_SYMBOL`, `GENE_NAME`
and `MGI_GENE_ACCESSION_ID` take 206 distinct values across 20 chromosomes — a
lesion spanning many genes.

Alleles per strain: 37,418 strains have 1, 19,699 have 2, 12,150 have 3, and
121 have 4–7. 28,329 strains are a single row.

**Consequence: any per-strain count taken off the raw frame is inflated.** Use
`catalog.unique(subset=["STRAIN/STOCK_ID"])` for strain-level questions, and
treat the raw frame as the strain↔gene edge list.

`(STRAIN/STOCK_ID, ALLELE_SYMBOL, GENE_SYMBOL)` is *nearly* a key — 97 violating
groups covering 229 rows. All of them have both allele and gene null; within
those groups only `MUTATION_TYPE` and `CHROMOSOME` ever differ, and 82 rows are
byte-identical duplicates.

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
