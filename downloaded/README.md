# Downloaded source data

Upstream files this project analyses, committed here so results stay reproducible.
**Both sources evolve** — MMRRC updates its catalog continuously and MP cuts
releases every few weeks — so anything derived from them is a snapshot. Record the
download date whenever you refresh these.

Files are stored gzipped: together they are ~21 MB compressed against ~248 MB raw (~237 MiB),
which is what makes committing them practical. Nothing needs unpacking —
`polars.read_csv` reads `.csv.gz` natively, and the ontology is parsed through
`gzip.open`.

## `mmrrc_catalog_data.csv.gz`

The full MMRRC strain catalog.

| | |
|---|---|
| Source | <https://www.mmrrc.org/about/mmrrc_catalog_data.csv> |
| Linked from | <https://www.mmrrc.org/methods/data_download.php> |
| Downloaded | **2026-07-27 16:53:58 UTC** |
| Size | 16,001,280 bytes gzipped · 147,286,875 bytes raw |
| SHA-256 (of the `.gz`) | `806985c4bcbca4d42d258b60149a0c6106782967d1ef8368d26290cade204605` |
| Contents | 588,980 data rows × 18 columns, 69,388 distinct strains |
| Publisher | Mutant Mouse Resource & Research Centers, <https://www.mmrrc.org/> |

**There is no version stamp in the file**, so the download date above is the only
thing identifying which snapshot this is. Do not drop it.

Columns: `STRAIN/STOCK_ID`, `STRAIN/STOCK_DESIGNATION`, `OTHER_NAMES`,
`STRAIN_TYPE`, `STATE`, `MGI_ALLELE_ACCESSION_ID`, `ALLELE_SYMBOL`, `ALLELE_NAME`,
`MUTATION_TYPE`, `CHROMOSOME`, `MGI_GENE_ACCESSION_ID`, `GENE_SYMBOL`, `GENE_NAME`,
`SDS_URL`, `ACCEPTED_DATE`, `MPT_IDS`, `PUBMED_IDS`, `RESEARCH_AREAS `.

Two quirks worth knowing before touching it, both documented in full in
[`../CLAUDE.md`](../CLAUDE.md):

- `RESEARCH_AREAS ` has a **trailing space** in the header.
- The labels in `MPT_IDS` are corrupt — phenotype names containing a comma are
  torn in two and the names at the end of each affected list are dropped. The MP
  identifiers are intact; read those and take labels from the ontology. Reported
  upstream as [issue #1](https://github.com/gaurav/mmrrc/issues/1).

## `mp.owl.gz`

The Mammalian Phenotype Ontology, used to resolve the `MPT_IDS` identifiers and to
group phenotypes under their parent categories.

| | |
|---|---|
| Source | <https://www.informatics.jax.org/downloads/reports/mp.owl> |
| Downloaded | **2026-07-27 18:16:43 UTC** |
| Ontology release | **2026-07-22** (`owl:versionIRI` → `http://purl.obolibrary.org/obo/mp/releases/2026-07-22/mp.owl`) |
| Size | 5,506,454 bytes gzipped · 101,124,429 bytes raw |
| SHA-256 (of the `.gz`) | `3faaf718c1f747afdf0cfb56c6a65f18207e5484cc63745dd5f74c04656d3a68` |
| Contents | 15,288 MP terms (457 obsolete), 28 top-level categories |
| Publisher | Mouse Genome Informatics, The Jackson Laboratory |
| Licence | CC BY 4.0 — see <https://www.informatics.jax.org/mgihome/other/copyright.shtml> |

Unlike the catalog this one *does* carry a version: prefer the `owl:versionIRI`
release date over the download date when citing it.

The file is RDF/XML and much larger than MP alone, because it merges its imports
(PATO, UBERON, GO, ChEBI, CL) — 126,454 `owl:Class` elements, of which only the
`MP_` ones are ours.

## Refreshing

```bash
curl -sSL https://www.mmrrc.org/about/mmrrc_catalog_data.csv \
  | gzip -9 > downloaded/mmrrc_catalog_data.csv.gz
curl -sSL https://www.informatics.jax.org/downloads/reports/mp.owl \
  | gzip -9 > downloaded/mp.owl.gz
```

Then update the dates, sizes, checksums and counts in this file. The figures above
are quoted throughout the notebook and `CLAUDE.md`, so re-run
`gene-mapper/gene-mapper-notebook.py` and check them after any refresh — several
of the findings are sensitive to how the catalog has changed.
