# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "marimo>=0.23.15",
#     "polars>=1.0",
#     # In WASM, marimo routes pl.read_csv through pyarrow -- polars' own CSV
#     # reader isn't available there.
#     "pyarrow; sys_platform == 'emscripten'",
# ]
# ///

import marimo

__generated_with = "0.23.15"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import polars as pl
    import gzip
    import io
    import re
    import sys
    import urllib.request
    import xml.etree.ElementTree as ET
    from collections import defaultdict, deque
    from pathlib import Path
    from urllib.parse import quote_plus

    return (
        ET,
        Path,
        defaultdict,
        deque,
        gzip,
        io,
        mo,
        pl,
        quote_plus,
        re,
        sys,
        urllib,
    )


@app.cell
def _(Path, pl, re, sys, urllib):
    DATA_URL = "https://raw.githubusercontent.com/gaurav/mmrrc/main/data/"


    async def data_bytes(name):
        """The named file from `data/`, local copy if there is one, else GitHub raw.

        The WASM build has neither the repo checkout nor a working urllib, so it
        fetches through the browser instead. raw.githubusercontent.com sends
        `access-control-allow-origin: *`, so the cross-origin fetch is allowed.
        """
        _local = Path("../data") / name
        if _local.exists():
            return _local.read_bytes()
        if sys.platform == "emscripten":
            from pyodide.http import pyfetch

            # pyfetch hands back the error page's body on a 404 rather than
            # raising, which would surface three cells later as "not a gzipped
            # file". urlopen already raises on its own.
            _resp = await pyfetch(DATA_URL + name)
            _resp.raise_for_status()
            return await _resp.bytes()
        return urllib.request.urlopen(DATA_URL + name).read()


    def extract_all(series, pattern):
        """Every match of `pattern` per row, as a List column.

        Not `pl.col(...).str.extract_all`: in WASM the catalog is read through
        pyarrow (polars' own CSV reader isn't built for it), and extract_all on a
        column backed by an arrow buffer panics with "capacity overflow". Both
        columns this is used on have a few thousand non-null rows, so Python's
        `re` is fast enough and behaves identically in both environments.
        """
        return pl.Series(
            series.name,
            [re.findall(pattern, _v) if _v is not None else None for _v in series],
            dtype=pl.List(pl.String),
        )
    return data_bytes, extract_all


@app.cell
async def _(data_bytes, gzip, mo, pl):
    CATALOG_CSV = "mmrrc_catalog_data.csv.gz"

    # Decompress here rather than handing polars the .gz: in WASM the read goes
    # through pyarrow, which won't sniff gzip out of an in-memory buffer.
    catalog = pl.read_csv(gzip.decompress(await data_bytes(CATALOG_CSV))).rename(
        str.strip
    )
    mo.md(f"**{CATALOG_CSV}**: {catalog.height:,} rows x {catalog.width} columns")
    return (catalog,)


@app.cell
def _(catalog, mo):
    catalog_table = mo.ui.table(catalog, page_size=20)
    catalog_table
    return


@app.cell
def _(catalog, mo, pl):
    _ci_strains = catalog.filter(pl.col("MUTATION_TYPE") == "CI")["STRAIN/STOCK_ID"].n_unique()
    _ci_rows = catalog.filter(pl.col("MUTATION_TYPE") == "CI").height
    _ci_genes = (
        catalog.filter(pl.col("MUTATION_TYPE") == "CI")
        .group_by("STRAIN/STOCK_ID")
        .agg(pl.col("GENE_SYMBOL").drop_nulls().n_unique().alias("g"))["g"]
        .mean()
    )

    # MMRRC mutation-type legend, from https://www.mmrrc.org/methods/data_download.php
    MUTATION_LABELS = {
        "Targeted mutation (TM)": "TM",
        "Gene trap (GT)": "GT",
        "Transgenic (TG)": "TG",
        "Deletion (DEL)": "DEL",
        "Spontaneous (SM)": "SM",
        "Insertion (INS)": "INS",
        "Inversion (INV)": "INV",
        "Duplication (DP)": "DP",
        "Transposition (TP)": "TP",
        "Chromosomal aberration (CH)": "CH",
        "Radiation induced (RAD)": "RAD",
        "Chromosomal segment (CS)": "CS",
        "Other (OTH)": "OTH",
        "No mutation type recorded": "(none)",
        "Chemically induced / ENU (CI)": "CI",
    }

    mutation_filter = mo.ui.multiselect(
        options=MUTATION_LABELS,
        value=[_k for _k in MUTATION_LABELS if MUTATION_LABELS[_k] != "CI"],
        label="Mutation types to include",
    )

    mo.vstack([
        mo.md(f"""
    # Gene exploration

    Navigate the catalog gene-first: which genes are manipulated across the most
    strains, how that breaks down by chromosome, and where to follow each gene up
    at MGI and NCBI.

    ## Why chemically-induced strains are excluded by default

    **{_ci_strains:,} chemically-induced (ENU) strains carry a mean of
    {_ci_genes:.1f} genes each**, and between them they produce {_ci_rows:,} of the
    catalog's {catalog.height:,} rows ({_ci_rows / catalog.height:.0%}). Every other
    mutation type averages one to three genes per strain.

    That is because an ENU strain's gene list is a set of *candidate variants found
    by sequencing*, not a set of deliberate manipulations. Longer genes accumulate
    more random ENU hits simply by being bigger targets, so ranking genes by strain
    count with `CI` included returns Ttn, Obscn, Neb, Macf1, Hmcn1 and Syne2 — the
    longest genes in the mouse genome. That ranking measures gene length, not
    research interest.

    So `CI` starts unticked and the tables below describe deliberately manipulated
    genes. **Tick it back on to see the ENU picture** — the ranking changes
    completely, which is itself the point.
    """),
        mutation_filter,
    ])
    return (mutation_filter,)


@app.cell
def _(catalog, mo, pl):
    MOUSE_CHROMS = [str(_i) for _i in range(1, 20)] + ["X", "Y", "MT"]

    # CHROMOSOME is free text: 50 distinct values for what should be 22, including
    # "Chr 1", "Chr11:4938754-4948064 bp", "8q21.13", "unk" and "N/A". Recover what
    # is recoverable and bucket the rest as "unmapped".
    _chrom = (
        pl.col("CHROMOSOME")
        .str.strip_chars()
        .str.to_uppercase()
        .str.replace(r"^CHR\s*", "")
        .str.strip_chars()
    )

    # Rows come in two disjoint kinds: gene rows (GENE_SYMBOL and
    # MGI_GENE_ACCESSION_ID set, allele columns null) and allele rows (the reverse).
    # Zero rows carry both, so a gene is linked to its alleles only via the strain.
    catalog_genes = catalog.filter(pl.col("GENE_SYMBOL").is_not_null()).with_columns(
        pl.when(_chrom.is_in(MOUSE_CHROMS))
        .then(_chrom)
        .otherwise(pl.lit("unmapped"))
        .alias("chrom"),
        pl.col("MUTATION_TYPE").fill_null("(none)").alias("mut"),
    )

    strain_alleles = (
        catalog.filter(pl.col("ALLELE_SYMBOL").is_not_null())
        .select(["STRAIN/STOCK_ID", "ALLELE_SYMBOL", "MGI_ALLELE_ACCESSION_ID"])
        .unique()
    )

    mo.md(
        f"`catalog_genes`: {catalog_genes.height:,} gene rows "
        f"({catalog_genes['GENE_SYMBOL'].n_unique():,} distinct symbols) &nbsp;·&nbsp; "
        f"`strain_alleles`: {strain_alleles.height:,} allele rows "
        f"({strain_alleles['ALLELE_SYMBOL'].n_unique():,} distinct alleles)"
    )
    return MOUSE_CHROMS, catalog_genes, strain_alleles


@app.cell
def _(catalog_genes, mo, mutation_filter, pl, quote_plus, strain_alleles):
    def mgi_url(mgi_id):
        """Link to an MGI marker detail page, e.g. MGI:98864 -> .../marker/MGI:98864."""
        return f"https://www.informatics.jax.org/marker/{mgi_id}"


    def ncbi_url(symbol, mgi_id=None):
        """Link to an NCBI Gene record via search.

        The catalog carries no Entrez id, but `Ttn[sym] AND "Mus musculus"[orgn]`
        resolves straight to the record page rather than a result list. Symbols with
        no MGI id are mostly non-mouse transgenes (SOD1, APP), which would not match
        the mouse organism clause, so it is dropped for them.
        """
        _term = f"{symbol}[sym]"
        if mgi_id is not None:
            _term += ' AND "Mus musculus"[orgn]'
        return "https://www.ncbi.nlm.nih.gov/gene/?term=" + quote_plus(_term)


    _selected_muts = mutation_filter.value
    _kept = catalog_genes.filter(pl.col("mut").is_in(_selected_muts))

    # Alleles reach a gene only through the strains that carry it (see catalog_genes).
    _gene_alleles = (
        _kept.select(["GENE_SYMBOL", "STRAIN/STOCK_ID"])
        .unique()
        .join(strain_alleles, on="STRAIN/STOCK_ID", how="left")
        .group_by("GENE_SYMBOL")
        .agg(pl.col("ALLELE_SYMBOL").drop_nulls().n_unique().alias("alleles"))
    )

    # Kept unfiltered so the gap against `strains` shows how much of a gene's
    # presence is ENU noise.
    _all_strains = catalog_genes.group_by("GENE_SYMBOL").agg(
        pl.col("STRAIN/STOCK_ID").n_unique().alias("strains_all")
    )

    gene_index = (
        _kept.group_by("GENE_SYMBOL")
        .agg(
            pl.col("STRAIN/STOCK_ID").n_unique().alias("strains"),
            pl.col("chrom").first().alias("chrom"),
            pl.col("GENE_NAME").drop_nulls().first().alias("gene_name"),
            pl.col("MGI_GENE_ACCESSION_ID").drop_nulls().first().alias("mgi_id"),
            pl.col("mut").unique().sort().str.join(", ").alias("mutation_types"),
        )
        .join(_gene_alleles, on="GENE_SYMBOL", how="left")
        .join(_all_strains, on="GENE_SYMBOL", how="left")
        .rename({"GENE_SYMBOL": "gene_symbol"})
        .with_columns(
            pl.struct(["gene_symbol", "mgi_id"])
            .map_elements(
                lambda r: ncbi_url(r["gene_symbol"], r["mgi_id"]), return_dtype=pl.String
            )
            .alias("ncbi")
        )
        .select([
            "gene_symbol", "gene_name", "chrom", "strains", "strains_all",
            "alleles", "mutation_types", "mgi_id", "ncbi",
        ])
        .sort("strains", descending=True)
    )

    mo.md(
        f"`gene_index`: **{gene_index.height:,} genes** under the current mutation-type "
        f"filter, out of {catalog_genes['GENE_SYMBOL'].n_unique():,} in the catalog."
    )
    return gene_index, mgi_url, ncbi_url


@app.cell
def _(MOUSE_CHROMS, catalog_genes, mo, mutation_filter, pl):
    _ck = catalog_genes.filter(pl.col("mut").is_in(mutation_filter.value))

    # Sort naturally (1..19, X, Y, MT, unmapped) rather than by count, so the table
    # reads like a karyotype; it is still click-sortable by any column.
    _order = {_c: _i for _i, _c in enumerate(MOUSE_CHROMS + ["unmapped"])}

    chrom_summary = (
        _ck.group_by("chrom")
        .agg(
            pl.col("GENE_SYMBOL").n_unique().alias("genes"),
            pl.col("STRAIN/STOCK_ID").n_unique().alias("strains"),
        )
        .join(
            _ck.group_by(["chrom", "GENE_SYMBOL"])
            .agg(pl.col("STRAIN/STOCK_ID").n_unique().alias("_s"))
            .sort("_s", descending=True)
            .group_by("chrom")
            .agg(pl.col("GENE_SYMBOL").first().alias("top_gene")),
            on="chrom",
            how="left",
        )
        .with_columns(
            pl.col("chrom").replace_strict(_order, default=99).alias("_o")
        )
        .sort("_o")
        .drop("_o")
    )

    chrom_table = mo.ui.table(
        chrom_summary,
        selection="multi",
        page_size=25,
        label="**Chromosome** — select to filter the gene table",
    )
    return (chrom_table,)


@app.cell
def _(chrom_table, gene_index, mgi_url, mo, pl):
    _csel = chrom_table.value
    _chroms = _csel["chrom"].to_list() if len(_csel) else None
    _shown = gene_index if _chroms is None else gene_index.filter(pl.col("chrom").is_in(_chroms))

    gene_table = mo.ui.table(
        _shown,
        selection="single",
        page_size=15,
        label=(
            f"**Genes** — {_shown.height:,} shown"
            + ("" if _chroms is None else f", chromosome {', '.join(_chroms)}")
        ),
        format_mapping={
            "mgi_id": lambda v: mo.md(f"[{v}]({mgi_url(v)})" if v else "—"),
            "ncbi": lambda v: mo.md(f"[NCBI Gene]({v})"),
        },
    )

    mo.hstack([chrom_table, gene_table], widths=[1, 2], align="start")
    return (gene_table,)


@app.cell
def _(
    catalog_genes,
    extract_all,
    gene_table,
    mgi_url,
    mo,
    mutation_filter,
    ncbi_url,
    pl,
    strain_alleles,
):
    _gsel = gene_table.value

    if not len(_gsel):
        _output = mo.md("_Select a gene above to see its strains and cross-references._")
    else:
        _g = _gsel["gene_symbol"][0]
        _mgi = _gsel["mgi_id"][0]
        _rows = catalog_genes.filter(
            pl.col("GENE_SYMBOL").eq(_g) & pl.col("mut").is_in(mutation_filter.value)
        )

        _links = [f"[NCBI Gene]({ncbi_url(_g, _mgi)})"]
        if _mgi:
            _links.insert(0, f"[{_mgi}]({mgi_url(_mgi)})")

        _by_mut = (
            _rows.group_by("mut")
            .agg(pl.col("STRAIN/STOCK_ID").n_unique().alias("strains"))
            .sort("strains", descending=True)
        )

        _base = _rows.unique(subset=["STRAIN/STOCK_ID"]).join(
            strain_alleles.group_by("STRAIN/STOCK_ID").agg(
                pl.col("ALLELE_SYMBOL").unique().sort().str.join(", ").alias("alleles")
            ),
            on="STRAIN/STOCK_ID",
            how="left",
        )

        _strains = (
            _base.with_columns(
                pl.col("OTHER_NAMES").str.extract(r"(RRID:MMRRC_[\w.-]+)").alias("rrid"),
                extract_all(_base["PUBMED_IDS"], r"\d{6,9}").alias("_pmids"),
            )
            .select([
                "STRAIN/STOCK_ID", "STRAIN/STOCK_DESIGNATION", "alleles",
                "mut", "STRAIN_TYPE", "STATE", "rrid", "_pmids", "SDS_URL",
            ])
            .rename({
                "STRAIN/STOCK_ID": "strain_id",
                "STRAIN/STOCK_DESIGNATION": "designation",
                "_pmids": "pubmed",
            })
            .sort("strain_id")
        )

        _output = mo.vstack([
            mo.md(
                f"### {_g} &nbsp; <small>{_gsel['gene_name'][0] or ''}</small>\n\n"
                f"Chromosome **{_gsel['chrom'][0]}** &nbsp;·&nbsp; "
                f"**{_gsel['strains'][0]:,}** strains under the current filter "
                f"({_gsel['strains_all'][0]:,} across all mutation types) &nbsp;·&nbsp; "
                f"**{_gsel['alleles'][0] or 0:,}** alleles\n\n"
                + " &nbsp;·&nbsp; ".join(_links)
            ),
            mo.md("**Strains by mutation type**"),
            mo.ui.table(_by_mut, selection=None, page_size=8),
            mo.md("**Strains carrying this gene**"),
            mo.ui.table(
                _strains,
                selection=None,
                page_size=10,
                format_mapping={
                    "rrid": lambda v: mo.md(f"[{v}](https://scicrunch.org/resolver/{v})" if v else "—"),
                    "SDS_URL": lambda v: mo.md(f"[data sheet]({v})" if v else "—"),
                    "pubmed": lambda v: mo.md(
                        ", ".join(f"[{p}](https://pubmed.ncbi.nlm.nih.gov/{p}/)" for p in v)
                        if v is not None and len(v)
                        else "—"
                    ),
                },
            ),
        ])

    _output
    return


@app.cell
def _(catalog, catalog_genes, mo, pl):
    _g = catalog_genes  # all gene rows, filter-independent
    _ci = _g.filter(pl.col("mut") == "CI")
    _non_ci = _g.filter(pl.col("mut") != "CI")
    _n_genes = _g["GENE_SYMBOL"].n_unique()
    _ci_only = _n_genes - _non_ci["GENE_SYMBOL"].n_unique()

    _caps = pl.col("GENE_SYMBOL").str.contains(r"[a-z]").not_() & pl.col(
        "GENE_SYMBOL"
    ).str.contains(r"[A-Z]")
    _sym = _g.select(["GENE_SYMBOL", "MGI_GENE_ACCESSION_ID"]).unique(subset=["GENE_SYMBOL"])

    _top_ci = _ci.group_by("GENE_SYMBOL").agg(
        pl.col("STRAIN/STOCK_ID").n_unique().alias("s")
    ).sort("s", descending=True)["GENE_SYMBOL"].head(5).to_list()
    _top_non_ci = _non_ci.group_by("GENE_SYMBOL").agg(
        pl.col("STRAIN/STOCK_ID").n_unique().alias("s")
    ).sort("s", descending=True)["GENE_SYMBOL"].head(5).to_list()

    _by_chrom = _g.group_by("chrom").agg(
        pl.col("GENE_SYMBOL").n_unique().alias("genes"),
        pl.col("STRAIN/STOCK_ID").n_unique().alias("strains"),
    )
    _placed = _by_chrom.filter(pl.col("chrom") != "unmapped")
    _most_genes = _placed.sort("genes", descending=True).row(0, named=True)
    _most_strains = _placed.sort("strains", descending=True).row(0, named=True)

    mo.md(f"""
    ## Things worth knowing about this data

    **1. The top of the raw ranking is a gene-length artifact.** Ranking by strain
    count with ENU included gives {", ".join(_top_ci)} — among the longest genes in the
    mouse genome, hit most often by random mutagenesis simply for being the biggest
    targets. Worse, **{_ci_only:,} of {_n_genes:,} genes
    ({_ci_only / _n_genes:.0%}) appear *only* in chemically-induced strains** — they
    have never been deliberately manipulated in this catalog at all.

    **2. The most-manipulated "genes" are reagents, not disease genes.** Excluding ENU,
    the ranking is {", ".join(_top_non_ci)} — the toolkit of mouse genetics.
    `Hprt1` is the ES-cell HAT-selection locus (and sits on the X); `EGFP`, `cre`,
    `lacZ` and `tTA` are cassettes with **no MGI id and no chromosome**, so they land
    in the `unmapped` bucket. `Gt(ROSA)26Sor` is the canonical safe-harbour locus.
    Genuine biology starts several rows down.

    **3. Human transgenes are hiding in the gene column.**
    {_sym.filter(_caps).height:,} symbols are ALL-CAPS
    ({_sym.filter(_caps & pl.col("MGI_GENE_ACCESSION_ID").is_null()).height:,} of them
    with no MGI id) — non-mouse genes carried as transgenes. The giveaway is that some
    sit on "chromosome" 20, 21 and 22, which mice do not have: `SOD1` and `APP` on
    chr21 are *human* coordinates, from ALS and Alzheimer models. The NCBI links above
    drop the mouse organism filter for these so they still resolve.

    **4. `CHROMOSOME` is free text, not a category.**
    {catalog["CHROMOSOME"].n_unique():,} distinct values for what should be 22,
    including `unknown`, `UN`, `unk`, `N/A`, `Chr 1`, `Chr11:4938754-4948064 bp`,
    `8q21.13` (a *human* cytoband), `919`, and one entry where chromosome 14 was
    typed with a stray backtick. Normalisation recovers most of it and leaves
    {_by_chrom.filter(pl.col("chrom") == "unmapped")["strains"].sum():,} strains
    unmapped — but that bucket is not only dirt, it is also where the reagent
    cassettes legitimately live.

    **5. Coverage is uneven across the genome.** Chromosome {_most_genes["chrom"]}
    carries the most distinct genes ({_most_genes["genes"]:,}) while chromosome
    {_most_strains["chrom"]} carries the most strains ({_most_strains["strains"]:,}) —
    gene density and research attention are not the same thing.
    """)
    return


@app.cell
def _(mo):
    mo.md("""
    # Phenotype exploration

    `MPT_IDS` holds Mammalian Phenotype annotations as pipe-separated
    `label [MP:id]` pairs:

    > `decreased bone mineral density [MP:0000063]| abnormal vertebrae morphology [MP:0000137]| …`

    This section cross-links them against `data/mp.owl`, the Mammalian Phenotype
    Ontology, to group phenotypes under their parent categories and see which areas
    of mouse biology the collection actually covers.

    **Only the ids in that column are usable.** MMRRC's export re-splits the joined
    label string on `", "`, so any phenotype name containing a comma is torn in two
    and the tail of each affected list is lost. The ids are untouched, so the cells
    below read ids only and take every label from the ontology. The highlights at the
    end of the section show the evidence
    ([issue #1](https://github.com/gaurav/mmrrc/issues/1)).
    """)
    return


@app.cell
async def _(ET, data_bytes, defaultdict, deque, gzip, io, mo, pl):
    MP_OWL = "mp.owl.gz"
    MP_ROOT = "MP:0000001"

    _OBO = "http://purl.obolibrary.org/obo/"
    _RDF = "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}"
    _RDFS = "{http://www.w3.org/2000/01/rdf-schema#}"
    _OWL = "{http://www.w3.org/2002/07/owl#}"

    mp_labels = {}
    mp_obsolete = set()
    mp_parents = defaultdict(list)

    # One streaming pass over 101 MB of RDF/XML (5 MB gzipped), ~1.5s -- no ontology
    # library needed. Only clear owl:Class elements: clearing every element wipes
    # child text before the parent's end event can read it.
    with gzip.open(io.BytesIO(await data_bytes(MP_OWL))) as _fh:
        for _ev, _el in ET.iterparse(_fh, events=("end",)):
            if _el.tag != _OWL + "Class":
                continue
            _about = _el.get(_RDF + "about", "")
            if _about.startswith(_OBO + "MP_"):
                _cid = "MP:" + _about.rsplit("MP_", 1)[1]
                _label = _el.findtext(_RDFS + "label")
                if _label:
                    mp_labels[_cid] = _label
                if _el.findtext(_OWL + "deprecated") == "true":
                    mp_obsolete.add(_cid)
                for _sub in _el.findall(_RDFS + "subClassOf"):
                    # Named parents only; anonymous owl:Restriction children carry no
                    # rdf:resource and are skipped.
                    _r = _sub.get(_RDF + "resource")
                    if _r and _r.startswith(_OBO + "MP_"):
                        mp_parents[_cid].append("MP:" + _r.rsplit("MP_", 1)[1])
            _el.clear()

    # The 28 children of "mammalian phenotype" -- the body-system grouping.
    mp_categories = {
        _c: mp_labels[_c]
        for _c in sorted(
            (_k for _k, _ps in mp_parents.items() if MP_ROOT in _ps),
            key=lambda _k: mp_labels[_k],
        )
    }


    def mp_ancestors(mp_id):
        """Every ancestor of a term. MP is a DAG, so a term can have several parents."""
        _seen, _out, _q = {mp_id}, set(), deque([mp_id])
        while _q:
            for _p in mp_parents.get(_q.popleft(), ()):
                if _p not in _seen:
                    _seen.add(_p)
                    _out.add(_p)
                    _q.append(_p)
        return _out


    _cat_ids = set(mp_categories)
    mp_rollup = pl.DataFrame(
        [
            {"mp_id": _c, "category_id": _k}
            for _c in mp_labels
            for _k in ({_c} | mp_ancestors(_c)) & _cat_ids
        ],
        schema={"mp_id": pl.String, "category_id": pl.String},
    )

    mo.md(
        f"`mp.owl`: **{len(mp_labels):,} MP terms** ({len(mp_obsolete)} obsolete), "
        f"**{len(mp_categories)} top-level categories**, "
        f"{mp_rollup.height:,} term→category edges "
        f"({mp_rollup.height / mp_rollup['mp_id'].n_unique():.2f} categories per term — "
        f"MP is a DAG, so these overlap)."
    )
    return mp_categories, mp_labels, mp_obsolete, mp_rollup


@app.cell
def _(catalog, mo, mp_labels, pl, re):
    # Take the MP ids and nothing else. The label half of MPT_IDS is unusable:
    # MMRRC joins the phenotype names into one comma-separated string, re-splits it on
    # ", " and truncates to the number of ids -- so any label containing a comma is
    # torn in two and the tail of the list is dropped. The ids are untouched: complete,
    # correctly ordered, and every one a valid MP term. Labels come from the ontology.
    # See https://github.com/gaurav/mmrrc/issues/1.
    _annotated = catalog.filter(pl.col("MPT_IDS").is_not_null()).unique(
        subset=["STRAIN/STOCK_ID"]
    )

    # Built in Python rather than extract_all + explode: polars' extract_all panics
    # on this frame under WASM (see the WASM notes in CLAUDE.md), and 46k pairs is
    # nothing to iterate.
    strain_phenotypes = pl.DataFrame(
        [
            (_sid, _mp)
            for _sid, _ids in zip(
                _annotated["STRAIN/STOCK_ID"], _annotated["MPT_IDS"]
            )
            for _mp in re.findall(r"MP:\d+", _ids)
        ],
        schema=["strain_id", "mp_id"],
        orient="row",
    ).unique()

    # Diagnostics for the highlights below: how badly the label half is mangled, and
    # confirmation that the ids are the intact half.
    _dmg = dict(strains=0, ids=0, truncated_strains=0, lost_slots=0, pos=0, pos_match=0)
    for _row in (
        catalog.filter(pl.col("MPT_IDS").is_not_null())
        .unique(subset=["STRAIN/STOCK_ID"])
        .select(["MPT_IDS"])
        .iter_rows(named=True)
    ):
        _entries = [_e.strip() for _e in _row["MPT_IDS"].split("|")]
        _ids = re.findall(r"MP:\d+", _row["MPT_IDS"])
        _texts = [re.sub(r"\s*\[MP:\d+\]\s*$", "", _e) for _e in _entries]
        _dmg["strains"] += 1
        _dmg["ids"] += len(_ids)
        if any(_i not in mp_labels for _i in _ids):
            continue
        # Reconstruct what MMRRC's exporter did: join the real labels, re-split on
        # ", ", keep only as many pieces as there are ids.
        _expanded = ", ".join(mp_labels[_i] for _i in _ids).split(", ")
        if len(_expanded) > len(_ids):
            _dmg["truncated_strains"] += 1
            _dmg["lost_slots"] += len(_expanded) - len(_ids)
        for _a, _b in zip(_texts, _expanded):
            _dmg["pos"] += 1
            _dmg["pos_match"] += _a.lower() == _b.lower()

    mp_label_damage = _dmg

    mo.md(
        f"`strain_phenotypes`: **{strain_phenotypes.height:,} strain→phenotype edges** "
        f"over {strain_phenotypes['strain_id'].n_unique():,} strains and "
        f"{strain_phenotypes['mp_id'].n_unique():,} distinct MP terms, taken from the ids "
        f"alone. Reconstructing MMRRC's mangled label column from those ids reproduces "
        f"**{_dmg['pos_match'] / _dmg['pos']:.1%}** of its {_dmg['pos']:,} label slots — "
        f"the ids are the intact half."
    )
    return mp_label_damage, strain_phenotypes


@app.cell
def _(
    catalog,
    catalog_genes,
    mo,
    mp_categories,
    mp_labels,
    mp_obsolete,
    mp_rollup,
    pl,
    strain_phenotypes,
):
    def mp_url(mp_id):
        """Link to the MGI Mammalian Phenotype browser for a term."""
        return f"https://www.informatics.jax.org/vocab/mp_ontology/{mp_id}"


    _cat_names = pl.DataFrame(
        {"category_id": list(mp_categories), "category": list(mp_categories.values())}
    )

    _term_cats = (
        mp_rollup.join(_cat_names, on="category_id", how="inner")
        .group_by("mp_id")
        .agg(pl.col("category").unique().sort().str.join(", ").alias("categories"))
    )

    # Genes reach a phenotype through the strain, exactly as alleles do (see the
    # grain notes) -- catalog_genes unfiltered, so this does not inherit the gene
    # section's mutation-type filter.
    _strain_genes = (
        catalog_genes.select(["STRAIN/STOCK_ID", "GENE_SYMBOL"])
        .unique()
        .rename({"STRAIN/STOCK_ID": "strain_id"})
    )

    phenotype_index = (
        strain_phenotypes.group_by("mp_id")
        .agg(pl.col("strain_id").n_unique().alias("strains"))
        .join(
            strain_phenotypes.join(_strain_genes, on="strain_id", how="inner")
            .group_by("mp_id")
            .agg(pl.col("GENE_SYMBOL").n_unique().alias("genes")),
            on="mp_id",
            how="left",
        )
        .join(_term_cats, on="mp_id", how="left")
        .with_columns(
            pl.col("mp_id").replace_strict(mp_labels, default=None).alias("label"),
            pl.col("mp_id").is_in(list(mp_obsolete)).alias("obsolete"),
            pl.col("genes").fill_null(0),
        )
        .select(["label", "mp_id", "strains", "genes", "categories", "obsolete"])
        .sort("strains", descending=True)
    )

    mo.md(
        f"`phenotype_index`: **{phenotype_index.height:,} distinct phenotypes** used by "
        f"{strain_phenotypes['strain_id'].n_unique():,} of "
        f"{catalog['STRAIN/STOCK_ID'].n_unique():,} strains "
        f"({strain_phenotypes['strain_id'].n_unique() / catalog['STRAIN/STOCK_ID'].n_unique():.1%})."
    )
    return mp_url, phenotype_index


@app.cell
def _(mo, mp_categories, mp_rollup, pl, strain_phenotypes):
    _cat_names2 = pl.DataFrame(
        {"category_id": list(mp_categories), "category": list(mp_categories.values())}
    )

    category_summary = (
        strain_phenotypes.join(mp_rollup, on="mp_id", how="inner")
        .group_by("category_id")
        .agg(
            pl.col("mp_id").n_unique().alias("phenotypes"),
            pl.col("strain_id").n_unique().alias("strains"),
        )
        .join(_cat_names2, on="category_id", how="inner")
        .select(["category", "strains", "phenotypes", "category_id"])
        .sort("strains", descending=True)
    )

    category_table = mo.ui.table(
        category_summary,
        selection="multi",
        page_size=30,
        label=(
            "**MP category** — select to filter. A term can sit under several "
            "categories, so these columns overlap and do not sum to the total."
        ),
    )
    return (category_table,)


@app.cell
def _(category_table, mo, mp_rollup, mp_url, phenotype_index, pl):
    _csel = category_table.value
    _cats = _csel["category_id"].to_list() if len(_csel) else None

    if _cats is None:
        _shown = phenotype_index
    else:
        _ids = mp_rollup.filter(pl.col("category_id").is_in(_cats))["mp_id"].unique().to_list()
        _shown = phenotype_index.filter(pl.col("mp_id").is_in(_ids))

    phenotype_table = mo.ui.table(
        _shown,
        selection="single",
        page_size=15,
        label=(
            f"**Phenotypes** — {_shown.height:,} shown"
            + ("" if _cats is None else ", filtered by category")
        ),
        format_mapping={
            "mp_id": lambda v: mo.md(f"[{v}]({mp_url(v)})" if v else "—"),
        },
    )

    mo.hstack([category_table, phenotype_table], widths=[1, 2], align="start")
    return (phenotype_table,)


@app.cell
def _(
    catalog,
    catalog_genes,
    mgi_url,
    mo,
    mp_url,
    ncbi_url,
    phenotype_table,
    pl,
    strain_phenotypes,
):
    _psel = phenotype_table.value

    if not len(_psel):
        _output = mo.md("_Select a phenotype above to see its genes and strains._")
    else:
        _mp = _psel["mp_id"][0]
        _strain_ids = strain_phenotypes.filter(pl.col("mp_id") == _mp)["strain_id"].to_list()

        _genes = (
            catalog_genes.filter(pl.col("STRAIN/STOCK_ID").is_in(_strain_ids))
            .group_by("GENE_SYMBOL")
            .agg(
                pl.col("STRAIN/STOCK_ID").n_unique().alias("strains"),
                pl.col("MGI_GENE_ACCESSION_ID").drop_nulls().first().alias("mgi_id"),
                pl.col("chrom").first().alias("chrom"),
            )
            .rename({"GENE_SYMBOL": "gene_symbol"})
            .sort("strains", descending=True)
        )

        _strains = (
            catalog.filter(pl.col("STRAIN/STOCK_ID").is_in(_strain_ids))
            .unique(subset=["STRAIN/STOCK_ID"])
            .with_columns(
                pl.col("OTHER_NAMES").str.extract(r"(RRID:MMRRC_[\w.-]+)").alias("rrid")
            )
            .select([
                "STRAIN/STOCK_ID", "STRAIN/STOCK_DESIGNATION",
                "MUTATION_TYPE", "STRAIN_TYPE", "STATE", "rrid", "SDS_URL",
            ])
            .rename({
                "STRAIN/STOCK_ID": "strain_id",
                "STRAIN/STOCK_DESIGNATION": "designation",
                "MUTATION_TYPE": "mut",
            })
            .sort("strain_id")
        )

        _output = mo.vstack([
            mo.md(
                f"### {_psel['label'][0]}\n\n"
                f"[{_mp}]({mp_url(_mp)}) &nbsp;·&nbsp; "
                f"**{_psel['strains'][0]:,}** strains &nbsp;·&nbsp; "
                f"**{_psel['genes'][0]:,}** genes &nbsp;·&nbsp; "
                f"_{_psel['categories'][0] or 'no category'}_"
            ),
            mo.md("**Genes most associated with this phenotype**"),
            mo.ui.table(
                _genes,
                selection=None,
                page_size=8,
                format_mapping={
                    "mgi_id": lambda v: mo.md(f"[{v}]({mgi_url(v)})" if v else "—"),
                    "gene_symbol": lambda v: mo.md(f"[{v}]({ncbi_url(v)})"),
                },
            ),
            mo.md("**Strains showing it**"),
            mo.ui.table(
                _strains,
                selection=None,
                page_size=10,
                format_mapping={
                    "rrid": lambda v: mo.md(
                        f"[{v}](https://scicrunch.org/resolver/{v})" if v else "—"
                    ),
                    "SDS_URL": lambda v: mo.md(f"[data sheet]({v})" if v else "—"),
                },
            ),
        ])

    _output
    return


@app.cell
def _(
    catalog,
    mo,
    mp_label_damage,
    mp_labels,
    mp_rollup,
    phenotype_index,
    pl,
    strain_phenotypes,
):
    _d = mp_label_damage
    _all_strains = catalog["STRAIN/STOCK_ID"].n_unique()
    _ann_ids = strain_phenotypes["strain_id"].unique().to_list()
    _ann = len(_ann_ids)

    _mut = (
        catalog.filter(pl.col("STRAIN/STOCK_ID").is_in(_ann_ids))
        .group_by("STRAIN/STOCK_ID")
        .agg(pl.col("MUTATION_TYPE").drop_nulls().unique().sort().str.join("+").alias("mut"))
        .group_by("mut")
        .agg(pl.len().alias("n"))
    )
    _tm = _mut.filter(pl.col("mut") == "TM")["n"].sum()
    _ci = _mut.filter(pl.col("mut") == "CI")["n"].sum()
    _top2 = phenotype_index.head(2)

    mo.md(f"""
    ## Things worth knowing about the phenotype data

    **1. The labels in `MPT_IDS` are corrupt; the ids are fine. Use the ids.**
    MMRRC's exporter joins the phenotype names into one comma-separated string,
    re-splits that string on `", "`, and zips the pieces against the id list,
    truncating to its length. Because MP labels legitimately contain commas
    (`decreased CD4-positive, alpha-beta T cell number`), every such label is torn in
    two, everything after it slides by one, and the tail of the label list falls off
    the end.

    Reconstructing the column from the ids alone — join the ontology's labels, split
    on `", "`, truncate — reproduces **{_d["pos_match"]:,} of {_d["pos"]:,}
    ({_d["pos_match"] / _d["pos"]:.1%})** of the label slots in the file, which is
    what identifies the mechanism. The residual is ordinary label drift (point 2).

    {_d["truncated_strains"]:,} of {_d["strains"]:,} annotated strains are affected,
    losing {_d["lost_slots"]:,} label slots. All {_d["ids"]:,} ids are valid MP terms.

    `MMRRC:011644-UNC` is the shape of it — four ids, but only enough label text for
    the first four fragments of three of them:

    | Entry in the file | The id is right | The label beside it is not |
    |---|---|---|
    | `abnormal trophoblast giant cell morphology [MP:0005033]` | abnormal trophoblast giant cell morphology | ✅ |
    | `embryonic lethality between implantation and somite formation [MP:0011096]` | …, **complete penetrance** | truncated at the comma |
    | `complete penetrance [MP:0011100]` | preweaning lethality, complete penetrance | the other half of the line above |
    | `preweaning lethality [MP:0012113]` | **decreased inner cell mass proliferation** | label ran out; text is a leftover |

    The strain really does have `MP:0012113 decreased inner cell mass proliferation` —
    its name never appears in the file because the label expansion was cut short.

    **2. The labels are also a stale snapshot.** Independent of the mangling, MMRRC's
    text lags the ontology — `hypoactivity` is now `decreased locomotor activity`,
    `retinal degeneration` is now `retina degeneration`, `thyroid inflammation` is now
    `thyroid gland inflammation`, `aorta dilation` is now `dilated aorta`. Another
    reason to take ids and render labels from `mp.owl`, as every table above does.

    **3. Phenotype coverage is thin and skewed.** Only **{_ann:,} of
    {_all_strains:,} strains ({_ann / _all_strains:.1%})** carry any phenotype, and
    they lean toward deliberately characterised lines — {_tm:,} purely
    targeted-mutation strains against {_ci:,} chemically-induced, in a catalog whose
    largest single group is gene traps. Only
    {strain_phenotypes["mp_id"].n_unique():,} of {len(mp_labels):,} MP terms
    ({strain_phenotypes["mp_id"].n_unique() / len(mp_labels):.0%}) are used at all.
    Absence of a phenotype here means absence of *characterisation*, never absence of
    an effect.

    **4. The second most common "phenotype" is the absence of one.**
    `{_top2["label"][1]}` ({_top2["mp_id"][1]}) sits on {_top2["strains"][1]:,}
    strains, just behind `{_top2["label"][0]}` at {_top2["strains"][0]:,}. It is a
    negative result, not a phenotype. It is left in the table — the ontology files it
    under *normal phenotype*, so the category column flags it — but any ranking that
    treats it as a finding is wrong.

    **5. Categories overlap by design.** MP is a DAG, not a tree:
    {mp_rollup.height:,} term→category edges across
    {mp_rollup["mp_id"].n_unique():,} terms, a mean of
    {mp_rollup.height / mp_rollup["mp_id"].n_unique():.2f} categories per term and up
    to 5. One strain with one phenotype can count toward several categories, so the
    category table's columns never sum to the totals.
    """)
    return


if __name__ == "__main__":
    app.run()
