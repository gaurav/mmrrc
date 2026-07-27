import marimo

__generated_with = "0.23.15"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import polars as pl
    from pathlib import Path
    from urllib.parse import quote_plus

    return Path, mo, pl, quote_plus


@app.cell
def _(Path, mo, pl):
    CATALOG_CSV = Path("../data/mmrrc_catalog_data.csv")

    catalog = pl.read_csv(CATALOG_CSV).rename(str.strip)
    mo.md(f"**{CATALOG_CSV.name}**: {catalog.height:,} rows x {catalog.width} columns")
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

        _strains = (
            _rows.unique(subset=["STRAIN/STOCK_ID"])
            .join(
                strain_alleles.group_by("STRAIN/STOCK_ID").agg(
                    pl.col("ALLELE_SYMBOL").unique().sort().str.join(", ").alias("alleles")
                ),
                on="STRAIN/STOCK_ID",
                how="left",
            )
            .with_columns(
                pl.col("OTHER_NAMES").str.extract(r"(RRID:MMRRC_[\w.-]+)").alias("rrid"),
                pl.col("PUBMED_IDS").str.extract_all(r"\d{6,9}").alias("_pmids"),
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


if __name__ == "__main__":
    app.run()
