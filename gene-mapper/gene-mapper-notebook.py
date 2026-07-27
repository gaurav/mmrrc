import marimo

__generated_with = "0.23.15"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import polars as pl
    import re
    import xml.etree.ElementTree as ET
    from collections import defaultdict, deque
    from pathlib import Path
    from urllib.parse import quote_plus

    return ET, Path, defaultdict, deque, mo, pl, quote_plus, re


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

    **The MP ids in this column cannot be used as-is.** MMRRC's export splits labels
    that contain a comma, and from the first split onwards every remaining label in
    that strain's list is paired with the *next* term's id. Roughly one in seven
    annotations is attached to the wrong term. The cells below re-derive each id from
    its label instead; the highlights at the end of the section show the evidence.
    """)
    return


@app.cell
def _(ET, Path, defaultdict, deque, mo, pl, re):
    MP_OWL = Path("../data/mp.owl")
    MP_ROOT = "MP:0000001"

    _OBO = "http://purl.obolibrary.org/obo/"
    _RDF = "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}"
    _RDFS = "{http://www.w3.org/2000/01/rdf-schema#}"
    _OWL = "{http://www.w3.org/2002/07/owl#}"
    _OIO = "{http://www.geneontology.org/formats/oboInOwl#}"


    def mp_normalise(text):
        """Fold a phenotype label to a comparison key (case, punctuation, spacing)."""
        return re.sub(r"[^a-z0-9]+", " ", text.lower()).strip()


    mp_labels = {}
    mp_obsolete = set()
    mp_parents = defaultdict(list)
    mp_label_index = defaultdict(set)

    # One streaming pass over 101 MB of RDF/XML, ~1.5s -- no ontology library needed.
    # Only clear owl:Class elements: clearing every element wipes child text before
    # the parent's end event can read it.
    for _ev, _el in ET.iterparse(MP_OWL, events=("end",)):
        if _el.tag != _OWL + "Class":
            continue
        _about = _el.get(_RDF + "about", "")
        if _about.startswith(_OBO + "MP_"):
            _cid = "MP:" + _about.rsplit("MP_", 1)[1]
            _label = _el.findtext(_RDFS + "label")
            if _label:
                mp_labels[_cid] = _label
                mp_label_index[mp_normalise(_label)].add(_cid)
                for _syn in _el.findall(_OIO + "hasExactSynonym"):
                    if _syn.text:
                        mp_label_index[mp_normalise(_syn.text)].add(_cid)
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
    return (
        mp_categories,
        mp_label_index,
        mp_labels,
        mp_normalise,
        mp_obsolete,
        mp_rollup,
    )


@app.cell
def _(catalog, mo, mp_label_index, mp_labels, mp_normalise, pl, re):
    _PAIR = re.compile(r"^\s*(.*?)\s*\[(MP:\d+)\]\s*$")
    _MAX_JOIN = 4  # longest comma-split label seen in this data is 3 parts


    def mp_resolve(label):
        """The single MP term matching a label or exact synonym, else None."""
        _hits = mp_label_index.get(mp_normalise(label), set())
        return next(iter(_hits)) if len(_hits) == 1 else None


    def mp_is_fragment(label, mp_id):
        """True when `label` is only the text before the first comma of the term's label.

        Used to *locate* the first split for the diagnostic below. Not used for repair:
        it trusts the file's id, which is the half the split breaks.
        """
        _full = mp_labels.get(mp_id, "")
        return bool(_full) and "," in _full and label == _full.split(",")[0].strip() and label != _full


    # ponytail: repair for MMRRC's comma-split export. Any label containing a comma is
    # broken into separate entries, and from the first split onward each label is paired
    # with the *next* term's id. Repair is id-independent -- rejoin the longest run of
    # consecutive entries whose comma-joined label resolves to exactly one term, then
    # resolve by label -- because the ids are the broken half. Two exceptions: before any
    # split has occurred the file's id is still reliable, so a truncated label whose id
    # names "<label>, ..." keeps the id and its qualifier (penetrance, mostly); and a
    # label that resolves to nothing falls back to the id. See github.com/gaurav/mmrrc
    # issue #1. Delete all of this if MMRRC fixes the export -- mp_repair_stats will show
    # when it has stopped doing anything.
    _edges = []
    _st = dict(
        strains=0, split_strains=0, joined=0, entries=0,
        by_label=0, by_id_specific=0, by_id=0, unresolved=0, corrected=0,
        pre_ok=0, pre_n=0, post_ok=0, post_n=0,
    )

    for _row in (
        catalog.filter(pl.col("MPT_IDS").is_not_null())
        .unique(subset=["STRAIN/STOCK_ID"])
        .select(["STRAIN/STOCK_ID", "MPT_IDS"])
        .iter_rows(named=True)
    ):
        _m = [_PAIR.match(_p) for _p in _row["MPT_IDS"].split("|")]
        _pairs = [(_x.group(1), _x.group(2)) for _x in _m if _x]
        _st["strains"] += 1

        # Diagnostic: how often the file's label matches the id the file gives it,
        # before vs. after the first comma-split in this strain's list.
        _split_at = next(
            (_i for _i, (_l, _p) in enumerate(_pairs) if mp_is_fragment(_l, _p)), None
        )
        for _i, (_l, _p) in enumerate(_pairs):
            _ok = mp_labels.get(_p) == _l
            if _split_at is None or _i < _split_at:
                _st["pre_n"] += 1
                _st["pre_ok"] += _ok
            else:
                _st["post_n"] += 1
                _st["post_ok"] += _ok
        _st["split_strains"] += _split_at is not None

        _i, _shifted = 0, False
        while _i < len(_pairs):
            _take, _label = 1, _pairs[_i][0]
            for _k in range(min(_MAX_JOIN, len(_pairs) - _i), 1, -1):
                _cand = ", ".join(_p[0] for _p in _pairs[_i : _i + _k])
                if mp_resolve(_cand):
                    _take, _label = _k, _cand
                    break
            _given = _pairs[_i][1]
            _i += _take
            _st["joined"] += _take > 1
            _st["entries"] += 1

            _by_label = mp_resolve(_label)
            _full = mp_labels.get(_given, "")
            _truncated = (
                not _shifted
                and _take == 1
                and _by_label is not None
                and _by_label != _given
                and _full.startswith(_label + ",")
            )

            if _truncated:
                # The label lost its qualifier but the id is still aligned here.
                _final = _given
                _st["by_id_specific"] += 1
            elif _by_label is not None:
                _final = _by_label
                _st["by_label"] += 1
                if _take > 1 or _by_label != _given:
                    _shifted = True
            elif _given in mp_labels:
                _final = _given
                _st["by_id"] += 1
            else:
                _st["unresolved"] += 1
                continue

            _st["corrected"] += _final != _given
            _edges.append({"strain_id": _row["STRAIN/STOCK_ID"], "mp_id": _final})

    strain_phenotypes = pl.DataFrame(
        _edges, schema={"strain_id": pl.String, "mp_id": pl.String}
    ).unique()
    mp_repair_stats = _st

    mo.md(
        f"`strain_phenotypes`: **{strain_phenotypes.height:,} strain→phenotype edges** "
        f"over {strain_phenotypes['strain_id'].n_unique():,} strains and "
        f"{strain_phenotypes['mp_id'].n_unique():,} distinct MP terms. "
        f"Rejoined {_st['joined']:,} comma-split labels and re-pointed "
        f"**{_st['corrected']:,} of {_st['entries']:,}** entries "
        f"({_st['corrected'] / _st['entries']:.1%}) whose id disagreed with their label; "
        f"{_st['by_label'] / _st['entries']:.1%} resolved by label, "
        f"{_st['by_id_specific']:,} kept a more specific id, "
        f"{_st['unresolved']:,} unresolved."
    )
    return mp_repair_stats, strain_phenotypes


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
    mp_labels,
    mp_repair_stats,
    mp_rollup,
    phenotype_index,
    pl,
    strain_phenotypes,
):
    _s = mp_repair_stats
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

    **1. One in seven annotations points at the wrong term — this notebook repairs
    them.** MMRRC's export splits any label containing a comma, so `decreased
    CD4-positive, alpha-beta T cell number` arrives as two entries. From the first
    split onward, every remaining label in that strain's list is paired with the
    *next* term's id. The evidence is the collapse in label/id agreement either side
    of the split point:

    | Position in a strain's list | Label matches the id it was given |
    |---|---|
    | before the first split | {_s["pre_ok"]:,} of {_s["pre_n"]:,} ({_s["pre_ok"] / max(_s["pre_n"], 1):.1%}) |
    | **after** the first split | {_s["post_ok"]:,} of {_s["post_n"]:,} (**{_s["post_ok"] / max(_s["post_n"], 1):.1%}**) |

    It is an off-by-one chain, each label carrying the previous entry's id:

    | Label in the file | Id the file gives it | …which is really | Correct id |
    |---|---|---|---|
    | `increased pro-B cell number` | MP:0008547 | abnormal neocortex morphology | **MP:0008186** |
    | `abnormal neocortex morphology` | MP:0008869 | anovulation | **MP:0008547** |
    | `anovulation` | MP:0008882 | abnormal enterocyte physiology | **MP:0008869** |

    {_s["split_strains"]:,} of {_s["strains"]:,} annotated strains are affected.
    Rejoining the split labels and resolving each against the ontology re-pointed
    **{_s["corrected"]:,} of {_s["entries"]:,} entries
    ({_s["corrected"] / _s["entries"]:.1%})**. {_s["by_label"] / _s["entries"]:.1%}
    resolve by label or exact synonym; {_s["by_id_specific"]:,} keep the file's id
    because it is *more* specific than a truncated label (`prenatal lethality` where
    the file means `prenatal lethality, complete penetrance` — safe only before the
    first split, where the ids are still aligned); {_s["by_id"]:,} fall back to the
    file's id because the label no longer resolves at all (terms MP has retired), and
    {_s["unresolved"]:,} are unresolved.

    **2. MMRRC's labels are a stale snapshot.** Every id in the column is a real MP
    term, but the labels have drifted — `aggression towards males` is now `aggression
    towards male mice`, `reduced long term potentiation` is now `reduced long-term
    potentiation`, `altered response to myocardial infarction` is now `abnormal
    response to cardiac infarction`. Resolution therefore matches `rdfs:label` *and*
    `oboInOwl:hasExactSynonym`, and every table above shows the ontology's label
    rather than the catalog's.

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
