import marimo

__generated_with = "0.23.15"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import polars as pl
    from pathlib import Path

    return Path, mo, pl


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


if __name__ == "__main__":
    app.run()
