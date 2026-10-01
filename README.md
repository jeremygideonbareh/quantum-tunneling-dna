# Quantum tunnelling fingerprints of DNA sequence context

**Question:** does the neighbouring sequence change how easily protons tunnel inside a base pair, and does that predict where replication-error mutations happen?

Jeremy Bareh (computation) · Nazia (biology/research). Target: bioRxiv by March 2027.

- **Plan (read first):** [`docs/research_plan.md`](docs/research_plan.md): literature check, five changes to v1, week-by-week tasks
- Pre-registration draft: [`docs/preregistration.md`](docs/preregistration.md)
- Data: [`docs/data_sources.md`](docs/data_sources.md) · Papers: [`docs/literature.md`](docs/literature.md)

## Setup
```bash
conda env create -f environment.yml && conda activate qtdna
pip install -e ".[dev]"
pytest
```
`pip install -e ".[dev]"` alone also works (it installs `tblite` for GFN2-xTB). ORCA and 3DNA need free academic registration.

## What runs today
```bash
python scripts/smoke_test_gc_scan.py --out results/smoke_gc   # ~1-2 min on 4 cores
bash scripts/download_data.sh                                  # COSMIC, DNAkmerQM, genomes
```
The smoke test builds a methyl-capped G:C pair, optimises G:C and G\*:C\*, scans the double proton transfer, and reports the barrier, ZPE and WKB tunnelling factor. Result so far: G\*:C\* at +10.1 kcal/mol, barrier 12.5 kcal/mol, reverse barrier 2.5 kcal/mol (GFN2-xTB, gas phase; a pipeline check, not a result).

## Layout
```
src/qtdna/contexts.py    32 trinucleotide / 512 pentanucleotide contexts, CpG flags, SBS-96 labels
src/qtdna/basepair.py    builds a methyl-capped Watson-Crick G:C pair (no external tools)
src/qtdna/tunneling.py   Wigner, Skodje-Truhlar, WKB tunnelling factors
src/qtdna/mutations.py   per-context transition counts / opportunity -> per-site rate
scripts/                 smoke test, data download
results/smoke_gc/        first calculation outputs
```
