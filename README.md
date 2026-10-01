# Quantum tunnelling fingerprints of DNA sequence context

**Question:** does the neighbouring sequence change how easily protons tunnel inside a base pair, and does that predict where replication-error mutations happen?

Jeremy Bareh (computation) · Nazia (biology/research). Target: bioRxiv by March 2027.

**Start here:** [`docs/research_plan.md`](docs/research_plan.md) (plan and status) · [`docs/needs_your_login.md`](docs/needs_your_login.md) (manual to-dos)

| Doc | What it is |
|---|---|
| [`known_vs_missing.md`](docs/known_vs_missing.md) | Literature summary / intro draft |
| [`preregistration.md`](docs/preregistration.md) | Hypotheses to post on OSF **before** fingerprint-vs-mutation tests |
| [`data_sources.md`](docs/data_sources.md) | Every dataset, its status and how to fetch it |
| [`outreach_emails.md`](docs/outreach_emails.md), [`one_page_summary.md`](docs/one_page_summary.md) | October emails + attachment |
| [`literature.md`](docs/literature.md) | Annotated references |

## Setup
```bash
conda env create -f environment.yml && conda activate qtdna && pip install -e ".[dev,dft]"
python -m venv .pymolenv && .pymolenv/bin/pip install pymol-open-source-whl "numpy<2"   # B-DNA builder
pytest
```

## Pipeline
```bash
bash scripts/download_data.sh                         # E. coli MA data, COSMIC, Al-Hashimi, DNAkmerQM (GitHub/PyPI only)
python scripts/ecoli_context_rates.py                 # -> data/processed/ecoli_rates_k{3,5}.csv
python scripts/cosmic_sitewise_rates.py               # -> data/processed/cosmic_sitewise_k3.csv
.pymolenv/bin/python scripts/build_bdna_clusters.py --k 3   # -> structures/k3/*.pdb (32 capped trimers)
python scripts/run_fingerprint.py --k 3 --env gas     # xTB DPT scans, ~5 min/context, resumable
python scripts/run_fingerprint.py --k 3 --env water   # same in ALPB water
python scripts/dft_check_scan.py results/smoke_gc     # DFT (PySCF) on an xTB scan
python scripts/smoke_test_gc_scan.py                  # single G:C pair sanity test
```

## Layout
```
src/qtdna/contexts.py    32 / 512 pyrimidine-centred contexts, CpG flags, SBS-96 labels
src/qtdna/mutations.py   per-context transition counts / genome opportunity -> per-site rate
src/qtdna/scan.py        DPT scan of the central pair in a frozen B-DNA cluster -> fingerprint
src/qtdna/tunneling.py   Wigner, Skodje-Truhlar, WKB tunnelling factors
src/qtdna/basepair.py    standalone methyl-capped G:C pair (smoke test)
structures/k3/           32 capped 3-bp B-DNA clusters (PyMOL fnab)
data/processed/          per-site rate tables (E. coli, COSMIC)
results/                 smoke test, fingerprint tables, DFT checks
```
