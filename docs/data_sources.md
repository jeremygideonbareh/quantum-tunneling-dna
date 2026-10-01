# Data sources

Status as of 2026-10-01. ✔ = downloaded and processed in this repo · ☐ = still needed.

| Dataset | Status | What we use | How we get it |
|---|---|---|---|
| **Jago et al. PNAS 2026 compilation** (Lagator group) | ✔ | 120,208 *E. coli* MA base substitutions with repair group, strand and ±10 bp motif. Includes **Foster 2018** (MMR⁻, 18 lines), **Niccum 2018** (mutD5 ± mutL), Lee 2012, Foster 2015, Long 2016, Tincher 2017 | `git clone https://github.com/Lagator-Group/extended-sequence-context` (GPL-3). Also ships the MG1655 genome (NC_000913.3) |
| *E. coli* MG1655 genome | ✔ | Opportunity counts | In the repo above |
| **COSMIC SBS v3.4, GRCh38** | ✔ | SBS6/14/15/20/21/26/44 (MMR), SBS10a–d (POLE/POLD1), SBS1/SBS5 controls | Bundled in the PyPI package `SigProfilerAssignment` (`data/Reference_Signatures/GRCh38/`) |
| GRCh38 and yeast trinucleotide (96) and pentanucleotide (1536) counts | ✔ | Per-site normalisation | Bundled in PyPI `SigProfilerMatrixGenerator` (`references/chromosomes/context_distributions/`) |
| **Al-Hashimi conformational fingerprints** | ✔ | G•T⁻ anion (16 T-centred contexts), A:T and G:C Hoogsteen, A:T base opening | `git clone https://github.com/alhashimilab/DNA-conformational-fingerprinting` (`Mutational-Signatures-JSD/JSD-python/Input-files/`). No licence file, so cite the paper and Zenodo 10.5281/zenodo.18459447 and don't redistribute |
| **DNAkmerQM** | ✔ downloaded | Baseline QM features (energies, Mulliken/ESP charges, Curves+ geometry) for all 6/7-mers, B/A/Z | `git clone https://github.com/SahakyanLab/DNAkmerQM` (CC-BY-4.0; 220 MB of zips) |
| Lujan et al. 2014 Genome Res 24:1751 | ☐ | Yeast pol variants × msh2Δ, about 40,000 mutations | Supplementary tables (journal site blocked here; see `needs_your_login.md` #4) |
| Zou et al. 2021 Nat Cancer | ☐ | Human iPSC MLH1/MSH2/MSH6/PMS2 knockouts | Check the data availability statement; raw data may need an EGA request |
| Hasenauer et al. 2025 NAR 53(2):gkae1196 | ☐ | MutS/MutL ChIP-seq mismatch hotspots (pre-repair) | GEO accession in the paper |
| Garushyants et al. 2024 GBE 16(4):evae035 | ☐ (new) | Wild-type *E. coli* natural mutation spectra (polymerase-error dominated) | Supplementary data |
| Roulette (Seplyarskiy 2023) | ☐ optional | Human germline per-site rates | genetics.bwh.harvard.edu/downloads/Vova/Roulette/ (large; blocked here) |

## Processed tables (`data/processed/`, committed)

- `ecoli_rates_k3.csv`, `ecoli_rates_k5.csv`: per group (`MMR_def`, `MMR_def_pol_exo`, `pol_exo`, `WT_like`) × context: `n_mut, n_sites, rate, rate_norm, cpg`
- `cosmic_sitewise_k3.csv`: per signature × context: `prob, n_sites, rate, rate_norm`

Regenerate with `bash scripts/download_data.sh && python scripts/ecoli_context_rates.py && python scripts/cosmic_sitewise_rates.py`.
