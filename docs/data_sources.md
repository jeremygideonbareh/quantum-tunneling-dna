# Data sources

| Dataset | What we use | Where | Script / status |
|---|---|---|---|
| COSMIC SBS v3.4 (GRCh38) | SBS6, 14, 15, 20, 21, 26, 44 (MMR); SBS10a–d (POLE/POLD1); SBS1 as CpG control | cancer.sanger.ac.uk/signatures/downloads | `download_data.sh` |
| DNAkmerQM | Baseline QM features, 7-mers → aggregate to 3/5-mers | github.com/SahakyanLab/DNAkmerQM | `download_data.sh` |
| Roulette | Human germline per-site rates (secondary) | genetics.bwh.harvard.edu/downloads/Vova/Roulette/ | manual (large) |
| E. coli MG1655 genome | Opportunity counts | NCBI GCF_000005845.2 | `download_data.sh` |
| S. cerevisiae R64 genome | Opportunity counts | NCBI GCF_000146045.2 | `download_data.sh` |
| Foster 2018 (Genetics) | E. coli mutS/mutL MA mutations | supplementary tables | ☐ Nazia: file names |
| Lagator lab PNAS 2026 | 100k+ E. coli mutations, 32 MA lines | data availability statement | ☐ Nazia |
| Hasenauer 2025 (NAR) | MutS/MutL ChIP-seq mismatch hotspots | GEO accession in paper | ☐ Nazia |
| Lujan 2014 (Genome Res) | Yeast pol variant × msh2Δ mutations | supplementary tables | ☐ Nazia |
| Zou 2021 (Nat Cancer) | Human iPSC MMR-knockout mutations | paper's data links | ☐ Nazia |
| Al-Hashimi 2026 | 16-context G•T⁻ propensities | supplementary / source data | ☐ Nazia: pull into `data/processed/alhashimi_gt.csv` |

**Format we convert everything to:** `chrom, pos (1-based), ref, alt, sample, strand_info` TSV in `data/processed/`. Then run `qtdna.mutations.count_transitions(df, genome, k)` to get per-context `n_mut / n_sites / rate`.
