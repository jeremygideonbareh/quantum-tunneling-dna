# What is known, what is missing (draft for Nazia to edit into the Introduction)

*Draft written 2026-10-01 from abstracts, public repositories and the data we downloaded. Nazia: please check each claim against the full text and mark ✔ when verified.*

## 1. Replication errors have strong, reproducible sequence-context preferences

- In mismatch-repair (MMR) deficient *E. coli*, per-site transition rates vary **54-fold across T-centred trinucleotides and 18-fold across non-CpG C-centred ones** (our computation from the 120,208-mutation compilation of Jago et al., PNAS 2026; `data/processed/ecoli_rates_k3.csv`). A 5′ G raises both C>T and T>C (GTA, GTC, GCC are the top contexts).
- Jago et al. 2026 show context effects reach ±6 bp and, in some cases, about 1 kb. Which contexts matter depends on which repair pathways are active and on leading vs lagging strand.
- Hasenauer et al. (NAR 2025) mapped mismatches *before* repair (MutS/MutL ChIP-seq). Hotspots sit in low-stability DNA, mononucleotide repeats and secondary-structure-prone sequence.
- Garushyants, …, Agashe, …, Gelfand (GBE 2024): even in wild-type *E. coli*, natural and lab mutations are explained mainly by **polymerase errors**, not repair defects.
- In humans, MMR-deficiency signatures (SBS6/15/20/21/26/44) are produced by differential misincorporation by the replicative polymerases (Zou et al., Nat Cancer 2021, iPSC knockouts).
- **Cross-species consistency (our exploratory check):** the *E. coli* MMR-deficient C>T context profile correlates with human MMR signatures after per-site normalisation (Spearman ρ = 0.80 SBS21, 0.69 SBS26, 0.65 SBS15, 0.59 SBS6). Context preferences shared by bacteria and human tumours suggest a cause intrinsic to DNA and polymerase chemistry, not to organism-specific biology.

## 2. Candidate physical explanations, and how far each has been tested

| Mechanism | Evidence for | Evidence against / gaps |
|---|---|---|
| **Rare tautomers via proton transfer (Löwdin)** | Slocombe et al. 2022 (Commun Phys): tunnelling dominates over classical hopping, tautomer occupation about 1.7×10⁻⁴. Slocombe 2022 (Commun Chem): strand separation can trap tautomers. Slocombe 2023 (JPCL): tunnelling in G–T wobble ⇄ tautomer | Soler-Polo et al. 2019 (JCTC, QM/MM): G\*:C\* unstable in B-DNA. **Gheorghiu, Coveney & Arabi 2020** (Interface Focus 10(6), 2020): G\*:C\* half-life is only picoseconds against milliseconds for unwinding, so double proton transfer "has a negligible contribution". JOC 2025: zero-point energy can erase the tautomer well. **No study has computed how this varies with sequence context across all contexts, or tested it against mutation data.** |
| **Anionic / tautomeric Watson–Crick-like G•T mismatches** | Kimsey et al. 2018 (Nature); Szekely, …, Al-Hashimi 2026 (Nat Commun): G•T⁻ population varies 50-fold across 16 contexts, driven by the 3′ neighbour | Their best COSMIC match was **SBS11 (temozolomide)**, not MMR or replication signatures. In our exploratory check, G•T⁻ does **not** predict MMR-deficient T>C rates (E. coli ρ = +0.32, p = 0.23; SBS6 ρ = −0.66; SBS44 ρ = −0.70). Measured in free duplexes, not in a polymerase |
| Generic QM descriptors of k-mers | Masuda & Sahakyan 2024 (Sci Data): DNAkmerQM features for all 7-mers predict A>C rates | Descriptive, not mechanistic: no reaction barriers or tautomer energetics |
| Duplex stability / breathing | Hasenauer 2025 hotspots enriched in low thermal stability | Not reaction-specific; must be a covariate in any test |

## 3. The gap this project fills

1. **No sequence-resolved map of proton-transfer energetics exists.** Previous calculations treat one isolated pair, one duplex (Gheorghiu 2020, a dodecamer) or a few stacked pairs. Our 32-context (then 512-context) fingerprint is new.
2. **The tautomer hypothesis has never been tested against mutation data.** The Surrey and Coveney groups argue the mechanism from physics alone; Al-Hashimi tests a different rare state (G•T⁻) against cancer signatures. We test the canonical-pair tautomer against replication-error data in three species with pre-registered hypotheses.
3. **A null result is informative.** If the fingerprint fails to predict replication-error context preferences, while those preferences are strong and conserved, that is quantitative evidence against canonical-pair tautomerism as a major source of context-dependent replication errors. It would address the Slocombe vs. Soler-Polo/Gheorghiu debate with data.

## 4. Points the Discussion must handle

- **Which errors?** G\*:C\* gives G\*•T and A•C\* mispairs, i.e. C>T on the pyrimidine strand; A\*:T\* gives T>C. Proofreading (mutD5 / POLE) and MMR remove different subsets, so compare the MMR-only and MMR+proofreading groups separately.
- **Polymerase context.** Free-duplex energetics may not reflect the active site. Discuss this, and add the G•T wobble → WC-like variant if time allows.
- **CpG.** Methyl-C deamination dominates CpG C>T in humans (SBS1) but not in *E. coli*, which lacks CpG methylation. *E. coli* lets us keep CpG contexts as a secondary test.
- **Strand asymmetry.** Leading vs lagging strand is available in the E. coli data (`strand` column), which allows a direction-specific test.
