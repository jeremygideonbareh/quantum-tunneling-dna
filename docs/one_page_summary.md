# Does DNA sequence context change proton-tunnelling risk, and does that predict replication errors?

**Jeremy Bareh** (computation) · **Nazia [surname]** (biology) · October 2026 · Bengaluru

**Question.** Replication errors (C>T and T>C transitions) are far more frequent in some sequence contexts than others. In mismatch-repair-deficient *E. coli*, per-site rates vary up to **54-fold** across trinucleotides, and the same context preferences appear in human MMR-deficient tumours (Spearman ρ up to 0.8 with SBS21). Löwdin's hypothesis says proton tunnelling inside base pairs creates rare tautomers that mispair. **Does the energetics of that tunnelling depend on neighbouring bases strongly enough to explain these context preferences?**

**Why now.** Physics papers argue for (Slocombe et al. 2022, 2023) and against (Soler-Polo 2019; Gheorghiu, Coveney 2020) a role for tautomers, but none tested sequence-resolved predictions against mutation data. Szekely, …, Al-Hashimi (Nat Commun 2026) fingerprinted a different rare state (anionic G•T⁻) across 16 contexts. It matched a chemotherapy signature (SBS11), not replication-error signatures.

**Approach.**
1. **Fingerprint:** for all 32 trinucleotide contexts, a 3-bp methyl-capped B-DNA cluster with the flanking pairs frozen. A relaxed double-proton-transfer scan (GFN2-xTB) gives the barrier, tautomer energy, reverse barrier, and WKB tunnelling factors in adiabatic and sudden limits. Validation: DFT (B3LYP-D3, PySCF) on the same paths. Extension: all 512 five-base contexts.
2. **Outcomes:** per-site transition rates (counts ÷ genome occurrences) from 120,208 *E. coli* MA mutations (Jago et al. PNAS 2026; Foster 2018; Niccum 2018), yeast msh2Δ (Lujan 2014), human MMR knockouts (Zou 2021) and COSMIC MMR/POLE signatures.
3. **Tests (pre-registered on OSF before any fingerprint–mutation comparison):** within-centre Spearman correlations, cross-species meta-analysis, and nested models against baselines (Al-Hashimi G•T⁻, DNAkmerQM QM features, duplex stability).

**Status (1 Oct 2026).** Pipeline running end to end: B-DNA cluster builder, xTB scans, WKB, DFT check, per-site rate tables for *E. coli* and COSMIC. Single G:C pair (xTB): G\*:C\* at +10 kcal/mol with only a 2–3 kcal/mol reverse barrier.

**Outputs.** bioRxiv preprint by March 2027; code and data on GitHub and Zenodo. A null result is publishable: it would be the first quantitative test of the tautomer hypothesis against mutation data.

**What we'd value from you:** feedback on the method or confounders, a DFT sanity check, or a modest compute allocation.
