# Pre-registration (draft, Nazia to finish; post on OSF before any correlation is run)

**Primary fingerprint:** κ-weighted tautomer population P = κ_sudden · exp(−ΔG‡/kT), per context, from 3-bp clusters, xTB with ALPB water. The DFT-validated version replaces it if rank ρ(xTB, DFT) < 0.7.

## Hypotheses
- **H1 (primary).** Across the 12 non-CpG C-centred trinucleotides, P correlates positively with the per-site C>T rate in MMR-deficient E. coli (Foster 2018 + Lagator 2026). One-sided Spearman, exact permutation p, α = 0.05.
- **H2.** Same for the 16 T-centred contexts vs T>C.
- **H3.** Replication: H1/H2 hold in yeast msh2Δ (Lujan 2014) and human MMR-KO (Zou 2021). Stouffer meta-analysis across species.
- **H4 (added value).** At the 5-mer level, adding P to a baseline (local stability + DNAkmerQM PC1–3) improves leave-one-out R² by more than 0.02, and the coefficient on P is positive.
- **H5 (mechanism check).** The fingerprint varies more with the 3′ neighbour than the 5′ neighbour (as Al-Hashimi found for G•T⁻).

## Fixed in advance
- CpG contexts are excluded from C>T tests. Proofreading-defective signatures (SBS10a–d) are analysed separately.
- Multiple testing: Holm across H1–H3 × datasets.
- Rates are per-site (divided by genome opportunity), never raw counts or signature weights.
- A null result is reported with confidence intervals on ρ.
