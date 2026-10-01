# Pre-registration (draft v2; Nazia to finish; post on OSF before any fingerprint-vs-mutation test)

## Fingerprint (fixed in advance)
Per pyrimidine-centred context, from capped 3-bp B-DNA clusters (flanking pairs and C1′ caps frozen), GFN2-xTB, relaxed scan along X = ξ₁ + ξ₂:

- **Primary:** forward tautomerisation rate proxy, `F = ln κ_sudden − ΔE‡ / RT` (T = 310 K).
- Secondary: (a) tautomer population `−ΔE / RT`; (b) reverse barrier ΔE‡_rev, i.e. tautomer survival; (c) the same three quantities in ALPB water.
- If the rank correlation between xTB and DFT values across contexts is below 0.7, the DFT values replace the xTB ones (decided before looking at mutation data).

## Outcomes
Per-site transition rate = transitions in context ÷ occurrences of that context in the genome (strands collapsed), computed by `qtdna.mutations` / `scripts/*_rates.py`.

## Hypotheses
- **H1 (primary).** *E. coli* MMR-deficient, proofreading-proficient lines (Jago et al. 2026 compilation, group `proofreading(+) MMR(-)`): across the 16 C-centred contexts, F correlates positively with the C>T per-site rate. *E. coli* lacks CpG methylation, so all 16 are kept. One-sided Spearman, exact permutation p, α = 0.05.
- **H2.** As H1 for the 16 T-centred contexts vs T>C.
- **H3 (replication).** H1/H2 in yeast msh2Δ (Lujan 2014) and human MMR signatures SBS6/15/21/26/44 (per-site normalised, the 4 CpG contexts excluded for C>T). Combine by Stouffer meta-analysis.
- **H4 (proofreading).** Repeat in proofreading-deficient groups (*E. coli* mutD5; SBS10a–d), reported separately with no directional prediction.
- **H5 (added value, 5-mers).** Adding F to a baseline (local GC/stability + DNAkmerQM PC1–3 + Al-Hashimi G•T⁻ for T-centred contexts) improves leave-one-out R² by more than 0.02 with a positive coefficient.
- **H6 (mechanism check).** The 3′ neighbour explains more of the variance in F than the 5′ neighbour (two-way ANOVA on the 16 contexts per centre).

Holm correction across H1–H3 × datasets. A null result is reported with 95% CIs on ρ.

## Analyses already run before registration (disclosed for transparency)
On 2026-10-01 we computed (i) per-site context rates for all datasets, (ii) correlations **between datasets** (E. coli vs COSMIC) and (iii) correlations between the **Al-Hashimi G•T⁻ baseline** and the outcomes. We did **not** correlate any proton-transfer fingerprint with any mutation outcome. Fingerprint values for some contexts were computed, but they were not compared with mutation data.
