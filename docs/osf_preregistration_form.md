# OSF Preregistration, ready to paste

Use osf.io → Registries → **OSF Preregistration** template. Each heading below is one field in that form, in order. Paste the text under it.

Final version: H2 = option A (chosen by the authors, 1 Oct 2026).

**How to post (about 15 minutes):**
1. Go to osf.io, sign up or log in, and click **Create new project**. Title: the title below. Add Nazia as a contributor (Contributors tab).
2. In the project, go to **Registrations** → **New registration** → template **"OSF Preregistration"**.
3. Paste each section below into the field with the same name. Leave fields not listed here blank.
4. Under files or supplements, attach `docs/preregistration.md` and link the GitHub repo.
5. Choose **"Make registration public immediately"** (or an embargo of up to 4 years if you'd rather keep it private until the preprint). Then click **Register**. Nazia must approve it by email if she is a contributor.
6. Copy the registration URL (osf.io/xxxxx) and send it to Claude.

---

## Title
Does sequence-context dependence of base-pair proton-transfer energetics predict replication-error mutation spectra? A pre-registered test across *E. coli*, yeast and human data

## Description
Replication errors (C>T and T>C transitions) occur at very different rates in different sequence contexts. Löwdin's hypothesis proposes that proton transfer within Watson–Crick base pairs creates rare tautomers that mispair during replication. We computed a "proton-transfer fingerprint" for each of the 32 pyrimidine-centred trinucleotide contexts: semiempirical (GFN2-xTB) double-proton-transfer scans of the central pair in capped 3-bp B-DNA clusters, with WKB tunnelling corrections. We will test whether this fingerprint predicts per-site transition rates in mismatch-repair-deficient mutation data. A null result will be reported.

## Hypotheses
- **H1 (primary).** Across the 16 C-centred contexts, the fingerprint F correlates positively with the per-site C>T rate in *E. coli* MMR-deficient, proofreading-proficient mutation-accumulation lines.
- **H2.** Across the 16 T-centred contexts, the energy at the half-transfer point of the A:T proton-transfer scan (X = 0, the single-proton-transfer region) correlates **negatively** with the per-site T>C rate in the same *E. coli* lines. Rationale: in our calculations A\*:T\* is never a stable double-proton-transfer tautomer, so the test targets the route the A:T pair can actually take.
- **H3 (replication).** The H1 and H2 relationships also hold in yeast *msh2Δ* lines (Lujan et al. 2014) and in human MMR-deficiency signatures SBS6, SBS15, SBS21, SBS26 and SBS44 (COSMIC v3.4, normalised per site), combined by Stouffer's method.
- **H4 (no directional prediction).** Same tests in proofreading-deficient data (*E. coli* mutD5 lines; COSMIC SBS10a–d), reported separately.
- **H5 (added value, 5-mers).** Adding F to a baseline model (nearest-neighbour duplex stability, DNAkmerQM PC1–3, and Al-Hashimi G•T⁻ for T-centred contexts) improves leave-one-out R² by more than 0.02, with a positive coefficient.
- **H6.** The fingerprint varies more with the 3′ neighbour than with the 5′ neighbour (two-way ANOVA within each centre).

## Study type
Observational study: secondary analysis of existing mutation datasets plus new computational predictions.

## Blinding
No blinding is involved. Analysts have not compared any fingerprint value with any mutation outcome (see "Explanation of existing data").

## Study design
Correlational: computed per-context predictors (fingerprint) vs per-context per-site mutation rates. Units are sequence contexts (16 per centre at 3-mer level, 256 per centre at 5-mer level).

## Existing data
**Registration prior to analysis of the data.** The mutation data exist and have been processed into per-context rates. No analysis relating the predictor to the outcome has been run.

## Explanation of existing data
Before registration we (i) computed per-site context rates for all outcome datasets, (ii) correlated outcome datasets with each other (*E. coli* MMR-deficient C>T vs human SBS21: ρ = 0.80), and (iii) correlated the Al-Hashimi G•T⁻ baseline with the outcomes (*E. coli* T>C ρ = +0.32). We computed fingerprint values for all 32 contexts but did not compare them with any outcome. One fingerprint rule (κ = 1 when the tautomer well is shallower than RT) was added after seeing the fingerprint values and before any outcome comparison. The analysis code (`scripts/preregistered_tests.py`) was written beforehand and refuses to run without this registration's URL.

## Data collection procedures
Outcome data are public:
- Jago et al. PNAS 2026 compilation of 120,208 *E. coli* MA substitutions (github.com/Lagator-Group/extended-sequence-context).
- Lujan et al. 2014 yeast supplementary tables.
- COSMIC v3.4 GRCh38 SBS signatures.
- Zou et al. 2021 human iPSC MMR-knockout mutation calls (Mendeley doi 10.17632/ymn3ykkmyx).

Per-site rate = transitions in a context ÷ occurrences of that context in the reference genome, with strands collapsed (code: `qtdna.mutations`, `scripts/*_rates.py`). The predictor is computed by `scripts/run_fingerprint.py` (GFN2-xTB via tblite, ASE; structures from PyMOL fnab).

## Sample size
32 trinucleotide contexts (16 per centre); 512 pentanucleotide contexts for H5. Outcome counts: about 37,700 transitions in the *E. coli* MMR-deficient group.

## Sample size rationale
Fixed by the genetic code: all contexts are used. With n = 16, a one-sided Spearman test at α = 0.05 detects ρ ≳ 0.43; this power limitation is acknowledged. The 5-mer analysis adds power.

## Stopping rule
Not applicable (all contexts, all listed datasets).

## Manipulated variables
None.

## Measured variables
- **Predictor F** (C-centred, H1) = ln κ_sudden − ΔE‡/RT at 310 K, from the gas-phase xTB scan. If no tautomer well is deeper than RT, κ = 1.
- **Predictor for H2** (T-centred) = −E_half, the negated scan energy at X = 0 (so that, as for F, larger means more errors).
- Secondary predictors: tautomer energy, reverse barrier, the same quantities in ALPB water, and DFT-corrected values.
- **Outcomes:** per-site C>T (C-centred) and T>C (T-centred) rates per dataset.

## Indices
F as defined above. The xTB-based F is replaced by DFT-based F if the rank correlation between xTB and DFT values across contexts is below 0.7, a decision made before any outcome comparison.

## Statistical models
One-sided Spearman rank correlation, exact or 200,000-permutation p-values (H1–H4). Stouffer combination across datasets (H3). Nested linear models compared by leave-one-out R² (H5). Two-way ANOVA (H6).

## Transformations
Rates are used as-is (rank-based tests). F is on a natural-log scale. For H5, predictors are standardised.

## Inference criteria
α = 0.05, one-sided where a direction is stated. Holm correction across H1–H3 × datasets. Effect sizes are reported with 95% bootstrap CIs.

## Data exclusion
Human SBS analyses exclude the 4 CpG contexts for C>T (methyl-C deamination). *E. coli* keeps all contexts because it has no CpG methylation. No other exclusions.

## Missing data
Contexts with zero genome occurrences are not possible at the 3-mer level. At the 5-mer level, contexts with fewer than 100 genomic occurrences are excluded.

## Exploratory analysis (not confirmatory)
Leading- vs lagging-strand-specific rates; extended context (±6 bp); the G•T wobble tautomer route; strand-separated (stretched) pair models.

## Other
Code and processed data: github.com/jeremygideonbareh/quantum-tunneling-dna (Zenodo DOI at publication). Authors: Jeremy Bareh, Nazia Sooting.
