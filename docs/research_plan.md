# Research plan v2: what we checked and what changes

**Authors:** Jeremy Bareh (computation) · Nazia (biology/research)
**Version:** 2026-10-01. Revises the 5-page plan of the same date (`docs/original_plan.md` holds the summary).

**Question.** Does the neighbouring sequence change how easily protons tunnel inside a base pair, and does that predict where replication-error mutations happen?

The question holds up and nobody has published this exact test. The literature check below turned up **five changes** to make before Phase 2. Each one is cheap now and expensive later.

---

## 1. What the literature check found

| Paper | What it means for us |
|---|---|
| **Szekely, …, Al-Hashimi**, *Nat Commun* 17:5615 (Apr 2026) | Used ¹⁹F NMR to measure the anionic Watson–Crick-like G•T⁻ state in all 16 triplet contexts. It varies up to **50-fold**, driven mostly by the **3′ neighbour**. They compared it with COSMIC and found plausible links to damage and therapy signatures. This is the closest competitor. Our angle has to stay distinct: we compute **canonical-pair proton transfer (G\*:C\*, A\*:T\*)** in all 32 contexts and test it on **replication-error data from MMR-deficient cells across three species**. Their 3′-neighbour result gives us a falsifiable check: does our fingerprint also depend mostly on the 3′ base? |
| **Masuda & Sahakyan**, *Sci Data* 2024 (DNAkmerQM) | QM features are already available for **every 7-mer** in B, A and Z DNA (CC-BY-4.0, on GitHub and Zenodo). We should use them as the baseline and not recompute generic QM descriptors. They also built an A>C mutation-rate model from them, so "QM features predict mutation" alone is not new. Our added value is the **reaction-specific** barrier, tautomer stability and tunnelling quantities. |
| **Lagator lab**, *PNAS* 123:e2601345123 (Jun 2026) | More than **100,000 E. coli mutations** from 32 MA experiments with varying proofreading and MMR. Context effects reach 6 bp and even about 1 kb. This is probably the best single E. coli dataset for us, and it was marked "not yet read" in v1. **Nazia should read it first.** |
| **Hasenauer et al.**, *NAR* 53(2):gkae1196 (2025) | MutS/MutL ChIP-seq maps of **mismatches before repair**, which is closer to the raw replication error than mutation calls are. Hotspots are enriched for **low thermal stability**, mononucleotide repeats and secondary structure, so local stability becomes a required covariate. (v1 cited this paper without volume and issue.) |
| **Soler-Polo et al.**, *JCTC* 15(12) (2019) | QM/MM free energies in B-DNA show **G\*:C\* is unstable under cellular conditions**. The main argument against our hypothesis is that the tautomer may not survive long enough to matter. Reviewers will cite it, so we must address it directly (Change 1). |
| **Slocombe et al.**, *Commun Phys* 2022 and *Commun Chem* 2022 | Open-quantum-systems tunnelling: tautomer occupation of about 1.7×10⁻⁴. Strand separation can trap tautomers. This is the counter-argument to Soler-Polo. |
| **Slocombe et al.**, *JPCL* 2023 | G–T wobble ⇄ tautomer misincorporation inside the polymerase. Tunnelling dominates G\*T ⇄ GT\*. |
| PCCP 2026 (*Solvation reshapes the classical but not the quantum mechanism…*) and JOC 2025 (*ZPE quantum barriers suppress double-proton tunnelling*) | Two points. (a) The tunnelling path is not the classical minimum-energy path. (b) Zero-point energy can erase the shallow tautomer well. Both affect how we define the fingerprint. |

Links are in `docs/literature.md`. Nature, PMC and Zenodo were blocked from this session, so I could not read full texts. **Nazia, please confirm two things from the Al-Hashimi paper and its supplement:** (i) whether they tested the MMR-deficient signatures (SBS6/14/15/20/21/26/44) specifically, and (ii) whether their 16 context values are in a supplementary table.

---

## 2. The five changes

### Change 1: the fingerprint is a vector, not one number
For each context, report:

| Feature | Why |
|---|---|
| ΔE(G:C → G\*:C\*) and forward/reverse barrier | Core energetics |
| **ZPE-corrected reverse barrier** and whether the tautomer is still a minimum | Answers Soler-Polo: if the well disappears, the tautomer lifetime is femtoseconds |
| Tunnelling factor κ in two limits: **adiabatic** (whole frame relaxes) and **sudden** (protons only, frame frozen) | These bracket the real answer. Slocombe-type tunnelling sits near the sudden limit. Our smoke test already shows a 2.5× spread between them |
| Equilibrium tautomer population exp(−ΔG/kT) × κ-weighted rate | The single number we correlate first (**pre-registered as primary**) |
| *(stretch goal)* Same quantities with the pair pulled 0.5 Å apart | Strand-separation regime (Slocombe 2022) |
| *(stretch goal)* G•T wobble → WC-like tautomer in context | Direct comparison with Al-Hashimi's 16 values |

### Change 2: fix the model before scaling to 32 contexts
- **Cluster:** 3 stacked base pairs. Sugars are replaced by methyl caps at N9/N1. The cap carbons and the flanking pairs' heavy atoms are frozen at B-DNA positions (from 3DNA `fiber`). Only the central pair relaxes. Without the freeze, xTB will unstack the trimer.
- **Environment:** run each scan twice, gas phase and ALPB water (xtb CLI `--alpb water`). Report both. A fingerprint that flips between them is not robust.
- **Reaction coordinate:** constrain **X = ξ₁ + ξ₂**, where ξ = d(D–H) − d(H···A). This lets the path be stepwise or concerted. *We learned this the hard way:* a plain N–H bond-length scan pulled the acceptor away and gave a meaningless 90 kcal/mol ramp. A concerted ξ₁ = ξ₂ scan gave a broken profile. The summed coordinate gives a smooth curve (section 4).
- **DFT validation, increased from 5 contexts to all 32 at stationary points:** reactant, tautomer and scan maximum (96 single points), plus full DFT scans for 5 contexts. Suggested level: ωB97X-D3/def2-TZVP single points on xTB geometries, with B3LYP-D3 as a second functional. Report the **rank correlation of DFT vs xTB across contexts**, not only absolute errors. Rankings are what feed the statistics.
- **Possible shortcut:** a universal ML potential trained on ωB97M-V data (e.g., Meta's UMA/OMol25 models) might give near-DFT scans at xTB cost. Check the licence and test it on the 5 DFT contexts before relying on it.

### Change 3: compute the 5-mers directly; don't extrapolate with ML from 32 points
There are **512** pyrimidine-centred 5-mers, not about 1,000. Training ML on 32 points and predicting 512 is extrapolation, and reviewers will reject it. xTB scans of a 5-bp cluster take minutes each, so all 512 can run in roughly a week on a workstation. The ML angle becomes: *learn the fingerprint from sequence plus DNAkmerQM features, check it against the 512 computed values, and then use it on 7-mers.* That framing also works for an AI-for-Science workshop paper.

### Change 4: statistics designed for n = 16
- Correlate **within each centre**: C-centred contexts vs C>T, T-centred contexts vs T>C. Mixing centres builds in a trivial G:C vs A:T difference.
- Exclude the 4 CpG contexts from C>T in the primary test (deamination). That leaves **n = 12**. With n = 12, Spearman ρ must reach about 0.58 for p < 0.05, so state up front that only a **large** effect is detectable at the 3-mer level. The 5-mer analysis (n = 256 per centre) carries the statistical power.
- Use exact permutation p-values. Combine E. coli, yeast and human by meta-analysis (Fisher or Stouffer), not by pooling counts.
- Test the "beyond baselines" question (RQ3) with **nested models** compared by leave-one-out CV: baseline (local GC/stability + DNAkmerQM PC1–3, + Al-Hashimi value) vs baseline + fingerprint. With 12–16 points, allow at most 1–2 baseline predictors.
- Use **per-site rates, not signature weights.** COSMIC SBS profiles reflect genome trinucleotide frequency, so divide by opportunity (`qtdna.mutations` does this).

### Change 5: fix the dataset list
- MMR-deficient COSMIC set: **SBS6, 14, 15, 20, 21, 26, 44**. v1 was missing SBS14 (POLE + MMR). Polymerase-proofreading set: **SBS10a/10b** (POLE) and **SBS10c/10d** (POLD1). Treat proofreading signatures as a *separate* hypothesis, because exonuclease defects change which errors survive.
- Add the Lagator PNAS 2026 compilation and Hasenauer 2025 ChIP-seq (mismatches before repair).
- Leading vs lagging strand: E. coli and yeast data allow strand-resolved spectra. Keep this as a secondary analysis, since tautomer-induced errors should depend on which base is in the template.

---

## 3. Updated timeline (today = week 1)

### Phase 1, weeks 1–4 (October)
**Jeremy**
- [x] Repo, Python environment, GFN2-xTB via `tblite` (2026-10-01)
- [x] First end-to-end calculation: single G:C pair, DPT scan, ZPE, WKB κ (`scripts/smoke_test_gc_scan.py`)
- [x] Context enumeration, per-site mutation-rate counter, stability baseline, with unit tests
- [x] ~~Register for ORCA and 3DNA~~ **No longer needed:** PySCF (pip) replaces ORCA and PyMOL `fnab` replaces 3DNA
- [x] Capped 3-bp B-DNA clusters for all 32 contexts (`structures/k3/`), frozen-frame DPT scans for both G:C and A:T centres (`scripts/run_fingerprint.py`)
- [x] Data: E. coli (Jago 2026, 120k mutations), COSMIC v3.4 + GRCh38/yeast context counts, Al-Hashimi fingerprints, DNAkmerQM, all fetched by `scripts/download_data.sh`
- [x] Per-site rate tables: `data/processed/ecoli_rates_k{3,5}.csv`, `cosmic_sitewise_k3.csv`; baselines `dnakmerqm_B_k{3,5}.csv`
- [x] DFT (B3LYP-D3/def2-SVP, PySCF) on the smoke-test scan path (`scripts/dft_check_scan.py`)
- [x] Pre-registered analysis code written and locked behind the OSF URL (`scripts/preregistered_tests.py`)
- [ ] Finish the 32-context gas-phase run, then repeat in ALPB water (`--env water`)
- [ ] DFT single points at reactant / maximum / tautomer for all 32 clusters (about 100 atoms each; needs the HPC allocation or about 2 days locally at def2-SVP)
- [ ] 5-mer clusters (512) after the HPC reply
- [ ] Tighten tautomer convergence (2 imaginary modes in the single-pair smoke test)

**Nazia**
- [ ] Read Jago 2026 first, then Al-Hashimi 2026 (confirm the MMR-signature question in `literature.md`), Gheorghiu 2020, Soler-Polo 2019 and Slocombe ×3
- [ ] Edit `docs/known_vs_missing.md` (drafted) into the 2-page summary; verify each claim
- [ ] Get the Lujan 2014 yeast tables and the Zou 2021 data link (`needs_your_login.md` #4–5)
- [ ] Confounder rules: CpG, strand, local stability or GC, repeats, transcription
- [ ] **Finish `docs/preregistration.md` (v2 drafted) and post it on OSF** before anyone runs `preregistered_tests.py`

**Handoff, end of October:** dataset list + confounder rules + pre-registration → Jeremy.

### Phase 2, weeks 5–10 (Nov to mid-Dec)
32 × 3-bp clusters × {gas, ALPB} with xTB. DFT at stationary points for all 32 and full DFT scans for 5. Mutation counting for all datasets. **Handoff:** fingerprint table (CSV + `results/fingerprint_v1.md`) → Nazia sanity-checks it, including the 3′-neighbour check against Al-Hashimi.

### Phase 3, weeks 11–16 (mid-Dec to Jan)
Run the pre-registered tests at 3-mer level, then 512 5-mers computed directly. Fit nested models against the baselines. Train the ML surrogate, check it on the 512, and extend to 7-mers. **Joint meeting:** agree on the result, positive or null.

### Phase 4, weeks 17–24 (Feb to Mar)
Unchanged: methods, figures, Zenodo DOI, bioRxiv, submission.

---

## 4. Smoke test result (2026-10-01, GFN2-xTB, gas phase, single methyl-capped G:C)

| Quantity | Value |
|---|---|
| H-bond heavy-atom distances O6···N4 / N1···N3 / N2···O2 | 2.68 / 2.83 / 2.79 Å (xTB runs slightly short) |
| G\*:C\* relative energy | **+10.1 kcal/mol** (ZPE-corrected +10.0) |
| Forward / reverse barrier (scan maximum) | **12.5 / 2.5 kcal/mol** |
| Single-proton ion pair G⁻C⁺ | Not a minimum; it relaxes back to G:C |
| WKB κ (310 K), adiabatic / sudden | 1.04 / 2.56 |

**Reading:** the tautomer sits in a well only about 2.5 kcal/mol deep. This is qualitatively what the DFT literature reports (a shallow, high-energy tautomer), and it is exactly why Change 1 tracks well depth and ZPE. The numbers are a pipeline check, **not a result**: one pair, no stacking, no solvent, semiempirical. The raw outputs are in `results/smoke_gc/`.

---

## 5. Risks

| Risk | Mitigation |
|---|---|
| Tautomer is not a minimum once ZPE and environment are included | That is itself a result. Fall back to "barrier/κ as an error-rate proxy" and to the strand-separation and G•T-wobble variants |
| xTB ranks contexts differently from DFT | DFT at stationary points for all 32 (Change 2); switch to DFT or an MLIP if the rank correlation is below about 0.7 |
| No significant correlation at n = 12–16 | Pre-registered, so a null is publishable (as planned). The 5-mer analysis supplies power |
| Al-Hashimi lab publishes a computational follow-up | Post the bioRxiv preprint early; the cross-species replication-error angle stays distinct |
| No HPC access | Ask for an allocation in the October professor emails. xTB work fits a laptop, DFT does not |

## 6. Admin to settle in October
- Author order (equal first authorship) in writing
- Code licence MIT, data CC-BY-4.0 (Zenodo)
- Professor and Surrey emails: attach a one-page summary of this plan and ask specifically for an **HPC allocation** and a **DFT sanity check**
