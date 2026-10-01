# Literature

✔ = existence and key claim checked on 2026-10-01 · ☐ = still to verify (Nazia)

## Closest competitors
- ✔ Szekely O., …, Al-Hashimi H.M. (2026). *Assessing the contribution of rare DNA states to cancer mutational signatures using sequence-specific conformational fingerprinting.* **Nat Commun** 17:5615. https://www.nature.com/articles/s41467-026-71596-5 · PubMed 41333427
  - ¹⁹F NMR, anionic WC-like G•T⁻ in all 16 triplet contexts. Up to 50-fold variation, driven by the 3′ neighbour. Compared with COSMIC.
  - ✔ Data and code: https://github.com/alhashimilab/DNA-conformational-fingerprinting (Zenodo 10.5281/zenodo.18459447). The 16 G•T⁻ values plus A:T/G:C Hoogsteen and A:T base-opening fingerprints are in `Mutational-Signatures-JSD/JSD-python/Input-files/`.
  - ✔ Method: Jensen–Shannon divergence between the fingerprint and each COSMIC signature, with a random-distribution null. Best G•T⁻ match: **SBS11 (temozolomide)**. ☐ Nazia: confirm from the full text whether MMR signatures were tested and how they scored.
- ✔ Masuda K., Sahakyan A.B. (2024). *Quantum mechanical electronic and geometric parameters for DNA k-mers as features for machine learning.* **Sci Data** 11. https://www.nature.com/articles/s41597-024-03772-5 · data https://github.com/SahakyanLab/DNAkmerQM (CC-BY-4.0) · code https://github.com/SahakyanLab/NucleicAcidsQM
  - All 7-mers in B/A/Z DNA, used for an A>C mutation-rate model.
- ✔ Jago M.J., Green R., Czernuszka M.R., Denisov S., Krašovec R., Knight C.G., Lagator M. (2026). *Extended sequence context shapes mutational bias in Escherichia coli.* **PNAS** 123(23):e2601345123. https://www.pnas.org/doi/10.1073/pnas.2601345123 · data and code https://github.com/Lagator-Group/extended-sequence-context (GPL-3)
  - More than 100k mutations, 32 MA experiments, context effects out to 6 bp and about 1 kb.

## Proton transfer / tunnelling mechanism
- ✔ Slocombe L., Sacchi M., Al-Khalili J. (2022). *An open quantum systems approach to proton tunnelling in DNA.* **Commun Phys** 5. https://www.nature.com/articles/s42005-022-00881-8 · arXiv 2110.00113
- ✔ Slocombe L. et al. (2022). *Proton transfer during DNA strand separation as a source of mutagenic guanine-cytosine tautomers.* **Commun Chem** 5. https://www.nature.com/articles/s42004-022-00760-x
- ✔ Slocombe L., Winokan M., Al-Khalili J., Sacchi M. (2023). *Quantum tunnelling effects in the guanine-thymine wobble misincorporation via tautomerism.* **J Phys Chem Lett.**
- ✔ Soler-Polo D. et al. (2019). *Proton transfer in guanine-cytosine base pairs in B-DNA.* **JCTC** 15(12). https://pubs.acs.org/doi/10.1021/acs.jctc.9b00757 (**counter-evidence**: G\*C\* unstable in cellular conditions)
- ✔ (2026) *Solvation reshapes the classical but not the quantum mechanism of double proton transfer in guanine–cytosine.* **PCCP**, doi 10.1039/d6cp01847e (B3LYP/def2-TZVP 2D surfaces + variational WKB, which is close to our method)
- ✔ (2025) *Unexpected suppression of double-proton tunneling induced by quantum barriers from zero-point energy.* **J Org Chem**, doi 10.1021/acs.joc.5c00827
- ✔ (2024) *Electronic and nuclear quantum effects on proton transfer reactions of G-T mispairs using QM/MM and ML potentials.* **Molecules** 29:2703
- ✔ Gheorghiu A., Coveney P.V., Arabi A.A. (2020). *The influence of base pair tautomerism on single point mutations in aqueous DNA.* **Interface Focus** 10(6). https://discovery.ucl.ac.uk/id/eprint/10113646/ (**counter-evidence**: G\*:C\* half-life is only picoseconds, so the contribution is "negligible")
- ✔ Kimsey I.J. et al. (2018). Dynamic basis for dG•dT misincorporation via tautomerization and ionization. **Nature**, doi 10.1038/nature25487 (G•T tautomer/anion populations, R1ρ)

## Mutation data
- ✔ Hasenauer F.C. et al. (2025). *Genome-wide mapping of spontaneous DNA replication error-hotspots using mismatch repair proteins in rapidly proliferating E. coli.* **NAR** 53(2):gkae1196
- ✔ Foster P.L. et al. (2018). *Determinants of base-pair substitution patterns revealed by whole-genome sequencing of DNA mismatch repair defective E. coli.* **Genetics** (bioRxiv 10.1101/346874)
- ✔ Garushyants S.K., Sane M., Selifanova M.V., Agashe D., Bazykin G.A., Gelfand M.S. (2024). *Mutational signatures in wild type Escherichia coli strains reveal predominance of DNA polymerase errors.* **GBE** 16(4):evae035
- ✔ Lujan S.A. et al. (2014). *Heterogeneous polymerase fidelity and mismatch repair bias genome variation and composition.* **Genome Res** 24:1751–1764
- ✔ Zou X. et al. (2021). *A systematic CRISPR screen defines mutational mechanisms underpinning signatures caused by replication errors and endogenous DNA damage.* **Nat Cancer**
- ✔ Seplyarskiy V. et al. (2023). *A mutation rate model at the basepair resolution identifies the mutagenic effect of polymerase III transcription.* **Nat Genet**. https://github.com/vseplyarskiy/Roulette
