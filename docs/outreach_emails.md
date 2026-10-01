# Outreach email drafts (October)

Send from your own addresses. Attach `docs/one_page_summary.md` exported as a PDF. Keep each email under 200 words; busy people answer short emails.

---

## 1. Dr Deepa Agashe, NCBS Bengaluru (mutation bias in *E. coli*)

**Why her:** her lab studies *E. coli* mutation spectra and bias (PNAS 2023; PLoS Biol 2025; GBE 2024 with Gelfand on polymerase errors). She is local, and our main dataset is *E. coli*.

> **Subject:** Student project: does DNA sequence context change proton-tunnelling risk? (E. coli MA data)
>
> Dear Dr Agashe,
>
> We are two students (Jeremy Bareh, computation; Nazia [surname], biology) working on a short research project: we compute a quantum-chemical "tunnelling fingerprint" (proton-transfer energetics of G:C and A:T pairs in each sequence context) and test whether it predicts where replication errors occur.
>
> Using the Jago et al. (PNAS 2026) compilation of about 120,000 *E. coli* MA mutations, we find per-site transition rates in MMR-deficient lines vary up to 54-fold across trinucleotide contexts. The pattern also correlates with human MMR signatures (ρ ≈ 0.8 for SBS21). Your work on mutation bias and polymerase errors (Garushyants et al., GBE 2024) is very close to this question.
>
> Could we have 20 minutes of your time for feedback on our analysis plan (one page attached), especially on confounders in the MA data? We will pre-register our hypotheses before testing them.
>
> Thank you,
> Jeremy Bareh and Nazia [surname]

---

## 2. Dr Marco Sacchi / Dr Louie Slocombe, University of Surrey (proton tunnelling in DNA)

**Why them:** they wrote the open-quantum-systems tunnelling papers we build on. Ask for **method validation**, not supervision.

> **Subject:** Sequence-context dependence of G:C / A:T proton transfer: a quick methods question
>
> Dear Dr Sacchi and Dr Slocombe,
>
> We are students running a systematic test of proton-transfer tautomerism against mutation data. We compute double-proton-transfer profiles for the central pair of capped 3-bp B-DNA clusters in all 32 trinucleotide contexts (GFN2-xTB scans, DFT single-point checks, WKB tunnelling in adiabatic and sudden limits). We then test them against replication-error spectra from MMR-deficient *E. coli*, yeast and human cells, with pre-registered hypotheses.
>
> Your 2022 papers showed that tunnelling and strand separation matter. Two quick questions:
> (1) Is a 1D WKB estimate on the scanned profile a defensible proxy for your open-quantum-systems rates, if we only need the **ranking** across contexts?
> (2) Would you expect neighbouring pairs to shift the barrier enough to matter (in our gas-phase xTB clusters the G:C barrier spans 11.7–13.5 kcal/mol across the 16 contexts, and A\*:T\* is never a minimum)?
>
> One-page summary attached. Thank you for any pointers.
> Jeremy Bareh and Nazia [surname]

---

## 3. Local computational-chemistry professor (Bengaluru), for an HPC allocation and a DFT check

**Who:** pick one faculty member in IISc's Solid State & Structural Chemistry Unit or JNCASR's Theoretical Sciences Unit whose page lists DFT or biomolecular simulation. Nazia: shortlist two and check their recent papers.

> **Subject:** Request: brief DFT advice and possible compute access for a student project on DNA proton transfer
>
> Dear Prof. [Name],
>
> We are students computing proton-transfer energetics of DNA base pairs in all sequence contexts, to test whether they explain context-dependent replication errors (one-page plan attached). Semiempirical (xTB) scans for 32 contexts already run on a laptop, and DFT single points (B3LYP-D3/def2-SVP, PySCF) work for single pairs. To validate at def2-TZVP and extend to 512 five-base contexts, we need a modest CPU allocation (about 5,000 core-hours).
>
> Could we meet for 15 minutes to ask whether our DFT protocol is sound, and whether a small allocation or a student account on your group's cluster would be possible? We would acknowledge your help, or offer co-authorship if you contribute substantively.
>
> Thank you,
> Jeremy Bareh and Nazia [surname]
