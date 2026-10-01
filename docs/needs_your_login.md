# Things only Jeremy or Nazia can do (they need your login)

Everything else in the plan has been started or done; see `research_plan.md` §3. These are the steps an automated session cannot do.

| # | What | Where | Time | Why it matters |
|---|---|---|---|---|
| 1 | ~~Create the GitHub repo~~ ✔ done | | | |
| 2 | **Widen this cloud environment's network access** (optional). Add these domains so Claude can read full papers and supplementary files: `pmc.ncbi.nlm.nih.gov`, `www.ncbi.nlm.nih.gov`, `ftp.ncbi.nlm.nih.gov`, `www.nature.com`, `zenodo.org`, `genome.cshlp.org`, `academic.oup.com`, `www.pnas.org`, `cancer.sanger.ac.uk`, `www.biorxiv.org`, `osf.io` | Session title bar → cloud environment → Edit → Network access | 3 min | Today most journal sites were blocked; only GitHub, PyPI and web search worked |
| 3 | **Create an OSF account and paste `docs/osf_preregistration_form.md`** field by field (choose the H2 option first). **Post it before anyone runs the fingerprint-vs-mutation correlation.** Then send Claude the URL | osf.io → Registries → "OSF Preregistration" | 20 min | Makes a null result publishable and blocks the p-hacking objection |
| 4 | Download the **Lujan 2014 yeast mutation tables** (Supplementary) | genome.cshlp.org/content/24/11/1751 → Supplemental Material | 10 min | Yeast replication test (H3). Put the files in `data/raw/yeast/` |
| 5 | Download the **Zou 2021 mutation calls** from Mendeley Data (no EGA request needed for calls) and upload them to the chat or `data/raw/human/` | https://doi.org/10.17632/ymn3ykkmyx | 10 min | Human replication test (H3) |
| 6 | Read Szekely/Al-Hashimi 2026 in full and confirm their method (JSD on signature profiles, no per-site normalisation?) | nature.com/articles/s41467-026-71596-5 | 1 h | So we describe the difference from our approach correctly |
| 7 | Fill in Nazia's surname, then send the three emails in `docs/outreach_emails.md` with `docs/pdf/one_page_summary.pdf` | your email | 20 min | HPC access, a DFT sanity check, an *E. coli* expert |
| 8 | Agree author order: send each other `docs/author_agreement.md` and reply "agreed" | message | 5 min | Avoids problems later |
| 9 | *(Later, February)* Zenodo account linked to GitHub for a DOI; bioRxiv account | zenodo.org, biorxiv.org | 15 min | Phase 4 |

**No longer needed:** ORCA and 3DNA registration. We now use **PySCF** (free, pip-installable) for DFT and **PyMOL** `fnab` for B-DNA building. Register for ORCA only if a reviewer asks for a specific functional PySCF lacks.
