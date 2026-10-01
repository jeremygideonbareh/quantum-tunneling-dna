# Things only Jeremy or Nazia can do (they need your login)

Everything else in the plan has been started or done; see `research_plan.md` §3. These are the steps an automated session cannot do.

| # | What | Where | Time | Why it matters |
|---|---|---|---|---|
| 1 | **Create an empty private repo `quantum-tunneling-dna`** (no README), then make sure the Claude GitHub App can access it | github.com/new, then https://claude.ai/connect-github | 2 min | Then Claude can push all the code, data tables and docs |
| 2 | **Widen this cloud environment's network access** (optional). Add these domains so Claude can read full papers and supplementary files: `pmc.ncbi.nlm.nih.gov`, `www.ncbi.nlm.nih.gov`, `ftp.ncbi.nlm.nih.gov`, `www.nature.com`, `zenodo.org`, `genome.cshlp.org`, `academic.oup.com`, `www.pnas.org`, `cancer.sanger.ac.uk`, `www.biorxiv.org`, `osf.io` | Session title bar → cloud environment → Edit → Network access | 3 min | Today most journal sites were blocked; only GitHub, PyPI and web search worked |
| 3 | **Create an OSF account** and start a pre-registration from `docs/preregistration.md`. **Post it before anyone runs the fingerprint-vs-mutation correlation** | osf.io → Registries → "OSF Preregistration" | 30 min | Makes a null result publishable and blocks the p-hacking objection |
| 4 | Download the **Lujan 2014 yeast mutation tables** (Supplementary) | genome.cshlp.org/content/24/11/1751 → Supplemental Material | 10 min | Yeast replication test (H3). Put the files in `data/raw/yeast/` |
| 5 | Find the **Zou 2021 human iPSC knockout data** link (data availability statement) | nature.com/articles/s43018-021-00200-0 | 10 min | Human replication test (H3); the raw data may need an EGA access request |
| 6 | Read Szekely/Al-Hashimi 2026 in full and confirm their method (JSD on signature profiles, no per-site normalisation?) | nature.com/articles/s41467-026-71596-5 | 1 h | So we describe the difference from our approach correctly |
| 7 | Send the three emails in `docs/outreach_emails.md` | your email | 30 min | HPC access, a DFT sanity check, an *E. coli* expert |
| 8 | Agree author order in writing (a message to each other is enough) | — | 5 min | Avoids problems later |
| 9 | *(Later, February)* Zenodo account linked to GitHub for a DOI; bioRxiv account | zenodo.org, biorxiv.org | 15 min | Phase 4 |

**No longer needed:** ORCA and 3DNA registration. We now use **PySCF** (free, pip-installable) for DFT and **PyMOL** `fnab` for B-DNA building. Register for ORCA only if a reviewer asks for a specific functional PySCF lacks.
