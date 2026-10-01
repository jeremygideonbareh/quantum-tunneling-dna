#!/usr/bin/env bash
# Fetch the public inputs into data/raw/. Re-runnable; skips what is already present.
# Everything here comes from GitHub or PyPI, so it works even on restricted networks.
#   bash scripts/download_data.sh
set -euo pipefail
cd "$(dirname "$0")/../data/raw"

clone() {  # repo dir
  [ -d "$2" ] && { echo "have $2"; return; }
  git clone --depth 1 "https://github.com/$1.git" "$2" && rm -rf "$2/.git"
}

# E. coli: Jago et al. PNAS 2026 compilation (120k MA mutations) + MG1655 genome. GPL-3.
clone Lagator-Group/extended-sequence-context extended-sequence-context

# Al-Hashimi lab conformational fingerprints (Szekely et al. Nat Commun 2026).
clone alhashimilab/DNA-conformational-fingerprinting alhashimi-repo
mkdir -p alhashimi
cp alhashimi-repo/Mutational-Signatures-JSD/JSD-python/Input-files/conf-sig-input-files-for-python-JSDs/*.csv alhashimi/

# DNAkmerQM (Masuda & Sahakyan 2024): QM features for all 6/7-mers. CC-BY-4.0. ~220 MB.
clone SahakyanLab/DNAkmerQM DNAkmerQM

# COSMIC v3.4 signatures + genome context counts, bundled in SigProfiler wheels on PyPI.
if [ ! -s cosmic/COSMIC_v3.4_SBS_GRCh38.txt ]; then
  tmp=$(mktemp -d)
  pip download -q --no-deps -d "$tmp" "SigProfilerAssignment==1.1.5" "SigProfilerMatrixGenerator==1.3.6"
  mkdir -p cosmic
  python - "$tmp" <<'EOF'
import sys, zipfile, glob, os
tmp = sys.argv[1]
want = {
    "Reference_Signatures/GRCh38/COSMIC_v3.4_SBS_GRCh38.txt",
    "context_distributions/context_counts_GRCh38_96.csv",
    "context_distributions/context_counts_GRCh38_1536.csv",
    "context_distributions/context_counts_yeast_96.csv",
    "context_distributions/context_counts_yeast_1536.csv",
}
for whl in glob.glob(os.path.join(tmp, "*.whl")):
    with zipfile.ZipFile(whl) as z:
        for n in z.namelist():
            if any(n.endswith(w) for w in want):
                open(os.path.join("cosmic", os.path.basename(n)), "wb").write(z.read(n))
                print("got ", os.path.basename(n))
EOF
  rm -rf "$tmp"
fi

# Not scriptable from here (journal supplementary files; see docs/needs_your_login.md):
#   Lujan 2014 (yeast), Zou 2021 (human iPSC), Hasenauer 2025 (GEO), Garushyants 2024.
echo "done."
