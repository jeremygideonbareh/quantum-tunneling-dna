#!/usr/bin/env bash
# Fetch the public inputs into data/raw/. Re-runnable; skips files already present.
# Usage: bash scripts/download_data.sh [--with-roulette]
set -euo pipefail
cd "$(dirname "$0")/../data/raw"

fetch() {  # url outfile
  [ -s "$2" ] && { echo "have $2"; return; }
  echo "get  $2"; curl -fL --retry 3 -o "$2.part" "$1" && mv "$2.part" "$2"
}

# COSMIC SBS v3.4 reference signatures (GRCh38). If the URL has moved, download
# by hand from https://cancer.sanger.ac.uk/signatures/downloads/
fetch https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.4_SBS_GRCh38.txt \
      COSMIC_v3.4_SBS_GRCh38.txt

# DNAkmerQM (Masuda & Sahakyan 2024): QM features for all 7-mers, B/A/Z DNA. CC-BY-4.0.
[ -d DNAkmerQM ] || git clone --depth 1 https://github.com/SahakyanLab/DNAkmerQM.git

# Reference genomes for opportunity counts
fetch "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz" \
      ecoli_MG1655.fna.gz
fetch "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz" \
      yeast_S288C_R64.fna.gz

# Roulette germline rates (Seplyarskiy et al. 2023): per-chromosome VCFs, very large.
# Check sizes on http://genetics.bwh.harvard.edu/downloads/Vova/Roulette/ before pulling.
if [[ "${1:-}" == "--with-roulette" ]]; then
  mkdir -p roulette
  echo "Download the per-chromosome files listed at the URL above into data/raw/roulette/"
fi

# Not scriptable (supplementary tables; Nazia confirms exact files in Phase 1):
#   Foster et al. 2018 Genetics      - E. coli MMR-defective MA mutations
#   Lujan et al. 2014 Genome Res     - yeast pol/msh2 mutations
#   Zou et al. 2021 Nature Cancer    - human iPSC repair-knockout mutations
#   Hasenauer et al. 2025 NAR        - E. coli MutS/MutL ChIP-seq hotspots (GEO)
#   Lagator lab, PNAS 2026           - 100k+ E. coli MA mutations, 32 experiments
echo "done. See docs/data_sources.md for the manual downloads."
