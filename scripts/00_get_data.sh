#!/usr/bin/env bash
# 00_get_data.sh: download the demo dataset and the SILVA reference
# Demo data: Mothur MiSeq SOP (Schloss lab) mouse gut 16S V4, 2x250 Illumina MiSeq.
# 19 fecal samples from one mouse over 150 days post-weaning, plus a 20-strain mock community.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/raw ref

if [ ! -f data/raw/F3D0_S188_L001_R1_001.fastq ]; then
  curl -L -o data/MiSeqSOPData.zip https://mothur.s3.amazonaws.com/wiki/miseqsopdata.zip
  unzip -q -o data/MiSeqSOPData.zip -d data/tmp
  find data/tmp -name "*.fastq" -exec mv {} data/raw/ \;
  find data/tmp -name "HMP_MOCK.v35.fasta" -exec mv {} data/ \;
  rm -rf data/tmp data/MiSeqSOPData.zip data/__MACOSX
fi

# SILVA 138.2 formatted for DADA2 (Callahan, Zenodo record 14169026)
[ -f ref/silva_nr99_v138.2_toGenus_trainset.fa.gz ] || \
  curl -L -o ref/silva_nr99_v138.2_toGenus_trainset.fa.gz \
  "https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1"
[ -f ref/silva_v138.2_assignSpecies.fa.gz ] || \
  curl -L -o ref/silva_v138.2_assignSpecies.fa.gz \
  "https://zenodo.org/records/14169026/files/silva_v138.2_assignSpecies.fa.gz?download=1"

# Sample metadata: days 0-9 post-weaning = Early, days 141-150 = Late
python3 - <<'EOF'
import os, re
rows = ["SampleID\tGroup\tDay"]
for f in sorted(os.listdir("data/raw")):
    if not f.endswith("_R1_001.fastq"): continue
    s = f.split("_")[0]
    if s == "Mock":
        rows.append(f"{s}\tMock\tNA"); continue
    day = int(re.match(r"F3D(\d+)", s).group(1))
    rows.append(f"{s}\t{'Early' if day <= 9 else 'Late'}\t{day}")
open("data/metadata.tsv", "w").write("\n".join(rows) + "\n")
print(f"metadata.tsv: {len(rows)-1} samples")
EOF
echo "Data ready."
