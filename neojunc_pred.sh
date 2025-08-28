#!/usr/bin/env bash
set -euo pipefail

####################################################################
# Neo-Junction Pipeline - Optimized for STAR 2.7.10a
# Parameters: --alignSJoverhangMin 8 --alignIntronMin 20 --outSJfilterOverhangMin 15 20 20 20
####################################################################

# ------------------ Define and activate venv ------------------

WORKDIR=$(pwd)
source "${WORKDIR}/venv/bin/activate"
PYTHON="${WORKDIR}/venv/bin/python3"

module load bedtools/2.29.2

# ------------------------- Configuration --------------------------

export INPUTDIR="/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/nextflow_rnaseq/output_hartwig/star_salmon/log/"
export OUTDIR="${WORKDIR}/results"
export ANNOT_DIR="${WORKDIR}/annotations"
export SAMPLES_TXT="${WORKDIR}/samples.txt"
export CANON_SET="${ANNOT_DIR}/canonical_junctions.txt"
export GTEX_SET="${ANNOT_DIR}/gtex_junctions.txt"
export PROTEIN_BED="${ANNOT_DIR}/protein_coding.bed"  # Must be defined!

# STAR-aligned filters
export MIN_JUNC_OVERHANG=8        # --alignSJoverhangMin
export MIN_INTRON_SIZE=20         # --alignIntronMin
export MAX_INTRON_SIZE=1000000    # --alignIntronMax

# Neo-junction filters
export MIN_READS_PER_JUNC=10      # Per-sample read support
export MIN_TOTAL_READS=20         # Cohort-wide depth
export MIN_FREQUENCY=0.01         # 1% spliced frequency
export MIN_SAMPLES_PCT=0.10       # 10% samples for "public" junctions
export MAX_GTEX_FREQUENCY=0.01    # 1% GTEx occurrence

# ------------------------- Sample List Generation --------------------------
if [[ ! -s "${SAMPLES_TXT}" ]]; then
  echo "[INFO] Generating samples.txt from *.SJ.out.tab files..."
  find "${INPUTDIR}" -type f -name "*.SJ.out.tab" \
    | awk -F'/' '{fname=$NF; sub(/\.SJ\.out\.tab$/, "", fname); print fname}' \
    | sort -u > "${SAMPLES_TXT}"
  echo "[INFO] Created ${SAMPLES_TXT} with $(wc -l < "${SAMPLES_TXT}") samples."
fi

# ------------------------- Sample Processing ----------------------

process_sample() {
  local sample=$1
  local sj_file="${INPUTDIR}/${sample}.SJ.out.tab"
  local norm_bed="${OUTDIR}/${sample}.juncs.bed"

  echo "[PROC] Processing ${sample}"
  awk -v MIN_OVERHANG=$MIN_JUNC_OVERHANG \
      -v MIN_INTRON=$MIN_INTRON_SIZE \
      -v MAX_INTRON=$MAX_INTRON_SIZE \
      'BEGIN {OFS="\t"} 
       $6 == 0 && $9 >= MIN_OVERHANG && ($3-$2) >= MIN_INTRON && ($3-$2) <= MAX_INTRON && ($7+$8) >= ENVIRON["MIN_READS_PER_JUNC"] {
         strand = ($4 == 1) ? "+" : ($4 == 2) ? "-" : ".";
         print $1, $2-1, $3, "junction_"NR, $7+$8, strand
       }' "${sj_file}" > "${norm_bed}"
}

while read -r sample; do
  process_sample "$sample"
done < "$SAMPLES_TXT"

# ------------------------- Aggregate & Filter ---------------------
echo "[AGG] Aggregating junctions..."

"${PYTHON}" << 'PYCODE'
import pandas as pd
import os

# Load config
outdir = os.environ['OUTDIR']
samples = [s.strip() for s in open(os.environ['SAMPLES_TXT'])]
total_samples = len(samples)
min_samples = int(max(2, float(os.environ['MIN_SAMPLES_PCT']) * total_samples))

# Load all junctions
dfs = []
for sample in samples:
    bed_file = f"{outdir}/{sample}.juncs.bed"
    if os.path.exists(bed_file):
        df = pd.read_csv(bed_file, sep='\t',
                         names=['chr','start','end','name','reads','strand'])
        df['sample'] = sample
        dfs.append(df)

if not dfs:
    raise SystemExit("[ERROR] No per-sample junctions found to aggregate.")

# Concatenate and build 1-based key to match downstream and other sets
all_df = pd.concat(dfs, ignore_index=True)
all_df['key'] = all_df.apply(lambda r: f"{r.chr}:{r.start+1}:{r.end}:{r.strand}", axis=1)

# Collect contributing samples per junction (sorted, unique)
sample_map = (all_df.groupby('key')['sample']
                     .agg(lambda x: ",".join(sorted(set(x))))
                     .to_dict())

# Aggregate counts per junction
agg = (all_df
       .groupby(['chr','start','end','strand','key'])
       .agg(total_reads=('reads','sum'),
            samples_expressed=('sample','nunique'))
       .reset_index())

# Load canonical junctions set
with open(os.environ['CANON_SET']) as f:
    canon_juncs = set(line.strip() for line in f if line.strip())

# Canonical reads for frequency calculation (0 if not present)
canon_reads = {k: v for k, v in
               agg[agg['key'].isin(canon_juncs)][['key','total_reads']].itertuples(index=False)}
agg['canonical_reads'] = agg['key'].map(canon_reads).fillna(0).astype(int)
agg['frequency'] = agg['total_reads'] / (agg['total_reads'] + agg['canonical_reads'].replace(0, 1))

# Save temp file for bedtools intersect (0-based BED; key in col 4)
agg[['chr','start','end','key','total_reads','strand']].to_csv(
    f"{outdir}/temp_junctions.bed", sep='\t', index=False, header=False)

# Protein-coding filter
os.system(f"bedtools intersect -a {outdir}/temp_junctions.bed -b {os.environ['PROTEIN_BED']} -wa -u > {outdir}/protein_filtered.bed")

protein_keys = set()
protein_filtered_path = f"{outdir}/protein_filtered.bed"
if os.path.exists(protein_filtered_path) and os.path.getsize(protein_filtered_path) > 0:
    protein_keys = set(pd.read_csv(protein_filtered_path, sep='\t', header=None)[3])

# GTEx filter (optional file)
gtex_juncs = set()
gtex_path = os.environ['GTEX_SET']
if os.path.exists(gtex_path) and os.path.getsize(gtex_path) > 0:
    with open(gtex_path) as f:
        gtex_juncs = set(line.strip() for line in f if line.strip())

# Apply final filters
neo_juncs = agg[
    (agg['total_reads'] >= int(os.environ['MIN_TOTAL_READS'])) &
    (~agg['key'].isin(canon_juncs)) &           # Exclude canonical
    (agg['key'].isin(protein_keys)) &           # Protein-coding only
    (~agg['key'].isin(gtex_juncs)) &            # Not in GTEx normals
    (agg['frequency'] > float(os.environ['MIN_FREQUENCY'])) &
    (agg['samples_expressed'] >= min_samples)
].copy()

# Attach contributing sample names
neo_juncs['samples'] = neo_juncs['key'].map(sample_map).fillna("")

# Save results (1-based coordinates) with frequency and samples
cols = ['chr','start','end','strand','total_reads','samples_expressed','frequency','samples']
neo_juncs[cols].to_csv(f"{outdir}/neo_junctions.annotated.bed",
                       sep='\t', index=False, header=True)

print(f"[RESULT] Found {len(neo_juncs)} neo-junctions passing filters and annotated with sample names.")
PYCODE
echo "[DONE] Results saved to ${OUTDIR}/neo_junctions.bed"
