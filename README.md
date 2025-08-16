# Neojunction_prediction

Neojunction prediction from 307 transcriptomic dataset acquired from Hartwig (Breast sample)

Found 32 neo-junctions that passed all the filtering criteria

The neo-junctions identified are:

✅ Non-canonical (novel splice sites)
✅ Present in ≥10% of samples (public junctions)
✅ Have sufficient read support (≥20 total reads)
✅ Located in protein-coding regions
✅ Not found in GTEx normal samples
✅ Have >1% splicing frequency

**Goal** (Summary)

- Start with STAR sj.out.tab files.
- Remove junctions annotated in GENCODE(GRch38).
- Keep only junctions overlapping non-mitochondrial, protein-coding genes.
- Remove junctions with <10 spliced reads (per sample) or <20 total reads (cohort).
- Compute spliced frequency:
- Frequency = total target spliced reads / (target + canonical spliced reads).
- Retain junctions with frequency >1%.
- Define “public” junctions:
- Expressed in ≥10% of cohort, with above criteria.
- Remove junctions found in >1% of GTEx normal samples.
- Output: high-confidence, cancer-specific, non-annotated junctions.

**STAR aligner** 

STAR produces SJ.out.tab: Contains every splice junction STAR detected, with per-junction read counts (unique + multi-mapping reads) and annotation flags.

- srun nextflow run nf-core/rnaseq \
 
  --aligner star_salmon \
  --skip_quantification \
  --skip_qc \
   --star_align_args "--outSAMtype BAM SortedByCoordinate \
                     --quantMode TranscriptomeSAM \
                     --outSAMunmapped Within \
                     --outSJfilterOverhangMin 15 20 20 20 \
                     --alignSJoverhangMin 8 \
                     --alignIntronMin 20 \
                     --alignIntronMax 1000000" \
  -profile singularity \
  -c /data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/config_files/icr_alma.config \
  -r 3.18.0 \
  -resume

--outSJfilterOverhangMin 15 20 20 20
Controls the minimum overhang (length of mapped sequence flanking a splice junction) for different splice junction types to be reported. This increases stringency on junction detection:
  -  15 for canonical GT/AG junctions
  -  20 for other types (e.g., GC/AG, AT/AC, non-canonical)
This helps reduce false positives and ensures only well-supported junctions are output.

--alignSJoverhangMin 8
Sets the minimum overhang length required on each side of a splice junction to consider it a valid splice junction during alignment. A value of 8 bases means STAR must see at least 8 matching bases flanking the junction.

--alignIntronMin 20
Specifies the minimum allowed intron length (20 bases) for splice junctions. Very short introns are usually sequencing or alignment artifacts, so this filters those out.

--alignIntronMax 1000000
Sets the maximum allowed intron length (1,000,000 bases). This limits extremely long introns that could be biologically implausible or mapping artifacts.


####################################################################
# Neo-Junction Pipeline - Optimized for STAR 2.7.10a
# Parameters: --alignSJoverhangMin 8 --alignIntronMin 20 --outSJfilterOverhangMin 15 20 20 20
####################################################################

#!/usr/bin/env bash
set -euo pipefail

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
import pandas as pd, os

outdir       = os.environ['OUTDIR']
samples_file = os.environ['SAMPLES_TXT']
min_samples  = int(max(2, float(os.environ['MIN_SAMPLES_PCT'])*len(open(samples_file).read().splitlines())))
               
# 1. Load all junctions with sample labels
dfs = []
for sample in open(samples_file).read().splitlines():
    bed = f"{outdir}/{sample}.juncs.bed"
    if os.path.exists(bed):
        df = pd.read_csv(bed, sep='\t',
                         names=['chr','start','end','name','reads','strand'])
        df['sample'] = sample
        # Create key identical to later agg key
        df['key'] = df.apply(lambda r: f"{r.chr}:{r.start+1}:{r.end}:{r.strand}", axis=1)
        dfs.append(df)
all_df = pd.concat(dfs, ignore_index=True)

# 2. Build sample list per junction
sample_map = all_df.groupby('key')['sample'].unique().to_dict()

# 3. Aggregate reads and count samples
agg = (all_df
       .groupby(['chr','start','end','strand','key'])
       .agg(total_reads=('reads','sum'),
            samples_expressed=('sample','nunique'))
       .reset_index())

# 4. Load canonical & GTEx sets
canon = set(open(os.environ['CANON_SET']).read().splitlines())
gtex  = set(open(os.environ['GTEX_SET']).read().splitlines()) if os.path.exists(os.environ['GTEX_SET']) else set()

# 5. Protein filter via bedtools intersect file
prot_bed = f"{outdir}/protein_filtered.bed"
protein_keys = set(pd.read_csv(prot_bed, sep='\t', header=None)[3]) if os.path.exists(prot_bed) else set()

# 6. Compute frequency
#    (use 1 when no canonical_reads to avoid div0)
canon_reads = {row.key: row.total_reads for _,row in agg[agg.key.isin(canon)].iterrows()}
agg['canonical_reads'] = agg['key'].map(canon_reads).fillna(0).astype(int)
agg['frequency'] = agg['total_reads'] / (agg['total_reads'] + agg['canonical_reads'].replace(0,1))

# 7. Apply final filters
neo = agg[
    (agg.total_reads       >= int(os.environ['MIN_TOTAL_READS'])) &
    (agg.key.isin(protein_keys))              &  # protein only
    (~agg.key.isin(canon))                    &  # exclude canonical
    (~agg.key.isin(gtex))                     &  # exclude GTEx
    (agg.frequency        > float(os.environ['MIN_FREQUENCY'])) &
    (agg.samples_expressed>= min_samples)
].copy()

# 8. Annotate samples list
neo['samples'] = neo['key'].map(sample_map).apply(lambda arr: ",".join(sorted(arr)))

# 9. Write annotated output
cols = ['chr','start','end','strand','total_reads','samples_expressed','frequency','samples']
neo.to_csv(f"{outdir}/neo_junctions.annotated.bed",
           sep='\t', index=False, header=True)

print(f"[RESULT] Found {len(neo)} neo-junctions. Detailed output in neo_junctions.annotated.bed")
PYCODE

echo "[DONE] Annotated results saved to ${OUTDIR}/neo_junctions.annotated.bed"

