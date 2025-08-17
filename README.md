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



