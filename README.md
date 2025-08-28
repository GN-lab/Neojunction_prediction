# Neojunction_prediction

Neojunction prediction from 307 transcriptomic dataset acquired from Hartwig (Breast sample)

Found 32 neo-junctions observed in 282 samples that passed all the filtering criteria

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
# To generate Annotation files
####################################################################

#!/usr/bin/env bash

# For canonical_junctions.txt
zcat ../../../Ref/genome.gtf.gz | awk -F'\t' '
$3 == "exon" && $9 ~ "gene_biotype \"protein_coding\"" {
    # Extract gene_name
    gene_name = "";
    split($9, attrs, ";");
    for (i in attrs) {
        if (attrs[i] ~ /gene_name/) {
            split(attrs[i], gn, "\"");
            gene_name = gn[2];
            break;
        }
    }
    # Print in BED format
    if (gene_name != "") {
        print $1 "\t" $4-1 "\t" $5 "\t" gene_name "\t0\t" $7;
    }
}' | sort -k1,1 -k2,2n | uniq > canonical_junctions.txt

# For protein_coding.bed
zcat ../../../Ref/genome.gtf.gz | awk -F'\t' '
$3 == "gene" && $9 ~ "gene_biotype \"protein_coding\"" {
    # Extract gene_name
    gene_name = "";
    split($9, attrs, ";");
    for (i in attrs) {
        if (attrs[i] ~ /gene_name/) {
            split(attrs[i], gn, "\"");
            gene_name = gn[2];
            break;
        }
    }
    # Print in BED format
    if (gene_name != "") {
        print $1 "\t" $4-1 "\t" $5 "\t" gene_name "\t0\t" $7;
    }
}' | sort -k1,1 -k2,2n | uniq > protein_coding.bed

# For Gtex junction file
zcat GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz \
  | tail -n +4 \
  | cut -f1 \
  | awk -F'_' '{
      chrom=$1; sub(/^chr/, "", chrom);
      strand = ($4 == "1") ? "+" : "-";
      print chrom ":" $2 ":" $3 ":" strand;
    }' \
  > gtex_junctions_fixed.txt


####################################################################
# Results (plots)
####################################################################


<img width="4800" height="8700" alt="neo_junctions_heatmap" src="https://github.com/user-attachments/assets/9d088474-66aa-4beb-9276-5c32de8e20f8" />


<img width="2400" height="1800" alt="neo_junctions_distribution" src="https://github.com/user-attachments/assets/8ff0a045-0b5e-49ef-8281-d9313b33a3c2" />


<img width="3000" height="3000" alt="fixed_sample_chr_circos" src="https://github.com/user-attachments/assets/10efb18e-4886-4763-8ba2-d255a08d0e1f" />


<img width="4800" height="1800" alt="neo_junctions_summary" src="https://github.com/user-attachments/assets/10717431-1364-4e40-b764-7db7e17f275e" />


<img width="6000" height="1800" alt="neo_junctions_per_sample" src="https://github.com/user-attachments/assets/d7a009bc-b554-41d4-a187-67e27c62e7ce" />
