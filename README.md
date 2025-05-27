# chipseq-analysis
chipseq analysis and RNA seq analysis


ChIP-seq Data Analysis Pipeline Overview

This document outlines the typical ChIP-seq data analysis workflow, explaining the purpose, tools, functions, inputs, and outputs for each step.

1. Raw Data Quality Control

Purpose: Assess the quality of raw sequencing reads.

Tool: FastQCCommand: fastqc input.fastq.gz -o fastqc_results/

Input: Raw FASTQ filesOutput: HTML and ZIP reports summarizing base quality, adapter content, duplication levels, etc.

2. Read Trimming

Purpose: Remove low-quality bases and sequencing adapters.

Tool: CutadaptCommand: cutadapt -u 5 -m 30 -o trimmed/input_trimmed.fastq.gz input.fastq.gz

Input: Raw FASTQ filesOutput: Trimmed FASTQ files (reads with poor quality or adapter contamination removed)

3. Pseudoreplicate Generation (Optional)

Purpose: Evaluate reproducibility by splitting trimmed reads into pseudoreplicates.

Tool: seqtkCommand:

seqtk sample -s1 input_trimmed.fastq.gz 100000 | gzip > rep1.fastq.gz
seqtk sample -s2 input_trimmed.fastq.gz 100000 | gzip > rep2.fastq.gz

Input: Trimmed FASTQ fileOutput: Two pseudoreplicate FASTQ files with a random subset of reads

4. Read Mapping

Purpose: Align reads to the reference genome.

Tool: Bowtie2 + SamtoolsCommand:

samtools sort -@ 2 -o mapped/output.bam \
  <(samtools view -@ 2 -F 4 -q 1 -bS \
    <(bowtie2 -p 3 -x reference/genome -U input.fastq.gz 2> log.txt))

Input: FASTQ file (original or pseudoreplicate), reference genome indexOutput: Sorted BAM file and log file with alignment statistics

5. Peak Calling

Purpose: Identify regions with significant read enrichment (binding sites).

Tool: MACS3Command:

macs3 callpeak -t treatment.bam -c control.bam \
  -f BAM -g hs -n output_prefix -B --broad --outdir peak_calling/

Input: Treatment and control BAM filesOutput: Peak files (.broadPeak or .narrowPeak), signal tracks, model info

6. IDR Analysis (Optional for narrow peaks)

Purpose: Assess reproducibility between replicates using Irreproducible Discovery Rate.

Tool: IDRCommand:

idr --samples rep1_peaks.narrowPeak rep2_peaks.narrowPeak \
  --peak-list pooled_peaks.narrowPeak \
  --input-file-type narrowPeak \
  -o idr_output --plot --use-best-multisummit-IDR

Input: Peak files from replicates and pooled peaksOutput: IDR-ranked peaks, plots, log file

7. Peak Annotation

Purpose: Link peaks to genomic features (e.g., promoters, exons).

Tool: ChIPseeker (R package)Function: annotatePeak()

Input: Peak file (e.g., .broadPeak), TxDb annotation objectOutput: Annotated peaks with genomic region information

8. Blacklist Filtering

Purpose: Remove peaks overlapping problematic genomic regions.

Tool: bedtoolsCommand: bedtools intersect -v -a peaks.broadPeak -b blacklist.bed > peaks.filtered.broadPeak

Input: Peak file and ENCODE blacklist BED fileOutput: Filtered peak file

9. Visualization (IGV or JBrowse)

Purpose: Visually assess peak distribution and reproducibility.

Input: BAM, peak (broad/narrowPeak), DEG list (optional), annotation filesOutput: Genome browser tracks for manual inspection

10. Functional Enrichment (Optional)

Purpose: Identify enriched pathways/biological processes in peak-associated genes.

Tool: clusterProfiler or enrichR (R packages)Input: List of genes linked to peaksOutput: KEGG/GO enrichment results

This workflow provides a standard and reproducible approach to ChIP-seq data analysis, allowing for both discovery and validation of chromatin-associated events.

