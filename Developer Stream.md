# Developer Stream

the process

1. environment setting
    1. miniforge3
    2. initialize
        
        ```python
        source ~/miniforge3/bin/activate
        ```
        
    3. run installer
        
        ```python
        conda install -n base -c conda-forge mamba
        ```
        
        ```python
        mamba create -n chipseq_env python=3.10
        ```
        
        ```python
        eval "$(mamba shell hook --shell zsh)”
        ```
        
        ```python
        mamba activate chipseq_env
        ```
        
2. install pachages
    
    ```python
    mamba install bioconda::fastqc
    ```
    
    - `cutadapt`
    - `bowtie2`
    - `phantompeakqualtools`**(not availible locally, and will be used in step 6. for instead, run it in virtual machines)**
    - `macs3`
    - `bedtools`
3. read quality check
    
    ```python
    mkdir -p fastqc_results/raw
    ```
    
    ```python
    fastqc -t 4 -o fastqc_results/raw \
    data/control.fastq.gz \
    data/test1.fastq.gz \
    data/test2.fastq.gz \
    data/test3.fastq.gz
    ```
    
4. **Remove adapters**
    
    according to the fastqc report, no significant adapter contamination has been found.
    
    But for each threads, [Per base sequence content](https://www.notion.so/Developer-Stream-203a3e4eb09880d1b2c6e42eff90dab7?pvs=21) no work out
    
    and for test 1 and 3, [Per sequence GC content](https://www.notion.so/Developer-Stream-203a3e4eb09880d1b2c6e42eff90dab7?pvs=21) show attention sign -**It is not higher than test2 GC performance -  doesn’t matter**
    
    and for test 2, [Sequence Duplication Levels](https://www.notion.so/Developer-Stream-203a3e4eb09880d1b2c6e42eff90dab7?pvs=21) show attention sign
    

trim

```python
mkdir -p trimmed
```

```python
cutadapt -u 5 -m 0 -o trimmed/control_trimmed.fastq.gz data/control.fastq.g
cutadapt -u 5 -m 30 -o trimmed/test1_trimmed.fastq.gz data/test1.fastq.gz
cutadapt -u 5 -m 30 -o trimmed/test2_trimmed.fastq.gz data/test2.fastq.gz
cutadapt -u 5 -m 30 -o trimmed/test3_trimmed.fastq.gz data/test3.fastq.gz
```

and quality check

```python
fastqc trimmed/control_trimmed.fastq.gz \
trimmed/test1_trimmed.fastq.gz \
trimmed/test2_trimmed.fastq.gz \
trimmed/test3_trimmed.fastq.gz \
-o workshop2_quality_control_reports -t 4
```

1. create **pseudoreplicates**

```python
mamba install -c bioconda seqtk
mkdir -p pseudoreps

# test1
seqtk sample -s100 trimmed/test1_trimmed.fastq.gz 0.5 | gzip > pseudoreps/test1_rep1.fastq.gz
seqtk sample -s200 trimmed/test1_trimmed.fastq.gz 0.5 | gzip > pseudoreps/test1_rep2.fastq.gz

# test2
seqtk sample -s100 trimmed/test2_trimmed.fastq.gz 0.5 | gzip > pseudoreps/test2_rep1.fastq.gz
seqtk sample -s200 trimmed/test2_trimmed.fastq.gz 0.5 | gzip > pseudoreps/test2_rep2.fastq.gz

# test3
seqtk sample -s100 trimmed/test3_trimmed.fastq.gz 0.5 | gzip > pseudoreps/test3_rep1.fastq.gz
seqtk sample -s200 trimmed/test3_trimmed.fastq.gz 0.5 | gzip > pseudoreps/test3_rep2.fastq.gz

```

1. read mapping
    
    ```python
    mkdir -p mapped/bowtie2
    ```
    
    Notice! here is breaking down chain command!
    
    bowtie2- **the reference is hg19**
    
    ```python
    # control
    samtools sort -@ 2 -o mapped/bowtie2/control.bam \
      <(samtools view -@ 2 -F 4 -q 1 -bS \
        <(bowtie2 -p 3 -x reference/chromosomes/hg19 -U trimmed/control_trimmed.fastq.gz \
          2>mapped/bowtie2/control.log))
          
    # ---------- test1 ----------
    # test1
    samtools sort -@ 2 -o mapped/bowtie2/test1.bam \
      <(samtools view -@ 2 -F 4 -q 1 -bS \
        <(bowtie2 -p 3 -x reference/chromosomes/hg19 -U trimmed/test1_trimmed.fastq.gz \
          2>mapped/bowtie2/test1.log))
    
    # pseudorep1
    samtools sort -@ 2 -o mapped/bowtie2/test1_rep1.bam \
      <(samtools view -@ 2 -F 4 -q 1 -bS \
        <(bowtie2 -p 3 -x reference/chromosomes/hg19 -U pseudoreps/test1_rep1.fastq.gz \
          2>mapped/bowtie2/test1_rep1.log))
    
    # pseudorep2
    samtools sort -@ 2 -o mapped/bowtie2/test1_rep2.bam \
      <(samtools view -@ 2 -F 4 -q 1 -bS \
        <(bowtie2 -p 3 -x reference/chromosomes/hg19 -U pseudoreps/test1_rep2.fastq.gz \
          2>mapped/bowtie2/test1_rep2.log))
    
    # ---------- test2 ----------
    samtools sort -@ 2 -o mapped/bowtie2/test2.bam \
      <(samtools view -@ 2 -F 4 -q 1 -bS \
        <(bowtie2 -p 3 -x reference/chromosomes/hg19 -U trimmed/test2_trimmed.fastq.gz \
          2>mapped/bowtie2/test2.log))
    
    samtools sort -@ 2 -o mapped/bowtie2/test2_rep1.bam \
      <(samtools view -@ 2 -F 4 -q 1 -bS \
        <(bowtie2 -p 3 -x reference/chromosomes/hg19 -U pseudoreps/test2_rep1.fastq.gz \
          2>mapped/bowtie2/test2_rep1.log))
    
    samtools sort -@ 2 -o mapped/bowtie2/test2_rep2.bam \
      <(samtools view -@ 2 -F 4 -q 1 -bS \
        <(bowtie2 -p 3 -x reference/chromosomes/hg19 -U pseudoreps/test2_rep2.fastq.gz \
          2>mapped/bowtie2/test2_rep2.log))
    
    # ---------- test3 ----------
    samtools sort -@ 2 -o mapped/bowtie2/test3.bam \
      <(samtools view -@ 2 -F 4 -q 1 -bS \
        <(bowtie2 -p 3 -x reference/chromosomes/hg19 -U trimmed/test3_trimmed.fastq.gz \
          2>mapped/bowtie2/test3.log))
    
    samtools sort -@ 2 -o mapped/bowtie2/test3_rep1.bam \
      <(samtools view -@ 2 -F 4 -q 1 -bS \
        <(bowtie2 -p 3 -x reference/chromosomes/hg19 -U pseudoreps/test3_rep1.fastq.gz \
          2>mapped/bowtie2/test3_rep1.log))
    
    samtools sort -@ 2 -o mapped/bowtie2/test3_rep2.bam \
      <(samtools view -@ 2 -F 4 -q 1 -bS \
        <(bowtie2 -p 3 -x reference/chromosomes/hg19 -U pseudoreps/test3_rep2.fastq.gz \
          2>mapped/bowtie2/test3_rep2.log))
    ```
    
    check the amount of reads
    
    ```python
    samtools flagstat mapped/bowtie2/control.bam 
    samtools flagstat mapped/bowtie2/test1.bam 
    samtools flagstat mapped/bowtie2/test1_rep1.bam 
    samtools flagstat mapped/bowtie2/test1_rep2.bam 
    samtools flagstat mapped/bowtie2/test2.bam 
    samtools flagstat mapped/bowtie2/test2_rep1.bam 
    samtools flagstat mapped/bowtie2/test2_rep2.bam 
    samtools flagstat mapped/bowtie2/test3.bam
    samtools flagstat mapped/bowtie2/test3_rep1.bam 
    samtools flagstat mapped/bowtie2/test3_rep2.bam 
    ```
    
2.  **Peak quality check（not acceptable locally）**
    
    method 1 - try to use other dependency
    
    for the r-spp dependence, try to install manually
    
    ```python
    mamba install -n phantom_env r-spp -c bioconda
    ```
    
    method 2 - try to create another env only for phantompeakqualtools- failed
    
    phantpmpealqualtools only have packages for `osx-64` and `linux-64`, no package for M1
    
    ```python
    mamba create -n phantom_env r-base=4.1 phantompeakqualtools -c bioconda -c conda-forge
    ```
    
    method 3-  try to use the different R package - channel error
    
    ```python
    mamba create -n phantom_env r-base=4.0 phantompeakqualtools -c bioconda -c conda-forge 
    ```
    

1. **Peak calling**
    
    ```python
    # 创建所有输出目录
    mkdir -p peak_calling/test1 peak_calling/test1_rep1 peak_calling/test1_rep2
    mkdir -p peak_calling/test2 peak_calling/test2_rep1 peak_calling/test2_rep2
    mkdir -p peak_calling/test3 peak_calling/test3_rep1 peak_calling/test3_rep2
    
    # ---------- test1 ----------
    macs3 callpeak \
      -t mapped/bowtie2/test1.bam \
      -c mapped/bowtie2/control.bam \
      -f BAM -g hs -n test1 -B \
      --outdir peak_calling/test1
    
    # test1 rep1
    macs3 callpeak \
      -t mapped/bowtie2/test1_rep1.bam \
      -c mapped/bowtie2/control.bam \
      -f BAM -g hs -n test1_rep1 -B \
      --outdir peak_calling/test1_rep1
    
    # test1 rep2
    macs3 callpeak \
      -t mapped/bowtie2/test1_rep2.bam \
      -c mapped/bowtie2/control.bam \
      -f BAM -g hs -n test1_rep2 -B \
      --outdir peak_calling/test1_rep2
    
    # ---------- test2 ----------
    macs3 callpeak \
      -t mapped/bowtie2/test2.bam \
      -c mapped/bowtie2/control.bam \
      -f BAM -g hs -n test2 -B \
      --outdir peak_calling/test2
    
    # test2 rep1
    macs3 callpeak \
      -t mapped/bowtie2/test2_rep1.bam \
      -c mapped/bowtie2/control.bam \
      -f BAM -g hs -n test2_rep1 -B \
      --outdir peak_calling/test2_rep1
    
    # test2 rep2
    macs3 callpeak \
      -t mapped/bowtie2/test2_rep2.bam \
      -c mapped/bowtie2/control.bam \
      -f BAM -g hs -n test2_rep2 -B \
      --outdir peak_calling/test2_rep2
    
    # ---------- test3 ----------
    macs3 callpeak \
      -t mapped/bowtie2/test3.bam \
      -c mapped/bowtie2/control.bam \
      -f BAM -g hs -n test3 -B \
      --outdir peak_calling/test3
    
    # test3 rep1
    macs3 callpeak \
      -t mapped/bowtie2/test3_rep1.bam \
      -c mapped/bowtie2/control.bam \
      -f BAM -g hs -n test3_rep1 -B \
      --outdir peak_calling/test3_rep1
    
    # test3 rep2
    macs3 callpeak \
      -t mapped/bowtie2/test3_rep2.bam \
      -c mapped/bowtie2/control.bam \
      -f BAM -g hs -n test3_rep2 -B \
      --outdir peak_calling/test3_rep2
    
    ```
    
    ```python
    #board peaks
    # test1 
    macs3 callpeak \
      -t mapped/bowtie2/test1.bam \
      -c mapped/bowtie2/control.bam \
      -f BAM -g hs -n test1_broad -B \
      --broad --broad-cutoff 0.1 \
      --nomodel --extsize 147 \
      --outdir peak_calling/test1_broad
    
    # test1_rep1
    macs3 callpeak \
      -t mapped/bowtie2/test1_rep1.bam \
      -c mapped/bowtie2/control.bam \
      -f BAM -g hs -n test1_broad_rep1 -B \
      --broad --broad-cutoff 0.1 \
      --nomodel --extsize 147 \
      --outdir peak_calling/test1_broad_rep1
    
    # test1_rep2
    macs3 callpeak \
      -t mapped/bowtie2/test1_rep2.bam \
      -c mapped/bowtie2/control.bam \
      -f BAM -g hs -n test1_broad_rep2 -B \
      --broad --broad-cutoff 0.1 \
      --nomodel --extsize 147 \
      --outdir peak_calling/test1_broad_rep2
    
    # test2-------
    # test2
    macs3 callpeak \
      -t mapped/bowtie2/test2.bam \
      -c mapped/bowtie2/control.bam \
      -f BAM -g hs -n test2_broad -B \
      --broad --broad-cutoff 0.1 \
      --nomodel --extsize 147 \
      --outdir peak_calling/test2_broad
    
    # test2_rep1
    macs3 callpeak \
      -t mapped/bowtie2/test2_rep1.bam \
      -c mapped/bowtie2/control.bam \
      -f BAM -g hs -n test2_broad_rep1 -B \
      --broad --broad-cutoff 0.1 \
      --nomodel --extsize 147 \
      --outdir peak_calling/test2_broad_rep1
    
    # test2_rep2
    macs3 callpeak \
      -t mapped/bowtie2/test2_rep2.bam \
      -c mapped/bowtie2/control.bam \
      -f BAM -g hs -n test2_broad_rep2 -B \
      --broad --broad-cutoff 0.1 \
      --nomodel --extsize 147 \
      --outdir peak_calling/test2_broad_rep2
      
    # test3-------
    # test3
    macs3 callpeak \
      -t mapped/bowtie2/test3.bam \
      -c mapped/bowtie2/control.bam \
      -f BAM -g hs -n test3_broad -B \
      --broad --broad-cutoff 0.1 \
      --nomodel --extsize 147 \
      --outdir peak_calling/test3_broad
    
    # test3_rep1
    macs3 callpeak \
      -t mapped/bowtie2/test3_rep1.bam \
      -c mapped/bowtie2/control.bam \
      -f BAM -g hs -n test3_broad_rep1 -B \
      --broad --broad-cutoff 0.1 \
      --nomodel --extsize 147 \
      --outdir peak_calling/test3_broad_rep1
    
    # test3_rep2
    macs3 callpeak \
      -t mapped/bowtie2/test3_rep2.bam \
      -c mapped/bowtie2/control.bam \
      -f BAM -g hs -n test3_broad_rep2 -B \
      --broad --broad-cutoff 0.1 \
      --nomodel --extsize 147 \
      --outdir peak_calling/test3_broad_rep2
    ```
    

1. instal idr(not use)
    
    ```python
    mamba create -n idr_run
    mamba activate idr_run 
    CONDA_SUBDIR=osx-64 conda install -c bioconda -c conda-forge idr
    ```
    
2. run idr(not use)
    
    ```python
    # test1 IDR
    idr \
      --samples peak_calling/test1_rep1/test1_rep1_peaks.narrowPeak \
               peak_calling/test1_rep2/test1_rep2_peaks.narrowPeak \
      --peak-list peak_calling/test1/test1_peaks.narrowPeak \
      --input-file-type narrowPeak \
      -o peak_calling/test1.idr \
      --plot \
      --use-best-multisummit-IDR \
      --log-output-file peak_calling/test1.idr.log
    
    # test2 IDR
    idr \
      --samples peak_calling/test2_rep1/test2_rep1_peaks.narrowPeak \
               peak_calling/test2_rep2/test2_rep2_peaks.narrowPeak \
      --peak-list peak_calling/test2/test2_peaks.narrowPeak \
      --input-file-type narrowPeak \
      -o peak_calling/test2.idr \
      --plot \
      --use-best-multisummit-IDR \
      --log-output-file peak_calling/test2.idr.log
    
    # test3 IDR
    idr \
      --samples peak_calling/test3_rep1/test3_rep1_peaks.narrowPeak \
               peak_calling/test3_rep2/test3_rep2_peaks.narrowPeak \
      --peak-list peak_calling/test3/test3_peaks.narrowPeak \
      --input-file-type narrowPeak \
      -o peak_calling/test3.idr \
      --plot \
      --use-best-multisummit-IDR \
      --log-output-file peak_calling/test3.idr.log
    
    # board peaks--------
    
    # test1 IDR
    idr \
      --samples peak_calling/test1_broad_rep1/test1_broad_rep1_peaks.broadPeak \
               peak_calling/test1_broad_rep2/test1_broad_rep2_peaks.broadPeak \
      --peak-list peak_calling/test1_broad/test1_broad_peaks.broadPeak \
      --input-file-type broadPeak \
      -o peak_calling/test1_broad.idr \
      --plot \
      --use-best-multisummit-IDR \
      --log-output-file peak_calling/test1_broad.idr.log
    
    # test2 IDR
    idr \
      --samples peak_calling/test2_broad_rep1/test2_broad_rep1_peaks.broadPeak \
               peak_calling/test2_broad_rep2/test2_broad_rep2_peaks.broadPeak \
      --peak-list peak_calling/test2_broad/test2_broad_peaks.broadPeak \
      --input-file-type broadPeak \
      -o peak_calling/test2_broad.idr \
      --plot \
      --use-best-multisummit-IDR \
      --log-output-file peak_calling/test2_broad.idr.log
    
    # test3 IDR
    idr \
      --samples peak_calling/test3_broad_rep1/test3_broad_rep1_peaks.broadPeak \
               peak_calling/test3_broad_rep2/test3_broad_rep2_peaks.broadPeak \
      --peak-list peak_calling/test3_broad/test3_broad_peaks.broadPeak \
      --input-file-type broadPeak \
      -o peak_calling/test3_broad.idr \
      --plot \
      --use-best-multisummit-IDR \
      --log-output-file peak_calling/test3_broad.idr.log
    
    ```
    
3. **Filter out unreliable peaks(not use)**
    
    ```python
    # test1
    awk -F'\t' '$12 > -log(0.05)/log(10)' peak_calling/test1.idr | cut -f1-10 > peak_calling/test1_idr.narrowPeak
    
    # test2
    awk -F'\t' '$12 > -log(0.05)/log(10)' peak_calling/test2.idr | cut -f1-10 > peak_calling/test2_idr.narrowPeak
    
    # test3
    awk -F'\t' '$12 > -log(0.05)/log(10)' peak_calling/test3.idr | cut -f1-10 > peak_calling/test3_idr.narrowPeak
    
    # board peaks-------
    # test1
    awk -F'\t' '$12 > -log(0.05)/log(10)' peak_calling/test1_broad.idr | cut -f1-10 > peak_calling/test1_broad_idr.broadPeak
    
    # test2
    awk -F'\t' '$12 > -log(0.05)/log(10)' peak_calling/test2_broad.idr | cut -f1-10 > peak_calling/test2_broad_idr.broadPeak
    
    # test3
    awk -F'\t' '$12 > -log(0.05)/log(10)' peak_calling/test3_broad.idr | cut -f1-10 > peak_calling/test3_broad_idr.broadPeak
    
    ```
    
4. filter blacklist 

1. put index to bam
    
    ```python
    samtools index mapped/bowtie2/test1.bam
    samtools index mapped/bowtie2/test1_rep1.bam
    samtools index mapped/bowtie2/test1_rep2.bam
    
    samtools index mapped/bowtie2/test2.bam
    samtools index mapped/bowtie2/test2_rep1.bam
    samtools index mapped/bowtie2/test2_rep2.bam
    
    samtools index mapped/bowtie2/test3.bam
    samtools index mapped/bowtie2/test3_rep1.bam
    samtools index mapped/bowtie2/test3_rep2.bam
    ```
    

1. bigwig file
    
    ```python
    mkdir -p bigwig_output
    
    bamCoverage -b mapped/bowtie2/test1.bam -o bigwig_output/test1.bw --normalizeUsing CPM
    bamCoverage -b mapped/bowtie2/test1_rep1.bam -o bigwig_output/test1_rep1.bw --normalizeUsing CPM
    bamCoverage -b mapped/bowtie2/test1_rep2.bam -o bigwig_output/test1_rep2.bw --normalizeUsing CPM
    
    bamCoverage -b mapped/bowtie2/test2.bam -o bigwig_output/test2.bw --normalizeUsing CPM
    bamCoverage -b mapped/bowtie2/test2_rep1.bam -o bigwig_output/test2_rep1.bw --normalizeUsing CPM
    bamCoverage -b mapped/bowtie2/test2_rep2.bam -o bigwig_output/test2_rep2.bw --normalizeUsing CPM
    
    bamCoverage -b mapped/bowtie2/test3.bam -o bigwig_output/test3.bw --normalizeUsing CPM
    bamCoverage -b mapped/bowtie2/test3_rep1.bam -o bigwig_output/test3_rep1.bw --normalizeUsing CPM
    bamCoverage -b mapped/bowtie2/test3_rep2.bam -o bigwig_output/test3_rep2.bw --normalizeUsing CPM
    
    ```
    
2. intersection
    
    ```python
    bedtools slop -i Annotate DESeq2.bed -g hg19.genome -l 1000 -r 0 -s > deg_genes_upstream1kb.bed
    
    ```