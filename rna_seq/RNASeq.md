## <center>RNA-Seq<center>

### Preparation
- ***create env***
  ```bash
  mamba create -n RNASeq python=3.7
  ```

- ***fastqc***
  ```bash
  mamba inatall fastqc
  fastqc -t 24 -o ./FastQC_result.fix/ -q ./fix.fastq/test_*.gz &
  ```

- ***cutadapt***
  ```bash
  cutadapt -j 6 --times 1  -e 0.1  -O 3  --quality-cutoff 25  -m 55 \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
  -o fix.fastq/test_R1_cutadapt.fq.gz \
  -p fix.fastq/test_R2_cutadapt.fq.gz \
  raw.fastq/test_R1.fq.gz \
  raw.fastq/test_R2.fq.gz > fix.fastq/test_cutadapt.temp.log 2>&1 &
  ```

### build index
- ***bowtie2***
  ```bash
  mamba install bowtie2
  bowtie2-build -t 6 ref_hg38.fa ref_hg38.fa > bt2_index.log 2>&1 &
  ```

- ***STAR***
  ```bash
  STAR --runThreadN 12 --runMode genomeGenerate --genomeDir STAR_index --genomeFastaFiles STAR_index/ref_hg38.fa --sjdbGTFfile gtf/hg38_refseq.gtf --sjdbOverhang 100 &
  ```

- ***hisat2***
  ```bash
  # without fasta and snp
  hisat2-build -p 10 ref_hg38.fa ref_hg38.fa > hisat2_build.log 2>&1 &

  # snp  
  hisat2_extract_exons.py hg38_refseq.gtf > hg38_refseq.exon &  # make exon
  hisat2_extract_splice_sites.py hg38_refseq.gtf > hg38_refseq.ss &  # make splice site
  hisat2_extract_snps_haplotypes_UCSC.py ref_hg38.fa snp151Common.txt snp151Common &  # make snp and haplotype
  hisat2-build -p 6 --snp snp151Common.snp --haplotype snp151Common.haplotype --exon hg38_refseq.exon  --ss hg38_refseq.ss ref_hg38.fa ref_hg38.fa.snp_gtf > hisat2_build.log 2>&1 &  # build index
  ```

### alignment
- ***tophat2***
  ```bash
  mamba create -n python2 python=2
  mamba install tophat
  tophat2 -o ./tophat2 -p 20 -G reference/gtf/hg38_NCBI_RefSeq.gtf reference/\ bowtie2_index/ref_hg38.fa ./fix.fastq/test_R1_cutadapt.fq.gz ./fix.fastq/\ test_R2_cutadapt.fq.gz > ./tophat2/tophat2.log 2>&1 &
  # view header 
  samtools view -H accepted_hits.bam
  # view
  samtools view accepted_hits.bam | less -S 
  ```

- ***STAR***
  ```bash  
  STAR 
  --genomeDir reference/STAR_index 
  --runThreadN 10 --readFilesIn fix.fastq/test_R1_cutadapt.fq.gz fix.fastq/test_R2_cutadapt.fq.gz 
  --readFilesCommand zcat
  --outFileNamePrefix test_STAR 
  --outSAMtype BAM Unsorted 
  --outSAMstrandField intronMotif 
  --outSAMattributes All 
  --outFilterIntronMotifs RemoveNoncanonical > test_STAR/test_STAR.log 2>&1 &
  ```

- ***hisat2***
  ```bash
  # without snp and gtf
  hisat2 -p 12 \
  -x ../reference/hisat2_index/ref_hg38.fa \
  -1 ../fix.fastq/test_R1_cutadapt.fq.gz \
  -2 ../fix.fastq/test_R2_cutadapt.fq.gz \
  -S ./test_hisat2.sam > ./test_hisat2.log 2>&1 &
  # with snp and gtf
  hisat2 -p 6 \
  -x /Users/meng/ngs_course/reference/hisat2_index/ref_hg38.fa.snp_gtf \
  -1 ./fix.fastq/test_R1_cutadapt.fq.gz \
  -2 ./fix.fastq/test_R2_cutadapt.fq.gz \
  -S ./bam/test_hisat2.sam > ./bam/test_hisat2.log 2>&1 &
  ```
  ```bash
  # samtools convert sam to bam
  samtools sort -O BAM -o test_hisat2.sort.bam -@ 10 -m 2G -T test_hisat2.sort.bam.tmp test_hisat2.sam &
  # view
  samtools view -h test_hisat2.sort.bam | less -S
  # index
  samtools index test_hisat2.sort.bam
  # view region
  samtools view test_hisat2.sort.bam chr1:100000000-110000000 | less -S
  # idx view
  samtools idxstats test_hisat2.sort.bam
  ```

### count by HTSeq and featureCount
- ***htseq***
  ```bash
  mamba install htseq
  mamba install subread

  htseq-count -f bam -r pos \
  --max-reads-in-buffer 1000000 \
  --stranded no \
  --minaqual 10 \
  --type exon --idattr gene_id \
  --mode union \
  --nonunique none \
  --secondary-alignments ignore \
  --supplementary-alignments ignore \
  --counts_output ./count_result/test_count.tsv \
  --nprocesses 1 \
  ./test_hisat2_snp/test_hisat2.sort.bam \
  ./reference/gtf/hg38_refseq.gtf > ./count_result/test_count.HTSeq.log  2>&1 & 
  ```
  
- ***featureCount***
  ```bash
  featureCounts -t exon -g gene_id \
  -Q 10 --primary -s 0 -p -T 6 \
  -a ../reference/gtf/hg38_refseq.gtf \
  -o ./test_count.featureCounts \
  ../test_hisat2_snp/test_hisat2.sort.bam > ./test_count.featureCounts.log  2>&1 &
  ```

### RNA Quantification
- ***cuffdiff***
  ```bash
  mamba install cufflinks
  cuffdiff -o cuffdiff_result --labels Ctrl,METTL3_KD -p 12 \
  --min-alignment-count 10 \
  --min-reps-for-js-test 2 \
  --library-type fr-unstranded \
  --dispersion-method pooled \
  --library-norm-method geometric \
  ./reference/gtf/hg38_refseq_from_ucsc.rm_XM_XR.fix_name.gtf \
  ./bam.STAR/Ctrl_rep1.bam,./bam.STAR/Ctrl_rep2.bam \
  ./bam.STAR/KD_rep1.bam,./bam.STAR/KD_rep2.bam \
  > ./cuffdiff_result/cuffdiff_run.log 2>&1 &
  ```