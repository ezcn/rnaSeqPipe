#### 1. Download reference genome and index fasta file [hg38p12](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)

#### 2. Align reads to reference genome


```
~/bin/bwa mem -t 4 /mpbastudies3/IMMA/hg38p12/hg38.p12.fa /mpbastudies3/lonardo/raw/fastq/CTR-1/253_CTR_1_S1_R1_001.fastq.gz /mpbastudies3/lonardo/raw/fastq/CTR-1/253_CTR_1_S1_R2_001.fastq.gz | /mpba0/mpba-sw/samtools view -b > /mpbastudies3/lonardo/raw/bam/CTR_1.raw.bam

```

#### 3. Sort and index bam file

```
~/bin/sambamba sort -t 1 -m 32G -p --tmpdir /scratch -o /mpbastudies3/lonardo/raw/bam/${id}.${num}.bam /mpbastudies3/lonardo/raw/bam/${id}.${num}.raw.bam
```
#### 4. Download annotation gtf file [hg38 gtf](ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz)

#### 5. Count reads 

```
singularity exec -B /mpbastudies3/lonardo/:/data /mpba0/mpba-sw/biocontainers/featureCounts.img featureCounts -p -B -C -t exon -g gene_id -a /data/annotation/Homo_sapiens.GRCh38.95.gtf -F GTF -o /data/countsLonardo.txt /data/raw/bam/CTR.1.bam /data/raw/bam/CTR.2.bam /data/raw/bam/CTR.3.bam /data/raw/bam/NODAL.1.bam /data/raw/bam/NODAL.2.bam /data/raw/bam/NODAL.3.bam /data/raw/bam/POS.1.bam /data/raw/bam/POS.2.bam /data/raw/bam/POS.3.bam /data/raw/bam/NEG.1.bam /data/raw/bam/NEG.2.bam /data/raw/bam/NEG.3.bam

```

#### 6. Create samples and contrast files

```
samples.tsv
CTR	CTR.1
CTR	CTR.2
CTR	CTR.3
NODAL	NODAL.1
NODAL	NODAL.2
NODAL	NODAL.3
POS	POS.1
POS	POS.2
POS	POS.3
NEG	NEG.1
NEG	NEG.2
NEG	NEG.3
```
```
contrast.tsv
NODAL CTR
POS	NEG
```

#### 7. Run perl script (edgeR)

```
/usr/bin/perl /mpbastudies3/lonardo/run_DE_analysis.pl  --matrix /mpbastudies3/lonardo/Lonardo_counts.txt  --method edgeR --dispersion 0.1 --samples_file /mpbastudies3/lonardo/samples_lonardo.tsv --contrasts /mpbastudies3/lonardo/contrast_lonardo.tsv  --output /mpbastudies3/lonardo/results
```


#### 8. Post edgeR
