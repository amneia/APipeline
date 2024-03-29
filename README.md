# APipeline
A Simple Process for Combined Methylation and Transcriptome Analysis
## Introduction

## Installation
[Nextflow](https://nf-co.re/docs/usage/installation)(version >= 22.03.0) and [Docker](https://docs.docker.com/engine/install/) need to be installed to run the pipeline.In addition to this,all dependencies are automatically deployed.
1. Install `nextflow`
  ```
  conda install -c bioconda nextflow
  ```
2.  Pull `docker` Image from Dockerhub :

3.  Cloning this Repository
  ```
  git clone
  ```
4.  Printing Help Information
  ```
  nextflow run 
  ```
5.  Test it in a minimal test dataset

6.  Start Running 

## Samplesheet input
### `--Input_Fastq`
`Input Fastq Files` is just like the following table splited by `,` separated,which is `.csv` suffix file.
#### Example for Transcriptome Analysis
|Sample_ID|Reads1|Reads2|Barcode_File|
|---------|------|------|------------|
|Sample1|SampleDir/Sample1_R1.fq.gz|SampleDir/Sample1_R2.fq.gz|SampleDir/Barcode.txt|
|Sample2|SampleDir/Sample2_R1.fq.gz|SampleDir/Sample2_R2.fq.gz|SampleDir/Barcode.txt|
|Sample3|SampleDir/Sample3_R1.fq.gz|SampleDir/Sample3_R2.fq.gz|SampleDir/Barcode.txt|
|Sample4|SampleDir/Sample4_R1.fq.gz|SampleDir/Sample4_R2.fq.gz|SampleDir/Barcode.txt|

#### Example for Methylation Analysis
|Sample_ID|IP_Reads1|IP_Reads2|INPUT_Reads1|INPUT_Reads2|
|---------|---------|---------|------------|------------|
|Sample1|SampleDir/Sample1_IP_R1.fq.gz|SampleDir/Sample1_IP_R2.fq.gz|SampleDir/Sample1_INPUT_R1.fq.gz|SampleDir/Sample1_INPUT_R2.fq.gz|
|Sample2|SampleDir/Sample2_IP_R1.fq.gz|SampleDir/Sample2_IP_R2.fq.gz|SampleDir/Sample2_INPUT_R1.fq.gz|SampleDir/Sample2_INPUT_R2.fq.gz|
|Sample3|SampleDir/Sample3_IP_R1.fq.gz|SampleDir/Sample3_IP_R2.fq.gz|SampleDir/Sample3_INPUT_R1.fq.gz|SampleDir/Sample3_INPUT_R2.fq.gz|
|Sample4|SampleDir/Sample4_IP_R1.fq.gz|SampleDir/Sample4_IP_R2.fq.gz|SampleDir/Sample4_INPUT_R1.fq.gz|SampleDir/Sample4_INPUT_R2.fq.gz|

### `--CompareFiles`
#### Example
|Control|Treat|
|-------|-----|
|Sample1|Sample3|
|Sample2|Sample4|

## Parameters
### Default Parameters
Software Version and Parameter Selection

*Fastq Files Splitting(Optional)

|Tools|Version|Parameters|
|-----|-------|----------|
|Fastq Multx|1.4.2|`fastq-multx -B ${Barcode_File} -b ${Reads1} ${Reads2} -o %${Sample_ID}.R1_1.fq.gz -o %${Sample_ID}.R2_1.fq.gz > ${Sample_ID}_Fastq_Multx_Result_1.txt`<br>`fastq-multx -B ${Barcode_File} -b ${Reads2} ${Reads1} -o %${Sample_ID}.R1_2.fq.gz -o %${Sample_ID}.R2_2.fq.gz > ${Sample_ID}_Fastq_Multx_Result_2.txt`


* Quality Control

|Tools|Version|Parameters|
|-----|-------|----------|
|`fastp`|0.23.4|`fastp -g -x -Q -p 20 -w 16 -i $Reads1 -o ${Sample_Name}_fastq_R1.fq.gz -I $Reads2 -O ${Sample_Name}_fastq_R2.fq.gz -h ${Sample_Name}.html`|
|`Trim Galore`|0.6.10|`trim_galore -q 25 --phred33 --stringency 3 --length 36 --fastqc --paired $Reads1 $Reads2 --gzip --basename ${Sample_Name} -o .`|
|`FastQC`|0.12.1||

* Alignment

|Tools|Version|Parameters|
|-----|-------|----------|
|`hisat2`|2.2.1|`hisat2 --summary-file ${Sample_name}_aligned.txt -t -x ${params.Hisat2_index} -p ${task.cpu} -1 $Reads1 -2 $Reads2 -S ${Sample_Name}.sam`|
|`Bowtie2`|2.5.2|`bowtie2 -p ${task.cpu} -x ${params.Bowtie2_index} -1 $Reads1 -2 $Reads2 -S ${Sample_Name}.sam 2>${Sample_Name}_aligned.txt`|
|`bwa`|0.7.17-r1188|`bwa bwasw -t ${task.cpu} ${params.BWA_index} $Reads1 $Reads2 > ${Sample_Name}.sam`|
|`Samtools`|1.8|`samtools sort -@ 40 -o ${Sample_name}.bam ${Sample_Name}.sam`|


*PCR Duplication(Optional)

|Tools|Version|Parameters|
|-----|-------|----------|
|`Samtools`|1.8|`samtools sort -n -@ ${task.cpu} $Bam -o ${Sample_Name}_nsorted.bam`<br>`samtools fixmate -@ ${task.cpu} -m ${Sample_Name}_nsorted.bam ${Sample_Name}_fixmate.bam`<br>`samtools sort -@ ${task.cpu} ${Sample_Name}_fixmate.bam -o ${Sample_Name}_fixmatesorted.bam`<br>`samtools markdup -r -@ ${task.cpu} IP_${Sample_Name}_fixmatesorted.bam IP_${Sample_Name}_rmdup.bam`|

* Coverage Analysis

|Tools|Version|Parameters|
|-----|-------|----------|
|`bedtools`|2.31.0|`bedtools coverage -a ${params.Coverage_bed_file} -b $Bam > ${Sample_Name}_Coverage.txt`|

* Transcriptome Analysis

|Tools|Version|Parameters|
|-----|-------|----------|
|`featureCounts`|2.0.6|`featureCounts -T ${task.cpu} -f -p -B -C -t gene -g gene_id -a ${params.GTF} -o ${Sample_Name}_counts.txt $Bam`|

* Methylation Analysis

### Specify the parameters as follows:
1. Edit the `nextflow.config`
2. Specify on the command line
```
nextflow main.nf -c nextflow.config --Input_Fastq samplesheet.csv
```


###
