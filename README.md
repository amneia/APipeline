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
|Tools|Version|Parameters|
|-----|-------|----------|
|Qualty Control|
|--------------|
|`fastp`|0.23.4|fastp -g -x -Q -p 20 -w 16 -i $Reads1 -o ${Sample_Name}_fastq_R1.fq.gz -I $Reads2 -O ${Sample_Name}_fastq_R2.fq.gz -h ${Sample_Name}.html|
|`FastQC`|0.12.1|
|



### Specify the parameters as follows:
1. Edit the `nextflow.config`
2. Specify on the command line
```
nextflow main.nf -c nextflow.config --Input_Fastq samplesheet.csv
```


###
