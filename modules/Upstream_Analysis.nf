/*
 *-----------------------------------------------------
 *    Splitting the Fastq File According to Barcode
 *                     ------fastq-multx
 *-----------------------------------------------------
 */

//Fastq File Splitting
 
process Fastq_Splitting{
        debug true
        tag "$Sample_ID"
        label 'fastq_multx'

        input:
          tuple val(Sample_ID),path(reads1),path(reads2),path(barcode_file)

        output:
          path '*R1_1.fq.gz', emit:Reads1
          path '*R2_1.fq.gz', emit:Reads2
          path '*R1_2.fq.gz', emit:reads1
          path '*R2_2.fq.gz', emit:reads2
          tuple val(Sample_ID), path("${Sample_ID}_Fastq_Multx_Result_1.txt"), path("${Sample_ID}_Fastq_Multx_Result_2.txt"), emit: Decode_Files 

        when:
          !params.skip_Fastq_Splitting

        script:
            """
          fastq-multx -B ${barcode_file} -b ${reads1} ${reads2} -o %${Sample_ID}.R1_1.fq.gz -o %${Sample_ID}.R2_1.fq.gz > ${Sample_ID}_Fastq_Multx_Result_1.txt
          fastq-multx -B ${barcode_file} -b ${reads2} ${reads1} -o %${Sample_ID}.R1_2.fq.gz -o %${Sample_ID}.R2_2.fq.gz > ${Sample_ID}_Fastq_Multx_Result_2.txt
          rm unmatched*
        """
}


//Fastq Files Decode Satuation
process Decode_Results {
        label 'fastq_multx'
        tag "$Decode_Files1"
        publishDir path : { params.Fastq_Splitting ? params.Outdir : "${params.Outdir}/Analyzed_Results/00.Decode_Satuation" },
                   mode : 'link',
                   overwrite : true
        
        input:
          tuple val(Sample_ID), path(Decode_Files1), path(Decode_Files2)


        output:
          path("${Sample_ID}_Decode_Results.csv")

        when:
          !params.skip_Fastq_Splitting
    
        script:
          """
          sed -i '\$d' $Decode_Files1
          sed -i '\$d' $Decode_Files2
          Rscript $baseDir/bin/Decode_Results.R $Decode_Files1 $Decode_Files2
          mv ./Decode_Results.csv ${Sample_ID}_Decode_Results.csv
        """
}


//Sample Names Collation
process SampleNames{
        label 'fastq_multx'
        tag "$Sample_Name"
    
        input:
          path Reads1
          path Reads2
          path reads1
          path reads2

        output:
          tuple val("${Sample_Name}"), path("${Sample_Name}.1R1.fq.gz"), path("${Sample_Name}.1R2.fq.gz"), emit:Fastq_R1
          tuple val("${Sample_Name}"), path("${Sample_Name}.2R1.fq.gz"), path("${Sample_Name}.2R2.fq.gz"), emit:Fastq_R2

        when:
          !params.skip_Fastq_Splitting
      
        script:
          Sample_Name=Reads1.name.split("\\.")[0]
          """
            mv $Reads1 ${Sample_Name}.1R1.fq.gz
            mv $Reads2 ${Sample_Name}.1R2.fq.gz
            mv $reads1 ${Sample_Name}.2R1.fq.gz
            mv $reads2 ${Sample_Name}.2R2.fq.gz
        """
}


//Fastq File Incorporation
process Fastq_Merge{
        label 'fastq_multx'
        tag "$Sample_Name"
        publishDir path : { params.skip_Fastq_Splitting ? params.Outdir : "${params.Outdir}/Analyzed_Data/00.Splitted_Data" },
                   mode : 'copy',
                   overwrite : true

        input:
          tuple val(Sample_Name), path(Reads1), path(Reads2)
          tuple val(Sample_Name), path(reads1), path(reads2)     

        output:
          tuple val(Sample_Name), path("${Sample_Name}.R1.fq.gz"),path("${Sample_Name}.R2.fq.gz")

        when:
          !params.skip_Fastq_Splitting

        script:
          """
          cat $Reads1 $reads2 > ${Sample_Name}.R1.fq.gz
          cat $Reads2 $reads1 > ${Sample_Name}.R2.fq.gz
        """
}



/*
 *-----------------------------------------------------
 *                Data Quality Control
 *                     ------Fastp & Trim_Galore
 *                     ------FastQC & Multiqc
 *-----------------------------------------------------
 */

        
process Fastp {
        label 'Data_Quality_Control'
        tag "$Sample_Name"
        cache false
        publishDir path : { params.skip_Fastp ? params.Outdir : "${params.Outdir}/Analyzed_Data/01.Clean_Data" },
                   mode : 'link'，
                   overwrite : true

        input:
          tuple val(Sample_Name), path(Reads1), path(Reads2)

        output:
          tuple val(Sample_Name), path("${Sample_Name}_fastq_R1.fq.gz"), path("${Sample_Name}_fastq_R2.fq.gz")

        when:
          !params.skip_QC && !params.skip_Fastp

        script:
          """
          fastp -g -x -Q -p 20 -w 16 -i $Reads1 -o ${Sample_Name}_fastq_R1.fq.gz -I $Reads2 -O ${Sample_Name}_fastq_R2.fq.gz -h ${Sample_Name}.html
        """
}


process Trim_Galore {
        label 'Data_Quality_Control'
        tag "$Sample_ID"
        publishDir path : { params.skip_Trim_Galore ? params.Outdir : "${params.Outdir}/Analyzed_Data/01.Clean_Data" },
                   mode : 'link',
                   overwrite : true

        input:
          tuple val(Sample_Name), path(Reads1), path(Reads2)

        output:
          tuple val(Sample_Name), path("${Sample_Name}_val_1.fq.gz"), path("${Sample_Name}_val_2.fq.gz")

        when:
          !params.skip_QC && !params.skip_Trim_Galore

        script:
          """
          trim_galore -q 25 --phred33 --stringency 3 --length 36 --fastqc --paired $Reads1 $Reads2 --gzip --basename ${Sample_Name} -o .
        """
}


process FastQC {
        label 'Data_Quality_Control'
        tag "$Sample_Name"
        publishDir path : { params.skip_FastQC ? params.Outdir : "${params.Outdir}/Analyzed_Data/02.Quality_Report" },
                   mode : 'link',
                   overwrite : true

        input:
          tuple val(Sample_Name), path(Reads1), path(Reads2)

        output:
          path("*.HTML")
          path("*.zip")

        when:
          !params.skip_QC && !params.skip_FastQC

        script:
          """
          fastqc -o . -t 40 ${Reads1}
          fastqc -o . -t 40 ${Reads2}
        """
}

              
/*
 *-----------------------------------------------------
 *                    Reads Mapping
 *                     ------Hisat2 & Bowtie2
 *                    Coverage Depth
 *                     ------Samtools
 *-----------------------------------------------------
 */

process Hisat2Align {
        label 'aligners'
        tag "$Sample_Name"
        publishDir path : { params.skip_Alignment ? params.Outdir : "${params.Outdir}/Analyzed_Data/02.Aligned_Data" },
                   mode : 'link'


        input:
        tuple val(Sample_Name), path(Reads1), path(Reads2)


        output:
        tuple val("${Sample_Name}"), path("${Sample_Name}.bam"), emit : Bam_File
        path("*_aligned.txt")


        when:
          !params.skip_Alignment


        script:
          """
          hisat2 --summary-file ${Sample_Name}_aligned.txt -t -x ${params.Hisat2_index} -p 40 \
          -1 $Reads1 -2 $Reads2 | samtools sort -@ 40 -o ${Sample_Name}.bam
        """
}
          
process Coverage {
        debug true
        label 'aligners'
        tag "$Sample_Name"
        publishDir path : { params.skip_Coverage ? params.Outdir : "${params.Outdir}/" },
                   mode : 'copy'


        input:
          tuple val(Sample_Name), path(Bam_file)


        output:
          path("*.txt")
          

        when:
          !params.skip_Coverage && !params.skip_Alignment


        script:
          """
          bedtools coverage -a ${params.Coverage_bed_file} -b ${Bam_file} > ${Sample_Name}_Coverage.txt 
        """
}

process Saturation {
        debug true
        label 'aligners'
        tag "$Sample_Name" 


        input:
          tuple val(Sample_Name), path(Bam_file)


        output:
          tuple path("${Sample_Name}.Gene_Number.txt"), path("${Sample_Name}.Reads_Number.txt")


        when:
         !params.skip_Saturation && !params.skip_Alignment


        script:
          """
          bash $baseDir/bin/Saturation.sh $Sample_name $Bam_file $params.GTF
  
        """
}

process Saturation_PlotS {
        debug true
        label 'aligners'
        tag "$Sample_Name"
        publishDir path : { params.skip_Saturation ? params.Outdir : "${params.Outdir}/Analyzed_Results/01.Saturation_Analysis" },
                   mode : 'link'        


        input:
          tuple path(Gene_Number_File), path(Reads_Number_File)


        output:
          path("*.png")


        when:
         !params.skip_Saturation && !params.skip_Alignment


        script:
          Sample_Name = Gene_Number_File.name.split("\\.")[0]
          """
          Rscript $baseDir/bin/Saturation_PlotS.R ${Gene_Number_File} ${Reads_Number_File}
          mv Saturation.png ${Sample_Name}_Saturation.png
        """
}


process Saturation_PlotA {
        debug true
        label 'aligners'
        tag "Saturation_Assembly"
        publishDir path : { params.skip_Saturation ? params.Outdir : "${params.Outdir}/Analyzed_Results/01.Saturation_Analysis" },
                   mode : 'link'


        input:
          path Num_Files


        output:
          path("*.png")


        when:
         !params.skip_Saturation && !params.skip_Alignment


        script:
          """
          Rscript $baseDir/bin/Saturation_PlotA.R $Num_Files 

        """
}
/*
 *-----------------------------------------------------
 *                   Gene Quantification
 *                     ------FeatureCounts
 *-----------------------------------------------------
 */

process FeatureCounts{
        label 'Gene_Quantification'
        tag "$Sample_Name"
        publishDir path : { params.skip_Transcriptome_Analysis ? params.Outdir : "${params.Outdir}/Analyzed_Data/03.Expression_Data" },
                   mode : 'link'


        input:
          tuple val(Sample_Name), path(Bam_file)


        output:
          path("*_counts.txt")


        when:
          !params.skip_Transcriptome_Analysis


        script: 
         """
          featureCounts -T 40 -f -p -B -C -t gene -g gene_id -a ${params.GTF} -o ${Sample_Name}_counts.txt $Bam_file
        """
}
