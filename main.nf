//HelpMessage
def helpMessage() {
        log.info"""
        Usage:

          nextflow <path/to/MeRIPseqPipe/main.nf> -c <path/to/MeRIPseqPipe/nextflow.config [options]

        OR,

          nextflow  run <path/to/MeRIPseqPipe> [options]


        Options(defaults in parentheses):


        Mandatory arguments:
          --Input_Fastq                         .csv format(table splited by commas) : Sample_ID, Reads1, Reads2, Barcode_File(Optional)
          --Outdir                              The output directory where the results will be saved ( Default : $baseDir/results )
 


        References
          --rRNA_fasta                          Path to HISAT2 rRNA fasta reference
          --Fasta                               Path to Fasta reference
          --GTF                                 Path to GTF reference
          


        Index
          --rRNA_index                          Path to HISAT2 rRNA fasta reference
          --Hisat2_index                        Path to HISAT2 index, eg. 'PATH/TO/Hisat2_index/*'
          --Bowtie2_index                       Path to Bowtie2 index, eg. 'PATH/TO/Bowtie2_index/*'
          --BWA_index                           Path to BWA index, eg. 'PATH/TO/BWA_index/*'



        Skipping modes Options
          --skip_QC                             Skip all Quality Control steps
          --skip_Fastp                          Skip fastp
          --skip_Trim_Galore                    Skip Trim Galore
          --skip_Fastqc                         Skip Fastqc
          --skip_Multiqc                        Skip Multiqc
          --skip_filterrRNA                     Skip the process of rRNA cleaning before mapping
          --skip_dPCR                           Skip the process of filtering PCR duplicates of bam files
          --skip_expression                     Skip gene expression analysis
          --skip_diffexpression                 Skip all Differential methylation analysis

 
  """.stripIndent() 
}


// Show help message
if (params.Help) {
        helpMessage()
        exit 0
}

//check Input Fastq Files
def HEADER = ["Sample_ID","Reads1","Reads2","Barcode_File"]
def header = ["Sample_ID","Reads1","Reads2"]


if ( params.Input_Fastq ) {
        File Input_Fastq = new File( params.Input_Fastq )
        if ( !Input_Fastq.exists() ) {
          exit 1,
          println( "ERROR : Input Fastq csv file not found : ${params.Input_Fastq} ")
        } else {
            Input_Fastq.eachLine { line, line_num ->
                def list = line.split(",")
                if (list.size() < 3) {
                    exit 1,
                    println( "ERROR : Invalid number of columns: Line ${line_num} " )
                }
                if (line_num == 1) {
                    if (list.size() == 3) {
                        if (list != header) {
                            exit 1,
                            println( "ERROR : Please check samplesheet header\nThe correct format is as follows\n[Sample_ID,Reads1,Reads2]" )
                        } else {
                            barcode = false
                        }
                    } else if (list.size() == 4) {
                        if (list != HEADER) {
                            exit 1,
                            println( "ERROR : Please check samplesheet header\nThe correct format is as follows\n[Sample_ID,Reads1,Reads2,Barcode_File]" )
                        } else {
                            barcode = true
                        }
                    }
                } else {
                      list.eachWithIndex { it, i ->
                          if (i > 0){
                              File Fastq_File = new File(it)
                              if ( !Fastq_File.exists() ) {
                                  exit 1,
                                  println( "Input file not found : ${it}")
                              }
                          }
                      }
                }
            }
        }
}

//Parameter Settings
if( params.skip_QC ){
        params.skip_Fastp                       = true
        params.skip_Trim_Galore                 = true
        params.skip_FastQC                      = true
        params.skip_MultiQC                     = true
}
if( params.skip_Alignment ){
        params.skip_Coverage                    = true
        params.skip_Saturation                  = true
        params.skip_Transcriptome_Analysis      = true
        params.skip_Methylation_Analysis        = true
        params.skip_Combined_Analysis           = true
}

if( params.skip_Transcriptome_Analysis ){
        params.skip_Diffexpression_Analysis     = true
}

if ( params.skip_Diffexpression_Analysis ){
        params.skip_DESeq2                      = true
        params.skip_edgeR                       = true
        params.Diff_Expression                  = "none"
}


//Module Importing
include {Fastq_Splitting; Decode_Results; SampleNames; Fastq_Merge; Fastp; Trim_Galore; Hisat2Align; Coverage; Saturation; Saturation_PlotS; Saturation_PlotR; FeatureCounts;} from './modules/Upstream_Analysis'
include {Expressed_Gene_Counts; Counts_Matrix_Build; DESeq2; edgeR} from './modules/Expression_Analysis'


//Workflow Establishment
//FASTQ Files Splitting(When the Barcode File is specified)
workflow Fastq_Multx {
         take:
         samplesheet

         main:
          Fastq_Splitting(samplesheet)
          SampleNames(Fastq_Splitting.out.Reads1.flatten(),
                      Fastq_Splitting.out.Reads2.flatten(),
                      Fastq_Splitting.out.reads1.flatten(),
                      Fastq_Splitting.out.reads2.flatten())
          Fastq_Merge(SampleNames.out[0],SampleNames.out[1])
          Decode_Results(Fastq_Splitting.out.Decode_Files)

        emit:
          Raw_data=Fastq_Merge.out
}

//Date Preprocessing & Quality Control
workflow QC {
        take:
        Raw_data


        main:
        Fastp(Raw_data)
        Trim_Galore(Fastp.out)


        emit:
        Align_Input=Trim_Galore.out[0]

}

//Alignment & Saturation Analysis
workflow Align{
        take:
        Clean_data


        main:
        Hisat2Align(Clean_data)
        Saturation(Hisat2Align.out.Bam_File)
        Saturation_PlotA(Saturation.out.collect())
        Saturation_PlotS(Saturation.out)


        emit:
        Hisat2Align.out[0]
}

//Expression Quantification
workflow Gene_Quan{
        take:
          Bam_Files


        main:
          FeatureCounts(Bam_Files)
          Expressed_Gene_Counts(FeatureCounts.out.collect())


        emit:
          FeatureCounts.out
}

//Differential Expression Analysis
workflow Diff_Express{
        take:
          Counts_files

        main:
          Counts_Matrix_Build(Counts_files)
          if ( params.Diff_Expression == "DEseq2" ){
            DESeq2(Counts_Matrix_Build.out)
          } else if ( params.Diff_Expression == "edgeR" ){
            edgeR(Counts_Matrix_Build.out)
          }
}


//Overall Workflow
workflow {
        if (barcode) {
                Fastq_Multx(channel.fromPath(params.Input_Fastq)
                                   .splitCsv(header: true, sep: ",", strip: true)
                                    .map {
                                        row ->
                                     [ row.Sample_ID, row.Reads1, row.Reads2, row.Barcode_File ]
                                    }
                           .unique())
                QC(Fastq_Multx.out.Raw_data)
        } else {
                QC(channel.fromPath(params.Input_Fastq)
                          .splitCsv(header: true, sep: ",", strip: true)
                          .map {
                                row ->
                                 [ row.Sample_ID, row.Reads1, row.Reads2 ]
                                }
                          .unique())
                }
        Align(QC.out.Align_Input)
        Gene_Quan(Align.out)
        Diff_Express(Gene_Quan.out.collect())
}
