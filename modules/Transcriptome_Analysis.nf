*
 *-----------------------------------------------------
 *         Expression Quantification Analysis
 *-----------------------------------------------------
 */
process Expressed_Gene_Counts{
        label 'Transcriptome_Analysis'
        tag "$Counts_file"
        publishDir path : { params.skip_Transcriptome_Analysis ? params.Outdir : "${params.Outdir}/Analyzed_Results/02.Expression_Analysis/01.Gene_Quantification" },
                   mode : 'link'

        input:
          path Counts_file


        output:
          path("Gene_Counts.csv")         


        script:
          """
          Rscript $baseDir/bin/Expression_Quan.R $Counts_file
        """
}


/*
 *-----------------------------------------------------
 *          Differential Expression Analysis
 *-----------------------------------------------------
 */

process Counts_Matrix_Build {
        label 'Transcriptome_Analysis'
        tag "Counts_file"
        publishDir path : { params.skip_Transcriptome_Analysis ? params.Outdir : "${params.Outdir}/Analyzed_Results/02.Expression_Analysis/01.Gene_Quantification" },
                   mode : 'link'

        input:
          path Counts_file


        output:
          path("Counts_Matrix.csv")


         script:
          """
          Rscript $baseDir/bin/Counts_Matrix.R $Counts_file
        """
}
