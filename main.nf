#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//--------------------------------------------------------------------------
// Param Checking
//--------------------------------------------------------------------------

if(!params.fastaSubsetSize) {
  throw new Exception("Missing params.fastaSubsetSize")
}

if(params.inputFilePath) {
  seqs = Channel.fromPath( params.inputFilePath )
           .splitFasta( by:params.fastaSubsetSize, file:true  )
}
else {
  throw new Exception("Missing params.inputFilePath")
}

//--------------------------------------------------------------------------
// Main Workflow
//--------------------------------------------------------------------------

workflow {
    signalp5(seqs)

    filterInputFastaByResults(params.inputFilePath, signalp5.out.pred_summary.collectFile())

    // Now that we have our filtered result, we can split it further to run with signalp4 and signalp6
    filteredSeq = filterInputFastaByResults.out.fasta.splitFasta( by:200, file:true )

    signalp4(filteredSeq)
    signalp6(filteredSeq, filterInputFastaByResults.out.recordCount.first())

    collectedGff = signalp5.out.gff.mix(signalp4.out.gff).mix(signalp6.out.gff).collectFile(name: 'result.gff3')

    indexResults(collectedGff, params.outputFileName)
}

/**
 NOTE:  there is a potential optimization here.  We could trim the first 70-100 amino acids
 from the N terminus and get the unique set of results.  (the def line would need to be
 a list of proteins which match the sequence or map back to proteins in some other way
*/
process filterInputFastaByResults {
    container = 'bioperl/bioperl:stable'

    input:
    path fasta
    path predictions

    output:
    path "filtered.fasta", emit: fasta
    env recordCount, emit: recordCount

    script:
    """
    filterProteinsByScore.pl --fasta $fasta \
        --predictions $predictions \
        --score_cutoff ${task.ext.filter_score_cutoff} \
        --pct_proteins_cutoff ${task.ext.filter_min_protein_percent_cutoff} \
        --output_file filtered.fasta


    recordCount=\$(grep -c "^>" filtered.fasta)
    """
}


process signalp6 {
    label = "signalp"

    input:
    path subsetFasta
    val recordCount

    output:
    path "combined.gff3", emit: gff


    when:
    Integer.parseInt(recordCount) <= task.ext.protein_count_limit

    script:
    """
    signalp6 --fastafile $subsetFasta \
        --format none \
        --organism ${task.ext.org} \
        --mode fast \
        --output_dir .

    fixAndCombineGff.pl --gff output.gff3 --region_gff region_output.gff3 --sp_version 6 --output_file combined.gff3
    """
}


process signalp4 {
    label = "signalp"

    input:
    path subsetFasta

    output:
    path "signalp4.gff3", emit: gff

    script:
    """
    signalp4  -f short \
        -t ${task.ext.org} \
        -n signalp4.gff2 \
        $subsetFasta >sp4_prediction_summary.txt

    # make gff3 format (remove the group column)
    fixAndCombineGff.pl --gff signalp4.gff2 --sp_version 4 --output_file signalp4.gff3
    """
}




process signalp5 {
    label = "signalp"

    input:
    path subsetFasta

    output:
    path 'sp5_prediction_summary.txt',  emit: pred_summary
    path "subsetFasta.gff3", emit: gff


    script:
    """
    signalp -fasta $subsetFasta \
        -format short \
        -gff3 \
        -org ${task.ext.org} \
        -plot 'none' \
        -prefix signalp5 \
        -stdout >sp5_prediction_summary.txt

    fixAndCombineGff.pl --gff signalp5.gff3 --sp_version 5 --output_file subsetFasta.gff3
    """
}

// make this a module??
process indexResults {
  container = 'biocontainers/tabix:v1.9-11-deb_cv1'
  publishDir params.outputDir, mode: 'copy'

  input:
    path gff
    val outputFileName

  output:
    path '*.gff.gz'
    path '*.tbi'

  script:
  """
  sort -k1,1 -k4,4n $gff > $outputFileName
  bgzip $outputFileName
  tabix -p gff ${outputFileName}.gz
  """
}
