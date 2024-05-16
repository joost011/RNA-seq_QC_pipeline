nextflow.enable.dsl=2


process createCovariateMatrix {
  publishDir "${params.outDir}/1_create_covariate_matrix/", mode: 'copy'

  time '6h'
  memory '12 GB'
  cpus 1

  input:
  path inputDir
  
  output:
  path "${params.cohortName}_covariates.txt.gz", emit: matrix
  path "failed_samples.txt"
  
  script:
  """
  # 1. Create covariate matrix
  create_covariate_matrix.py --input_dir ${params.inputDir} \
  --cohort_name ${params.cohortName}

  # 2. Gzip the matrix
  gzip ${params.cohortName}_covariates.txt
  """
}

process normalizeCovariateMatrix {
  publishDir "${params.outDir}/2_normalize_covariate_matrix/", mode: 'copy'

  time '6h'
  memory '6 GB'
  cpus 1

  input:
  path covariateMatrix
  
  output:
  path "${params.cohortName}_covariates_normalized.txt.gz"
  
  script:
  """
  # 1. Normalize covariate matrix
  normalize_covariate_matrix.py --matrix_path ${covariateMatrix} \
  --cohort_name ${params.cohortName}
  
  # 2. Gzip the matrix
  gzip ${params.cohortName}_covariates_normalized.txt
  """
}

process createGeneCountsMatrix {
  publishDir "${params.outDir}/3_create_gene_counts_matrix/", mode: 'copy'

  time '6h'
  memory '12 GB'
  cpus 1

  input:
  path inputDir
  
  output:
  path "${params.cohortName}_gene_counts.txt.gz"
  
  script:
  """
  # 1. Create gene count matrix
  create_gene_counts_matrix.py --input_dir ${params.inputDir} \
  --cohort_name ${params.cohortName}

  # 2. Gzip the matrix
  gzip ${params.cohortName}_gene_counts.txt
  """
}

process geneCountsTMMNormalization {
  publishDir "${params.outDir}/4_gene_counts_tmm_normalization/", mode: 'copy'

  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path geneCountMatrix
  
  output:
  path "${params.cohortName}_gene_counts-TMM.txt.gz", emit: tmmNormalizedMatrix
  path "${params.cohortName}_gene_counts-CPM.txt.gz"
  
  script:
  """
  # 1. Load RPlus module
  ml RPlus

  # 2. Run TMM normalization
  calculate_TMM.R -c ${geneCountMatrix} \
  -t ./${params.cohortName}_gene_counts-TMM.txt \
  -o ./${params.cohortName}_gene_counts-CPM.txt

  # 3. Gzip output
  gzip *.txt
  """
}

process geneCountsPCAOutlierIdentification {
  publishDir "${params.outDir}/5_gene_counts_pca_outlier_identification/", mode: 'copy'

  time '6h'
  memory '12 GB'
  cpus 1

  input:
  path geneCountMatrix
  
  output:
  path 'pc1_2.txt'
  path 'pc1_2_no_outliers.txt'
  path 'passed_samples.txt', emit: passedSamples
  
  shell:
  '''
  # 1. Load Java module
  ml Java

  # 2. Do PCA
  java -Xmx20g -Xms20g -jar !{baseDir}/bin/eqtl-mapping-pipeline.jar	\
      --mode normalize \
      --in !{geneCountMatrix} \
      --out $(pwd) \
      --logtransform \
      --qqnorm \
      --adjustPCA \
      --maxnrpcaremoved 0 \
      --stepsizepcaremoval 0

  # 3. Write first PCs out to file  
  zcat *.QuantileNormalized.Log2Transformed.PCAOverSamplesEigenvectors.txt.gz \
  | awk 'BEGIN {OFS="\\t"} NR == 1 {print "Sample","Comp1","Comp2"} NR \
  > 1 {print $1,$2,$3}' > pc1_2.txt

  # 4. Remove uneccsary files 
  rm *gz

  # 5. Identify outliers based on z-score
  identify_pca_outliers.py --pc_file !{geneCountMatrix} \
  --z_score_threshold !{params.pcaOutlierZScoreThreshold}
  '''
}

process geneCountsPCAOutlierRemoval {
  publishDir "${params.outDir}/6_gene_counts_pca_outlier_removal/", mode: 'copy'

  time '6h'
  memory '12 GB'
  cpus 1

  input:
  path includedSamples
  path geneCountMatrix
  
  output:
  path "${params.cohortName}_gene_counts-TMM.SampleSelection.ProbesWithZeroVarianceRemoved.txt.gz", emit: matrixNoZeroVariance
  path "${params.cohortName}_gene_counts-TMM.SampleSelection.txt.gz"
  
  script:
  """
  # 1. Load Java module
  ml Java

  # 2. Remove outliers
  java -jar ${baseDir}/bin/eqtl-mapping-pipeline.jar \
    --mode normalize \
    --sampleInclude ${includedSamples} \
    --out ./ \
    --in ${geneCountMatrix}
  """
}

process geneCountsLog2Transform {
  publishDir "${params.outDir}/7_gene_counts_log_2_transform/", mode: 'copy'

  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path geneCountMatrix
  
  output:
  path "${params.cohortName}_gene_counts-TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.txt.gz"
  
  script:
  """
  # 1. Load Java module
  ml Java

  # 2. Perform Log2Transormation
  java -Xmx20g -Xms20g -jar ${baseDir}/bin/eqtl-mapping-pipeline.jar	\
      --mode normalize \
      --in ${geneCountMatrix} \
      --out ./ \
      --logtransform
  """
}

process geneCountsForceNormalDistribution {
  publishDir "${params.outDir}/8_gene_counts_force_normal_distribution/", mode: 'copy'

  time '6h'
  memory '12 GB'
  cpus 1

  input:
  path geneCountMatrix
  
  output:
  path "gene_counts-TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.forcenormal.txt.gz"
  
  script:
  """
  force_normal_distribution.py ${geneCountMatrix} \
  ./gene_counts-TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.forcenormal.txt.gz
  """
}

process geneCountsCovariateCorrection {
  publishDir "${params.outDir}/9_covariate_correction/", mode: 'copy'

  time '6h'
  memory '12 GB'
  cpus 1

  input:
  path geneCountMatrix
  path covariateMatrix
  
  
  output:
  path "*.txt"
  path "*.txt.gz"
  path "${params.cohortName}_gene_counts-TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.forcenormal.covariatecorrected.txt.gz.CovariatesRemovedOLS.txt.gz", emit: geneCountMatrix

  script:
  """
  # 1. Covariate PCA 
  pca.py ${covariateMatrix} covariates

  # 2. Extract and transpose PCs that explain variance greater than a certain threshold
  python3 extract_and_transpose_pcs.py --covariate_pcs covariates_PCs.txt \
    --covariates_explained_variance explainedVariance.txt \
    --explained_variance_threshold 0.99 

  # 3. Gzip the covariate matrix
  gzip covariates-pca_PCs-filtered-transpose.txt

  # 4. Load Java module
  ml Java

  # 5. Regress covariates
  java -jar ${baseDir}/bin/2023-12-08-Regression.jar \
    ${geneCountMatrix} \
    covariates-pca_PCs-filtered-transpose.txt.gz \
    ./${params.cohortName}_gene_counts-TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.forcenormal.covariatecorrected.txt.gz
  """
}

process finalPCA {
  publishDir "${params.outDir}/10_final_pca/", mode: 'copy'

  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path geneCountMatrix
  
  output:
  path "*.txt"
  path "*.png"
  
  script:
  """
  python pca.py ${geneCountMatrix} ./geneExp-PCA
  """
}

workflow {
  
  // Create and normalize covariate matrix
  createCovariateMatrix(params.inputDir)
  normalizeCovariateMatrix(createCovariateMatrix.out)

  // Create and normalize gene count matrix
  createGeneCountsMatrix(params.inputDir)
  geneCountsTMMNormalization(createGeneCountsMatrix.out)
  geneCountsPCAOutlierIdentification(geneCountsTMMNormalization.out.tmmNormalizedMatrix)
  geneCountsPCAOutlierRemoval(geneCountsPCAOutlierIdentification.out.passedSamples, geneCountsTMMNormalization.out.tmmNormalizedMatrix)
  geneCountsLog2Transform(geneCountsPCAOutlierRemoval.out.matrixNoZeroVariance)
  geneCountsForceNormalDistribution(geneCountsLog2Transform.out)
  geneCountsCovariateCorrection(geneCountsForceNormalDistribution.out, normalizeCovariateMatrix.out)
  finalPCA(geneCountsCovariateCorrection.out.geneCountMatrix)
}