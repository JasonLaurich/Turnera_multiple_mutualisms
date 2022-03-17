# Turnera_multiple_mutualisms
This is a repository for data and code pertaining to the manuscript "The genetic architecture of multiple mutualisms and mating system in Turnera ulmifolia"

Below are a description of each file included in this repository. For ease of interpretation, I have broken code and data files into conceptual sections.

1. "Keystone" data summarizing key information or used in multiple analyses.

  a. Trait and seedset data.xlsx: an excel file containing phenotype data for all _Turnera ulmifolia_ plants. Second tab includes data from seedset experiment.

  b. pheno.txt: text file containing data for all plants. This file is uploaded in the following R codes : GMatrix_Estimation_final.R, Population_Mean_Trait_Variation.R.

  c. seedset.txt: text file containing data for the seedest experiment. This file is uploaded in the following R codes : GMatrix_Estimation_Final.R.

  d. Pedigree.xlsx: excel file containing pedigree information for all _Turnera ulmifolia_ plants. 

  e. pedmat.txt: text file containing pedigree information. This file is used in the following R codes : GMatrix_Estimation_final.R, Bayesian_Heritability_final.R
  
2. Fitting G matrices and univariate Bayesian objects

  a. GMatrix_Estimation_final.R: (sample) code that fits G matrices and univariate Bayesian objects for all populations of _T. ulmifolia_. Also includes code for   construction of G matrix for seedset experiment. 
  
  b. b. Bayesian_Heritability_final.R: (sample) code that fits univariate Bayesian objects to randomized datasets for each population of _T. ulmifolia_ for the purposes of significance testing.
  
3. Miscellaneous
  a. Population_Mean_Trait_Variation.R: code that generates Figure 3, which displays variation in mean trait values among populations.
