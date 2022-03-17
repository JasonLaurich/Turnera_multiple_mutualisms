# Turnera_multiple_mutualisms
This is a repository for data and code pertaining to the manuscript "The genetic architecture of multiple mutualisms and mating system in Turnera ulmifolia"

Below are a description of each file included in this repository. For ease of interpretation, I have broken code and data files into conceptual sections.

######################################################################################

1. "Keystone" data summarizing key information or used in multiple analyses.

  a. Trait and seedset data.xlsx: an excel file containing phenotype data for all _Turnera ulmifolia_ plants. Second tab includes data from seedset experiment.

  b. pheno.txt: text file containing data for all plants. This file is uploaded in the following R codes : GMatrix_Estimation_final.R, Population_Mean_Trait_Variation.R.

  c. seedset.txt: text file containing data for the seedest experiment. This file is uploaded in the following R codes : GMatrix_Estimation_Final.R.

  d. Pedigree.xlsx: excel file containing pedigree information for all _Turnera ulmifolia_ plants. 

  e. pedmat.txt: text file containing pedigree information. This file is used in the following R codes : GMatrix_Estimation_final.R, Bayesian_Heritability_final.R
  
#######################################################################################
  
2. Fitting G matrices and univariate Bayesian objects

  a. GMatrix_Estimation_final.R: (sample) code that fits G matrices and univariate Bayesian objects for all populations of _T. ulmifolia_. Also includes code for   construction of G matrix for seedset experiment, the signigicance testing of diagonal elements, and the generation of randomized data sets for each population.
  
  b. Bayesian_Heritability_final.R: (sample) code that fits univariate Bayesian objects to randomized datasets for each population of _T. ulmifolia_ for the purposes of significance testing.
  
#####################################################################################
  
3. Comparison of G matrices

  a. Gcomp_final.R: code which assesses differences among our 5 population-specific estimates of G, and uploads randomly-generated G matrices for the purposes of significance testing.
  
  b. krzan.txt: text file containing the results of Krzanowski's analysis performed on null G matrices fitted per Morrisey et al., 2019, Evolution. This file is used in the following R codes : Gcomp_final.R
  
  c. fourth_order.txt: text file containing the results of the fourth order genetic covariance tensor performed on null G matrices fitted per Morrisey et al., 2019, Evolution. This file is used in the following R codes : Gcomp_final.R
  
  d. ranskew.txt: text file containing the results of the random skewer analysis performed on null G matrices fitted per Morrisey et al., 2019, Evolution. This file is used in the following R codes : Gcomp_final.R
  
  e-i. BT_ranG_final - SA_ranG_final: summary files containing the values for genetic variances and covariances fit based on randomized data, as per Morrisey et al., 2019, Evolution. This file is used in the following R codes : Gcomp_final.R, Gcomp_random.R.
  
  j. Dmatrix_estimation_Bayes.R: code used to estimate D and compare it to various estimates of G (Gw and population-specific G's)
  
  k. Gcom_random.R: code to be run on a server, that uploads the randomized G's and then runs comparisons for each set of null results. 
  
######################################################################################
  
4. Miscellaneous

  a. Population_Mean_Trait_Variation.R: code that generates Figure 3, which displays variation in mean trait values among populations.
