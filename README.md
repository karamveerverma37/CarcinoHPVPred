# CarcinoHPVPred

- Standalone version of CarcinoHPVPred is consist of scripts and models written and tested with perl-5 and python-3.9.9.
- It works with linux OS only.
- This tool can be used to predict the carcinogenic nature of new or already available HPV strains. This tool uses genomic information from two HPV core genes E2 and E6.

## Dependencies:

- perl-5
- python-3.9.9
- sklearn: 1.0.2
- numpy: 1.23.1
- scipy: 1.9.0 
- pandas: 1.4.3
- matplotlib: 3.5.2
- joblib: 1.1.0

## Installation:

`unzip CarcinoHPVPred.zip`

To run program simply type:

`cd CarcinoHPVPred/`

`perl CarcinoHPVPred.pl abc.fasta modelName type_of_input`

abc.fasta => genome or predicted genes in fasta format

modelName => logistic_regression or svm or knn or lda

type_of_input => full_genome or E2_E6 or E6
