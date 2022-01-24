# This folder contains source code used by the repository.

## R scripts 

+ [`parse_offsides_by_all_ae.R`](parse_offsides_by_all_ae.R) processes OFFSIDES dataset and generates positive/negative samples for all adverse events. Adverse events with fewer than 100 positive/negative samples were removed from the parsed result file. Positive samples are drugs with lower bound of Proportional Reporting Ratio (PRR) > 1. Negative samples are drugs with lower bound of Proportional Reporting Ratio (PRR) < 1, and upper bound of Proportional Reporting Ratio (PRR) > 1.

+ [`parse_tox21_sample.R`](parse_tox21_sample.R) parses Tox21 screening dataset and generates active/inactive compounds for all assays from Tox21.

+ [`parse_lincs.R`](parse_lincs.R) processes LINCS gene expression data after perturbation to obtain Tox21 compound-induced gene expression profile that matches input query time and dosage.

+ [`parse_tox21_atc.R`](parse_tox21_atc.R) maps compounds from Tox21 HepG2 viability assay to their ATC drug categories.

+ [`parse_dili_phenotype.R`](parse_dili_phenotype.R) identifies positive/negative DSSTox compounds for drug-induced liver injury (DILI) phenotypes, and tests the association between each DILI phenotype and HepG2 cell viability.

+ [`parse_diki_phenotype.R`](parse_diki_phenotype.R) identifies positive/negative DSSTox compounds for drug-induced kidney injury (DIKI) phenotypes, and tests the association between each DIKI phenotype and HEK293 cell viability.

## Executable shell scripts

+ [`run/parse_tox21_sample.sh`](run/parse_tox21_sample.sh) runs [`parse_tox21_sample.R`](parse_tox21_sample.R) on screening data of Tox21 assays.

+ [`run/parse_lincs.sh`](run/parse_lincs.sh) implements [`parse_lincs.R`](parse_lincs.R) to generate Tox21 compound-induced gene expression profile that matches 6 input time-dosage combinations.
