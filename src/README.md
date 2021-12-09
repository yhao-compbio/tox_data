# This folder contains source code used by the repository.

## R scripts 

+ [`parse_offsides_by_all_ae.R`](parse_offsides_by_all_ae.R) processes OFFSIDES dataset and generates positive/negative samples for all adverse events. Adverse events with fewer than 100 positive/negative samples were removed from the parsed result file. Positive samples are drugs with lower bound of Proportional Reporting Ratio (PRR) > 1. Negative samples are drugs with lower bound of Proportional Reporting Ratio (PRR) < 1, and upper bound of Proportional Reporting Ratio (PRR) > 1.

+ [`parse_tox21_sample.R`](parse_tox21_sample.R) parses Tox21 screening dataset and generates active/inactive compounds for all assays from Tox21.

## Executable shell scripts

+ [`run/parse_tox21_sample.sh`](run/parse_tox21_sample.sh) runs [`parse_tox21_sample.R`](parse_tox21_sample.R) on screening data of Tox21 assays.
