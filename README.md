# optimise-ar

This repository contains the code used to simulate results for the paper **A practical toolkit with recommendations for analysing and visualising patient-reported outcomes in early phase dose-finding oncology trials: OPTIMISE-AR**. 

## Background 
OPTIMISE-AR promotes appropriate analysis and data visualisation of PRO data, facilitating robust, patient-centred tolerability conclusions and supporting the broader development of tolerable and effective treatments. This repository provides R code for simulating synthetic data and generating the associated figures and case studies presented in this paper.

## Description of R files
* **for_paper** - this folder contains .R code and .pdf files required to generate each figure as per the OPTIMISE-AR toolkit. 

* **CheckMate066.csv** - this .csv file contains EORTC QLQ-C30 Global Health Status score data extracted from the CheckMate066 trial (NCT01721772)[1].
  
* **continuous_data.R** - this .R code generates dummy continuous PRO data. Functions required to run this code are defined in `functions_synthesise_data.R`.

* **dummy_continuous.csv** - this .csv file contains synthetic PRO data used to generate figures presented in this toolkit as per the file `continuous_data.R`.

* **dummy_ordinal.csv** - this .csv file contains synthetic PRO data used to generate figures presented in this toolkit as per the file `ordinal_data.R`.

* **functions_synthesise_data.R** - the .R code includes all functions required to simulate synthetic PRO data.

* **ordinal_data.R** - this .R code generates dummy ordinal PRO data. Functions required to run this code are defined in `functions_synthesise_data.R`.

[1] Long GV, Atkinson V, Ascierto PA, et al. Effect of nivolumab on health-related quality of life in patients with treatment-naïve advanced melanoma: Results from the phase iii checkmate 066 study. Annals of Oncology. 2016/10/01/ 2016;27(10):1940-1946. doi:https://doi.org/10.1093/annonc/mdw265
