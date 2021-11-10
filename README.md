# CancerDormancy

This repository contains code for the evaluation of cellular dormancy/quiescence in cancer.

![CancerDormancy_picture](https://user-images.githubusercontent.com/51481454/141082122-c3711ca2-1c96-4853-bb46-589356403996.png)

# Table of contents

## BRCA_Quiescence

This folder contains the code which explores the association between breast cancer dormancy and genomic features shown to be linked with this state by a breast cancer specific random forest model.

## CEP89_QuiescenceAssociation

This folder cotains code which explores the impact of CEP89 activity on quiescence-proliferation decisions.

## PanCancer_ElasticNetRegression

This folder contains the code for running the pan-cancer elastic net regression model designed to identify genomic features associated with cellular dormancy in cancer.

## Prognosis_and_treatment_response_analysis

This folder contains the code which explores the impact of cellular dormancy on patient prognosis as well as treatment response in scRNA-seq datasets.

## QuiescenceBiomarkerGeneIdentification

This folder contains code which identifies a list of generic quiescence biomarker genes, as well as genes specificially descriptive of 5 individual quiescence types, based on data published by Min and Spencer, 2019, PLoS Biol.

## QuiescenceScoreValidation

This folder contains code used to validate the combined z-score quiescence scoring metholody used to asssess the presence of a generic quiescence state as well as 5 specific quiescence types.

## QuiescenceTypeClassification

This folder contains code used to estimate the dominant form of quiescence across different TCGA cohort cancer types and in scRNA-seq datasets.

## TCGA_DataDownload_and_Processing

This folder contains code used to download and process TCGA expression, copy-number variation and mutational data.

## TCGA_DormancyEvaluation

This folder contains code used to profile quiescence levels of solid cancer, primary tumour samples, across the TCGA cohort. 
