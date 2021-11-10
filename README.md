# CancerDormancy

This repository contains code for the evaluation of cellular dormancy/quiescence in cancer.

![CancerDormancy_picture](https://user-images.githubusercontent.com/51481454/141082122-c3711ca2-1c96-4853-bb46-589356403996.png)

# Table of contents

## BRCA_Quiescence

This folder contains the code which explores the association between breast cancer dormancy and genomic features shown to be linked with this state by a breast cancer specific random forest model.

- **BRCA_quiescence_apobec_sbs2_association_heatmap.R:** shows the association between dormancy scores, APOBEC enrichment and breast cancer subtypes in the TCGA cohort 
- **BRCA_subtype_SBS2_prevelance.R:** compares levels of APOBEC-linked SBS2 mutational burden across breast cancer subtypes within the TCGA cohort.

## CEP89_QuiescenceAssociation

This folder cotains code which explores the impact of CEP89 activity on quiescence-proliferation decisions.

- **TCGA_CEP89_status_vs_CA20.R:** shows the association between CEP89 amplification status and CA20 scores (refelctive of centromere amplification).
- **TCGA_CoxProportionalHazardsModelCEP89_stratification.R:** calculates Cox proportional hazards analysis estimates for the impact of CEP89 expression on patient prognosis within individual cancer studies.
- **TCGA_QuiescenceScore_vs_CA20.R:** shows the pan-cancer relationship between CA20 scores (refelctive of centromere amplification) and tumour quiescence across the TCGA cohort.

## PanCancer_ElasticNetRegression

This folder contains the code for running the pan-cancer elastic net regression model designed to identify genomic features associated with cellular dormancy in cancer.

## Prognosis_and_treatment_response_analysis

This folder contains the code which explores the impact of cellular dormancy on patient prognosis as well as treatment response in scRNA-seq datasets.

- **BulkRNAseqDatasets_PalbociclibTreatmentResponse.R:** performs the comparison of quiescene programme estimates in cancer cell lines before and after Palbociclib treatment 
- **BulkRNAseqDatasets_TreatmentResponse.R:** performs the comparison of quiescence levels in bulkRNA-seq cancer cell line datasets in response to Palbociclib, 5-FU and BRAF inhibition.
- **GSE134839_treatment_response.R:** shows changes in quiescence estimates within PC9 NSCLC cells in response to EGFR inhibition.
- **GSE149224_treatment_response.R:** shows changes in quiescence estimates within TP53 WT and TP53 MT colorectal cancer cells in response to 5-FU treatment.
- **TCGA_CoxProportionalHazardsModel_CommonQS.R:** runs a Cox proportional hazards model to determine the effect of quiescence on disease specific survival while accounting for cancer type, tumour stage, mutational rate, as well as patient age and sex.
- **TCGA_CoxProportionalHazardsModel_IndividualCancerTypes.R:** runs a Cox proportional hazards model to determine the effect of quiescence on disease specific survival within specific cancer types while accounting for tumour stage.
- **TCGA_CoxProportionalHazardsModel_IndividualQuiescenceTypes.R:** runs a Cox proportional hazards model to determine the effect of different quiescence types on disease specific survival while accounting for cancer type, tumour stage, mutational rate, as well as patient age and sex.

## QuiescenceBiomarkerGeneIdentification

This folder contains code **(QuiescenceBiomarker_Identification.R)** which identifies a list of generic quiescence biomarker genes, as well as genes specificially descriptive of 5 individual quiescence types, based on data published by Min and Spencer, 2019, PLoS Biol.

## QuiescenceScoreValidation

This folder contains code used to validate the combined z-score quiescence scoring metholody used to asssess the presence of a generic quiescence state as well as 5 specific quiescence types.

## QuiescenceTypeClassification

This folder contains code used to estimate the dominant form of quiescence across different TCGA cohort cancer types and in scRNA-seq datasets.

## TCGA_DataDownload_and_Processing

This folder contains code used to download and process TCGA expression, copy-number variation and mutational data.

## TCGA_DormancyEvaluation

This folder contains code used to profile quiescence levels of solid cancer, primary tumour samples, across the TCGA cohort. 
