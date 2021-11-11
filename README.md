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

- **TCGA_CEP89_status_vs_CA20.R:** shows the association between CEP89 amplification status and CA20 scores (refelctive of centromere amplification) using data obtained from de Almeida et al 2019 (PLOS Computational Biology)
- **TCGA_CoxProportionalHazardsModelCEP89_stratification.R:** calculates Cox proportional hazards analysis estimates for the impact of CEP89 expression on patient prognosis within individual cancer studies.
- **TCGA_QuiescenceScore_vs_CA20.R:** shows the pan-cancer relationship between CA20 scores (refelctive of centromere amplification) and tumour quiescence across the TCGA cohort, using data obtained from de Almeida et al 2019 (PLOS Computational Biology)

## PanCancer_ElasticNetRegression

This folder contains the code for running the pan-cancer elastic net regression model designed to identify genomic features associated with cellular dormancy in cancer.

- **FishersExactTest_DDRmutation_vs_quiescence.R:** shows the depletion of mutations across the DNA damage repair pathways in TCGA samples with high levels of quiescence
- **MCF7_PTEN_mutation_vs_DoublingTime.R:** shows the association between quiescence scores, doubling time and PTEN mutational status in MCF7 cell line stains using data from Ben-David et al 2018, Nature
- **Quiescence_MSI_association.R:** shows the association between quiescence and microsatelite instability using data from Cortes-Ciriano et al 2017, Nature Communications
- **PanCancer_ElasticNetRegression.R:** runs 1000 iterations of the pan-cancer elastic net regression model to identuft genomic events significantly associated with dormancy.
- **PanCancerModel_InternalValidation_and_SHAP.R:** shows the correlation between predicted and observed quiescence scores in the interal TCGA validation cohort. It also illustrates the importance of individual features using a SHAP plots.

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
- **QuiescenceScoreValidation_ROC.R:** illustrates the performance of the combined z-score quiescence scoring methodology on separating actively proliferating and quiescent cells
- **QS_vs_ProliferationMarkers_Performance.R:** Compares the performance of the combined z-score quiescence scoring methodology on separating actvely proliferating and quiescent cells, compared to other markers of proliferating cells.
- **QuiescenceScore_vs_QuiescenceDepth.R:** Shows quiescence levels estimates of rat embryonic fibroblast cells under serum starvation for various amounts of time, using data from Fujimaki et al 2019, Proc Natl Acad Sci USA
- **CDK4_6_Inhibition_QS_Validation.R, Contact_Inhibition_QS_Validation.R, MEK_Inhibition_QS_Validation.R, SerumStarvation_QS_Validation.R, Spontaneous_QS_Validation.R"** show the correlation between quiescence score estimates for the specific quiescence programmes and the expression of genes associated with the corresponding form of quiescence in the literature.

## QuiescenceTypeClassification

This folder contains code used to estimate the dominant form of quiescence across different TCGA cohort cancer types and in scRNA-seq datasets.

- **TCGA_QuiescenceTypeClassification.R:** identifies the the main quiescence programme within highly dormant TCGA samples.
- **GSE134839_knn.R:** performs knn analysis to classify dormant NSCLC cells into individual quiescence type categories

## TCGA_DataDownload_and_Processing

This folder contains code used to download and process TCGA expression, copy-number variation and mutational data.

## TCGA_DormancyEvaluation

This folder contains code used to profile quiescence levels of solid cancer, primary tumour samples, across the TCGA cohort. 

- **TCGA_QS_vs_SenescenceMarkers.R:** investigates the relationship between quiescence and senescnece in cancer
- **TCGA_QS_vs_StemnessExpressionIndicies.R:** investigates the relationship between quiescence and stem cell marker expression in cancer using data Malta et al 2018, Cell
- **TCGA_QS_vs_TelomeraseActivity.R:** investigates the relationship between quiescence and telomerase activation in cancer using data from Noureen et al 2021, Nature Communications
- **TCGA_CDKN1A_expression_QS.R**: shows the association between quiescence levels and CDKN1A (encoding p21) expression in cancer
- **TCGA_Kmeans_QuiescenceGenes.R:** clusters TCGA samples into high quiescence and low quiescence groups, after the removal of tissue specific expression patterns
- **TCGA_PHATE_with_TissueTypeComBatCorrection.R:"** performs PHATE dimentionality reduction on TCGA samples based on the expression of quiescence biomarker genes, following the removal of tissue specific expression patterns.
- **TCGA_QS_vs_WGD.R** shows the association between quiescence levels of whole-genome duplication events in cancer
- **TCGA_QS_correlation_with_ProliferationMarkers.R:** shows the association between quiescence scores and the expression of previously reported markers of proliferating cells 
- **TCGA_PHATE_without_TissueTypeComBatCorrection.R:** performs PHATE dimentionality reduction on TCGA samples based on the expression of quiescence biomarker genes, without the removal of tissue specific expression patterns.
- **AvStemCellDivision_vs_meanTCGA_QS.R:** shows the correlation between mean quiescence scores across TCGA cancer types and reported stem cell division estimates for the corresponding tissues, using data from Tomasetti and Vogelstein, 2015, Science
