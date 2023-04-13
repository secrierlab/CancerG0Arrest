# Genomic hallmarks and therapeutic implications of G0 cell cycle arrest in cancer

## Authors: Anna Wiecek, Maria Secrier

This repository contains code for the evaluation of cell cycle arrest (G0) in cancer based on transcriptomics data.

![CancerQuiescence-workflow](https://user-images.githubusercontent.com/1819561/155572705-cd6097d9-eb47-404e-9bc6-c927cfec6090.png)


## System Requirements
Operating system(s): Unix (linux, mac)

Programming Language: R

## Data
The majority of the data used in this work can be found in the folder /Data.

## How to run a script (DEMO)

1) Check that all the depencies are installed
2) Open RStudio or an R-shell
3) Set the working directory to the path with the script of interest as follows

For example if the user wants to run the pan-cancer elastic net regression model (PanCancer_ElasticNetRegression.R) then these are the steps to follow:

```r
# Move to the respective directory 
setwd("CancerCellQuiescence/PanCancer_ElasticNetRegression/ElasticNetRegressionModel/Code/")

# run the script
source("PanCancer_ElasticNetRegression.R")
```

# Table of contents

## BRCA_Quiescence

This folder contains the code which explores the association between breast cancer cell quiescence and genomic features shown to be linked with this state by a breast cancer specific random forest model.

- **BRCA_quiescence_apobec_sbs2_association_heatmap.R:** shows the association between quiescence scores, APOBEC enrichment and breast cancer subtypes in the TCGA cohort 
- **BRCA_subtype_SBS2_prevelance.R:** compares levels of APOBEC-linked SBS2 mutational burden across breast cancer subtypes within the TCGA cohort.
- **random_forest-BRCA-code-for-paper.ipynb:** creates the breast cancer specifc random forest model for quiescence prediction

## CEP89_QuiescenceAssociation

This folder cotains code which explores the impact of CEP89 activity on quiescence-proliferation decisions.

- **TCGA_CEP89_status_vs_CA20.R:** shows the association between CEP89 amplification status and CA20 scores (reflective of centromere amplification) using data obtained from de Almeida et al 2019 (PLOS Computational Biology)
- **TCGA_CoxProportionalHazardsModelCEP89_stratification.R:** calculates Cox proportional hazards analysis estimates for the impact of CEP89 expression on patient prognosis within individual cancer studies.
- **TCGA_QuiescenceScore_vs_CA20.R:** shows the pan-cancer relationship between CA20 scores (reflective of centromere amplification) and tumour quiescence across the TCGA cohort, using data obtained from de Almeida et al 2019 (PLOS Computational Biology)

## PanCancer_ElasticNetRegression

This folder contains the code for running the pan-cancer elastic net regression model designed to identify genomic features associated with cellular quiescence in cancer.

- **FishersExactTest_DDRmutation_vs_quiescence.R:** shows the depletion of mutations across the DNA damage repair pathways in TCGA samples with high levels of quiescence
- **MCF7_PTEN_mutation_vs_DoublingTime.R:** shows the association between quiescence scores, doubling time and PTEN mutational status in MCF7 cell line stains using data from Ben-David et al 2018, Nature
- **Quiescence_MSI_association.R:** shows the association between quiescence and microsatelite instability using data from Cortes-Ciriano et al 2017, Nature Communications
- **PanCancer_ElasticNetRegression.R:** runs 1000 iterations of the pan-cancer elastic net regression model to identuft genomic events significantly associated with quiescence.
- **PanCancerModel_InternalValidation_and_SHAP.R:** shows the correlation between predicted and observed quiescence scores in the interal TCGA validation cohort. It also illustrates the importance of individual features using a SHAP plots.
- **ValidationBallonPlot.R:** shows the association between genomic features identifed by the elastic net model and quiescence across the TCGA cohort and cBioPortal studies

## Prognosis_and_treatment_response_analysis

This folder contains the code which explores the impact of cellular quiescence on patient prognosis as well as treatment response in bulk RNA and scRNA-seq datasets.

- **BulkRNAseqDatasets_PalbociclibTreatmentResponse.R:** performs the comparison of quiescene programme estimates in cancer cell lines before and after Palbociclib treatment 
- **BulkRNAseqDatasets_TreatmentResponse.R:** performs the comparison of quiescence levels in bulkRNA-seq cancer cell line datasets in response to Palbociclib, 5-FU and BRAF inhibition.
- **GSE134839_treatment_response.R:** shows changes in quiescence estimates within PC9 NSCLC cells in response to EGFR inhibition.
- **GSE134499_treatment_response.R:"** shows quiescence dynamics upon various treatment modalities in MCF7 and A549 cell lines.
- **GSE149224_treatment_response.R:** shows changes in quiescence estimates within TP53 WT and TP53 MT colorectal cancer cells in response to 5-FU treatment.
- **GSE134838_treatment_response.R:** shows changes in quiescence estimates within melanoma cells in response to BRAF inhibition
- **GSE134836_treatment_response.R:** shows changes in quiescence estimates within NSCLC cells in response to Tyrosine kinase inhibition.
- **GSE137912_treatment_response.R:** shows changes in quiescence estimates within lung adenocarcinoma cells in response to KRAS inhibition
- **TCGA_CoxProportionalHazardsModel_CommonQS.R:** runs a Cox proportional hazards model to determine the effect of quiescence on disease specific survival while accounting for cancer type, tumour stage, mutational rate, as well as patient age and sex.
- **TCGA_CoxProportionalHazardsModel_IndividualCancerTypes.R:** runs a Cox proportional hazards model to determine the effect of quiescence on disease specific survival within specific cancer types while accounting for tumour stage.
- **TCGA_CoxProportionalHazardsModel_IndividualQuiescenceTypes.R:** runs a Cox proportional hazards model to determine the effect of different quiescence types on disease specific survival while accounting for cancer type, tumour stage, mutational rate, as well as patient age and sex.
- **TCGA_QS_vs_TumourStage.R:** shows the relationship between quiescence scores and tumour stage across the TCGA cohort

## QuiescenceBiomarkerGeneIdentification

This folder contains code **(QuiescenceBiomarker_Identification.R)** which identifies a list of generic quiescence biomarker genes, as well as genes specificially descriptive of 5 individual quiescence types, based on data published by Min and Spencer, 2019, PLoS Biol.

## QuiescenceScoreValidation

This folder contains code used to validate the combined z-score quiescence scoring metholody used to asssess the presence of a generic quiescence state as well as 5 specific quiescence types.
- **QuiescenceScoreValidation_ROC.R:** illustrates the performance of the combined z-score quiescence scoring methodology on separating actively proliferating and quiescent cells
- **QS_vs_ProliferationMarkers_Performance.R:** Compares the performance of the combined z-score quiescence scoring methodology on separating actvely proliferating and quiescent cells, compared to other markers of proliferating cells.
- **QuiescenceScore_vs_QuiescenceDepth.R:** Shows quiescence levels estimates of rat embryonic fibroblast cells under serum starvation for various amounts of time, using data from Fujimaki et al 2019, Proc Natl Acad Sci USA
- **CDK4_6_Inhibition_QS_Validation.R, Contact_Inhibition_QS_Validation.R, MEK_Inhibition_QS_Validation.R, SerumStarvation_QS_Validation.R, Spontaneous_QS_Validation.R"** show the correlation between quiescence score estimates for the specific quiescence programmes and the expression of genes associated with the corresponding form of quiescence in the literature.
- **NSCLC_CellLine_QS_EuD_pRB_correlation.R:** shows the correlation between quiescence score estimates and Edu + pRB measurements across NSCLC cell lines.

## QuiescenceTypeClassification

This folder contains code used to estimate the dominant form of quiescence across different TCGA cohort cancer types and in scRNA-seq datasets.

- **TCGA_QuiescenceTypeClassification.R:** identifies the the main quiescence programme within highly dormant TCGA samples.
- **GSE134839_knn.R:** performs knn analysis to classify dormant NSCLC cells into individual quiescence type categories

## TCGA_DataDownload_and_Processing

This folder contains code used to download and process TCGA expression, copy-number variation and mutational data.

- **TCGA_combining_CNV_data.R:** combines CNV data from individual TCGA studies.
- **TCGA_Combining_FPKM_ExpressionData.R:** combines FPKM normalised expression data from individual TCGA studies.
- **TCGA_Expression_CNV_DataDownload.R:** downloads expression and CNV TCGA data for all solid tissue cancer samples from TCGABiolinks.
- **TCGA_ExpressionData_PurityScaling.R:** scales TCGA FPKM normalised expression data accoridng to tumour sample purity.

## TCGA_QuiescenceEvaluation

This folder contains code used to profile quiescence levels of solid cancer, primary tumour samples, across the TCGA cohort. 
- **TCGA_QuiescenceScore_Calculation.R:** calculates quiescence scores across the TCGA cohort
- **TCGA_QS_TissueVariation.R:** shows variation in quiescence score estimates across different cancer type studies within the TCGA cohort
- **TCGA_QS_vs_DREAM_MutationalStatus.R:** shows the association between quiescence levels and DREAM complex mutational status in cancer
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

# How to cite
At present, a version of this manuscript is available on bioRxiv: https://www.biorxiv.org/content/10.1101/2021.11.12.468410v3

# Copyright
This code is free and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY. See the GNU General Public License for more details.
