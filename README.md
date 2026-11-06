# Machine Learning Classification of Pick's Disease Using Single-Cell RNA Sequencing Data
## Project Overview

This project investigates the application of statistical learning and machine learning techniques to distinguish Pick's disease (a subtype of frontotemporal lobar degeneration) from neurologically healthy controls using single-cell RNA sequencing (scRNA-seq) data.

The primary objective is to evaluate whether gene expression patterns can provide meaningful predictive information about disease status while identifying genes that may contribute to the classification of Pick's disease.

By combining high-dimensional biological data with statistical modeling, this project explores the potential applications of machine learning in neurodegenerative disease research and biomarker discovery.

## Research Objectives

- Investigate differences in gene expression between Pick's disease and healthy control samples.
- Explore dimensionality reduction techniques to identify patterns in high-dimensional gene expression data.
- Develop predictive classification models using gene expression features.
- Compare model performance using cross-validation and classification metrics.
- Investigate genes selected by penalized regression models as potential disease-associated features.

## Dataset

The analysis utilizes single-cell RNA sequencing data containing gene expression measurements and associated sample metadata.

**Dataset characteristics:**
- Original dataset: 432,555 cells × 29,363 genes.
- Disease-specific subset: 196,837 cells × 29,363 genes.
- Observational units: Individual cells, with donor-level metadata.
- Outcome variable: Disease status (Pick's disease vs. healthy control).
- Predictors: Gene expression measurements.

To account for the hierarchical structure of the data, pseudobulk aggregation was explored by donor and cell type.

*Note: The analysis distinguishes between individual cells and biological donors. Donor-level aggregation is important when evaluating disease classification to avoid treating cells from the same individual as independent biological samples.*

## Methodology
### 1. Data Preprocessing
- Subsetting samples to include Pick's disease and healthy control observations.
- Processing gene expression matrices and associated metadata.
- Exploring cell-type-specific gene expression patterns.
- Aggregating expression measurements through pseudobulk analysis.

### 2. Exploratory Data Analysis
- Principal Component Analysis (PCA).
- Examination of gene expression distributions.
- Exploration of disease-associated patterns across biological samples.
- Investigation of potential differences between disease groups.

### 3. Statistical Modeling and Machine Learning
Several classification approaches were explored:

| Model | Purpose |
|---|---|
| Logistic Regression | Baseline disease classification |
| LASSO Regression | Feature selection and regularization |
| Support Vector Machine | Nonlinear classification |

### 4. Model Validation

Leave-One-Out Cross-Validation (LOOCV) was used to evaluate model performance across individual observations.

Additional model evaluation included:
- Receiver Operating Characteristic (ROC) curves.
- Area Under the ROC Curve (AUC).
- Deviance and Akaike Information Criterion (AIC).
- Comparison of regularized and unpenalized logistic regression.
- Examination of selected features across model iterations.

Special attention is given to the limitations of small donor sample sizes and the importance of avoiding data leakage between training and validation observations.

## Preliminary Findings

The preliminary analysis identified several genes of interest through different feature-selection and classification approaches.

Genes investigated include:

- **RELN**
- **CXCR4**
- **SERPINH1**
- **KIT**
- **GBP2**
- **HBA1 / HBA2**
- **HBB**
- **HSPB1**
- **TPD52L1**

Additional genes selected through LASSO regularization included QKI, PLP1, SLC44A1, SYT1, and EDIL3.

These findings represent exploratory statistical associations and predictive features. They should not be interpreted as validated diagnostic biomarkers or evidence of causal relationships.

## Technologies & Tools

- **R:** Statistical analysis and model development.
- **Python:** Data processing and machine learning.
- **scRNA-seq:** Single-cell gene expression analysis.
- **PCA:** Dimensionality reduction.
- **glmnet:** Penalized regression and LASSO.
- **pROC:** ROC curves and AUC evaluation.
- **edgeR:** Gene expression analysis.
- **Machine Learning:** Logistic regression, SVM, Random Forest, KNN, LDA, and QDA.

## Project Significance

This project demonstrates the application of biostatistical methods and machine learning to complex biomedical datasets.

By investigating gene expression patterns associated with Pick's disease, the analysis explores how computational methods can contribute to neurodegenerative disease research, predictive modeling, and the identification of potential molecular features for future investigation.

The project also highlights important statistical considerations when working with high-dimensional biological data, including feature selection, model validation, and the limitations of small sample sizes.

## Limitations and Future Directions

Potential future extensions include:

- Expanding the number of biological donors.
- Investigating cell-type-specific disease signatures.
- Conducting differential gene expression analysis.
- Comparing additional feature-selection methods.
- Evaluating models using donor-level nested cross-validation.
- Investigating external validation datasets.
- Exploring explainable machine learning techniques.

**Disclaimer:** This is an exploratory academic research project. Model predictions and identified genes are not intended for clinical diagnosis or treatment decisions.
