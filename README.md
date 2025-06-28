# Gene-Expression-Analysis-Breast-Cancer
Analysis and ML classification for breast cancer subtypes using GDS1329 dataset.

## 📖 Overview

This project provides an in-depth gene expression analysis of the GDS1329 breast cancer dataset. The primary objectives are to identify differentially expressed genes among three main tumor subtypes (luminal, basal, and apocrine) and to build machine learning models to classify these subtypes based on their gene expression profiles. 

## 📊 Dataset

- Dataset Name: GDS1329: Molecular apocrine breast tumors 
- Source: NCBI Gene Expression Omnibus (GEO)
- Content: 49 human breast tumor samples (Homo sapiens)
- Platform: Affymetrix Human Genome U133A Array

## ⚙️ Project Workflow

The analysis in this project is divided into two main stages using R and Python:

**1. Differential Gene Expression Analysis (using R):**
- Downloaded and pre-processed data from the GEO database.
- Filtered low-variance genes to reduce noise, decreasing the gene count from 22,283 to 6,321.
- Identified Differentially Expressed Genes (DEGs) using the limma package.
- Visualized analysis results using Heatmaps, Volcano Plots, and Boxplots for interpretation. 

**2. Classification Modeling (using Python & Scikit-learn):**
- Prepared data for modeling: scaled features using
- StandardScaler and split the data into training and testing sets (50:50 ratio).
- Trained and evaluated 7 different classification models.
- Assessed each model's performance using 5-fold cross-validation, confusion matrices, and ROC-AUC curves. 

## Key Results & Model Performance
**Gene Expression Analysis Findings**

The analysis successfully identified hundreds of genes that were significantly differentially expressed among the three tumor subtypes. The heatmap visualization of the top 50 DEGs showed clear clustering patterns, where samples from the same subtype grouped together, confirming that gene expression profiles can effectively distinguish the tumor classes. 

![image](https://github.com/user-attachments/assets/a35effe8-bc4f-429f-9b5c-a373fb52b0b3)

**Classification Model Performance**

SVM (Linear), Logistic Regression, and Lasso Logistic Regression models demonstrated the best performance, achieving a cross-validation accuracy of 98%. 

The following table compares model performance based on the mean 5-fold cross-validation accuracy:


| Model | Accuracy (Mean CV) | Std Dev |
|---|---|---|
| **SVM (Linear)** | **98.00%** | **4.00%** |
| **Logistic Regression** | **98.00%** | **4.00%** |
| **Lasso Logistic Regression**| **98.00%** | **4.00%** |
| K-Nearest Neighbors (KNN) | 94.00% | 8.00% |
| Random Forest | 92.00% | 7.48% |
| Gaussian Naive Bayes | 89.78% | 6.34% |
| Decision Tree | 85.78% | 10.12% |

The SVM (Linear) model was considered the best overall due to its highly stable performance and its excellent suitability for high-dimensional data like gene expression profiles.



