# 🫁 COPD Severity Analysis — Wales (2015–2022)

> **What factors are associated with the risk of severe COPD or death in Wales?**  
> A data-driven vs hypothesis-driven approach using machine learning and survival analysis.

---

## 📌 Overview

This repository contains the full analytical pipeline for my MSc dissertation in **Health Data Science** (Swansea University, 2024), extended with additional machine learning techniques for portfolio demonstration.

The study investigates clinical and lifestyle risk factors associated with **severe COPD (Chronic Obstructive Pulmonary Disease) or death** in Wales, using data from the **SAIL (Secure Anonymised Information Linkage) Databank** — one of the largest anonymised health data repositories in the UK.

> ⚠️ **Note on data:** The original SAIL Databank data is confidential and cannot be shared publicly. This repository uses a **synthetic dataset** generated to mirror the statistical properties and distributions of the original cohort (n = 42,005). All code, methods, and findings are directly reproducible.

---

## 🔬 Research Question

*"What are the factors associated with the risk of severe COPD (or death) in Wales: a data-driven vs hypothesis-driven approach?"*

---

## 📊 Dataset

| Property | Value |
|---|---|
| **Source** | SAIL Databank, Swansea University (original) |
| **Study period** | January 2015 – March 2022 |
| **Cohort size** | ~42,000 severe COPD patients |
| **Data types** | GP records, hospital admissions, death records |

**Risk factors studied:** Smoking, obesity, passive smoking, respiratory diseases, age, gender, genetic predisposition (alpha-1 antitrypsin deficiency), toxic gas exposure

**Comorbidities studied:** Asthma, anxiety, diabetes, heart disease, hypertension, osteoporosis, tuberculosis

---

## 🧪 Methods

| Technique | Purpose |
|---|---|
| **Logistic Regression** | Predict severe COPD outcome; estimate odds ratios |
| **Decision Tree** | Visual classification of high-risk groups |
| **Random Forest** *(extended)* | Ensemble learning; variable importance ranking |
| **Cox Proportional Hazards** | Survival analysis; hazard ratios by factor |
| **Kaplan-Meier** *(extended)* | Survival curves by smoking, age, gender |

---

## 📁 Repository Structure

```
copd-wales-analysis/
│
├── data/
│   └── copd_simulated.csv          # Synthetic dataset (mirrors SAIL cohort)
│
├── scripts/
│   ├── 00_simulate_data.R          # Generate synthetic dataset
│   ├── 01_descriptive_stats.R      # EDA and prevalence charts
│   ├── 02_logistic_regression.R    # Logistic regression + ROC curves
│   ├── 03_cox_survival.R           # Cox model + Kaplan-Meier curves
│   └── 04_random_forest.R          # Random Forest + model comparison
│
├── reports/
│   └── *.png                       # All generated visualisations
│
└── README.md
```

---

## 🔑 Key Findings

### Risk Factors
- **Smoking** is the dominant risk factor — smokers are **>2× more likely** to experience severe outcomes (HR = 2.14, Cox model)
- **Age** is consistently the second most important factor across all models
- **Genetic predisposition** (alpha-1 antitrypsin deficiency) is significant in logistic regression but not in the Cox model — suggesting it influences disease development more than progression speed
- **Toxic gas exposure** showed no statistical significance — possibly due to low prevalence in the dataset

### Comorbidities
- **Diabetes** and **heart disease** are the most dangerous comorbidities (significant across both logistic and Cox models)
- **Asthma** paradoxically showed a *reduced* hazard ratio — likely because asthma patients are already actively managing their respiratory health
- **Hypertension** and **osteoporosis** showed borderline significance — requiring deeper investigation

### Public Health Insight
- A significant gap exists between COPD deaths and formal diagnoses — many patients die with COPD as a contributing cause without ever being diagnosed, pointing to a **critical underdiagnosis problem** in Wales

---

## 📦 Requirements

```r
# Install required R packages
install.packages(c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "survival",
  "survminer",
  "randomForest",
  "caret",
  "pROC",
  "broom",
  "tibble"
))
```

---

## ▶️ How to Run

```r
# Run scripts in order from the project root
source("scripts/00_simulate_data.R")      # Creates data/copd_simulated.csv
source("scripts/01_descriptive_stats.R")  # EDA plots
source("scripts/02_logistic_regression.R") # Logistic + ROC
source("scripts/03_cox_survival.R")        # Cox + Kaplan-Meier
source("scripts/04_random_forest.R")       # Random Forest + comparison
```

---

## 👤 Author

**Oluwatosin Samson Adewole**  
PgD Health Data Science — Swansea University (2024)  
BSc Microbiology — Wesley University Ondo (2013)  

🔗 [LinkedIn](https://www.linkedin.com/in/oluwatosin-adewole-) | 📧 adewoleoluwatosin3@gmail.com

---

## 📚 References

- Adewole, O.S. (2024). *Factors associated with severe COPD (or death) in Wales: a data-driven vs hypothesis-driven approach.* MSc Dissertation, Swansea University.
- SAIL Databank: [https://saildatabank.com](https://saildatabank.com)
- WHO COPD Fact Sheet: [https://www.who.int/news-room/fact-sheets/detail/chronic-obstructive-pulmonary-disease-(copd)](https://www.who.int/news-room/fact-sheets/detail/chronic-obstructive-pulmonary-disease-(copd))
- dos Santos, N.C. et al. (2022). Prevalence and Impact of Comorbidities in Individuals with COPD. *Tuberculosis and Respiratory Diseases, 85*(3), 205–220.

---

*This repository is for academic and portfolio purposes. The synthetic dataset is not a substitute for the original SAIL data.*
