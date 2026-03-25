# =============================================================================
# 04_random_forest.R
# Random Forest — Extending the Thesis with Ensemble Learning
# Author: Oluwatosin Samson Adewole
# Note: This extends the original dissertation (logistic regression + decision
#       tree) by adding Random Forest for improved predictive performance.
# =============================================================================

library(randomForest)
library(ggplot2)
library(dplyr)
library(caret)
library(pROC)

# --- Load data ---
copd <- read.csv("data/copd_simulated.csv")

cat("=== RANDOM FOREST ANALYSIS ===\n\n")

# --- Prepare features ---
features_risk <- c("age", "gender_sex", "genecode", "obese_code",
                   "poisongas_code", "respi_code", "smkcode")

features_comorbid <- c("age", "gender_sex", "anxiety_code", "asthma_code",
                        "diabetes_code", "heart_code", "hyper_code", "osteo_code")

copd$severe_tag <- as.factor(copd$severe_tag)

# Train/test split (80/20)
set.seed(123)
train_idx   <- createDataPartition(copd$severe_tag, p = 0.8, list = FALSE)
train_data  <- copd[train_idx, ]
test_data   <- copd[-train_idx, ]

cat(sprintf("Training set : %d patients\n", nrow(train_data)))
cat(sprintf("Test set     : %d patients\n\n", nrow(test_data)))

# ─────────────────────────────────────────────────────────────────────────────
# MODEL 1: Random Forest — Risk Factors
# ─────────────────────────────────────────────────────────────────────────────
cat("--- RF Model 1: Risk Factors ---\n")

rf_risk <- randomForest(
  x        = train_data[, features_risk],
  y        = train_data$severe_tag,
  ntree    = 500,
  mtry     = 3,
  importance = TRUE,
  seed     = 42
)

cat("\nRandom Forest (Risk Factors):\n")
print(rf_risk)

# Predictions
pred_risk_prob <- predict(rf_risk, test_data[, features_risk], type = "prob")[, 2]
pred_risk_class <- predict(rf_risk, test_data[, features_risk])

# Confusion matrix
cm_risk <- confusionMatrix(pred_risk_class, test_data$severe_tag, positive = "1")
cat("\nConfusion Matrix — Risk Factors:\n")
print(cm_risk$table)
cat(sprintf("\nAccuracy  : %.3f\n", cm_risk$overall["Accuracy"]))
cat(sprintf("Sensitivity: %.3f\n", cm_risk$byClass["Sensitivity"]))
cat(sprintf("Specificity: %.3f\n", cm_risk$byClass["Specificity"]))

# AUC
roc_rf_risk <- roc(as.numeric(test_data$severe_tag) - 1,
                    pred_risk_prob, quiet = TRUE)
cat(sprintf("AUC        : %.3f\n", auc(roc_rf_risk)))

# Variable importance plot
imp_risk <- as.data.frame(importance(rf_risk)) %>%
  tibble::rownames_to_column("Variable") %>%
  arrange(desc(MeanDecreaseGini)) %>%
  mutate(Variable = recode(Variable,
    age            = "Age",
    gender_sex     = "Gender",
    genecode       = "Genetic (Alpha-1 AT)",
    obese_code     = "Obesity",
    poisongas_code = "Toxic Gas Exposure",
    respi_code     = "Respiratory Disease",
    smkcode        = "Smoking"
  ))

p_imp_risk <- ggplot(imp_risk,
                      aes(x = MeanDecreaseGini,
                          y = reorder(Variable, MeanDecreaseGini))) +
  geom_col(fill = "#C0392B", alpha = 0.85, width = 0.6) +
  geom_text(aes(label = round(MeanDecreaseGini, 1)),
            hjust = -0.1, size = 3.5, fontface = "bold") +
  labs(
    title    = "Random Forest: Variable Importance — Risk Factors",
    subtitle = sprintf("500 trees | AUC = %.3f | Simulated SAIL Databank Wales cohort",
                       auc(roc_rf_risk)),
    x        = "Mean Decrease Gini",
    y        = NULL,
    caption  = "Source: Adewole (2024) — Swansea University MSc Dissertation (Extended)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey40")
  )

ggsave("reports/12_rf_importance_risk.png", p_imp_risk,
       width = 9, height = 5, dpi = 180)
cat("✅ Saved: reports/12_rf_importance_risk.png\n")

# ─────────────────────────────────────────────────────────────────────────────
# MODEL 2: Random Forest — Comorbidities
# ─────────────────────────────────────────────────────────────────────────────
cat("\n--- RF Model 2: Comorbidities ---\n")

rf_comorbid <- randomForest(
  x          = train_data[, features_comorbid],
  y          = train_data$severe_tag,
  ntree      = 500,
  mtry       = 3,
  importance = TRUE,
  seed       = 42
)

pred_comorbid_prob  <- predict(rf_comorbid, test_data[, features_comorbid], type = "prob")[, 2]
pred_comorbid_class <- predict(rf_comorbid, test_data[, features_comorbid])

cm_comorbid <- confusionMatrix(pred_comorbid_class, test_data$severe_tag, positive = "1")
roc_rf_comorbid <- roc(as.numeric(test_data$severe_tag) - 1,
                        pred_comorbid_prob, quiet = TRUE)

cat(sprintf("Accuracy  : %.3f\n", cm_comorbid$overall["Accuracy"]))
cat(sprintf("AUC        : %.3f\n", auc(roc_rf_comorbid)))

imp_comorbid <- as.data.frame(importance(rf_comorbid)) %>%
  tibble::rownames_to_column("Variable") %>%
  arrange(desc(MeanDecreaseGini)) %>%
  mutate(Variable = recode(Variable,
    age           = "Age",
    gender_sex    = "Gender",
    anxiety_code  = "Anxiety",
    asthma_code   = "Asthma",
    diabetes_code = "Diabetes",
    heart_code    = "Heart Disease",
    hyper_code    = "Hypertension",
    osteo_code    = "Osteoporosis"
  ))

p_imp_comorbid <- ggplot(imp_comorbid,
                          aes(x = MeanDecreaseGini,
                              y = reorder(Variable, MeanDecreaseGini))) +
  geom_col(fill = "#1A5276", alpha = 0.85, width = 0.6) +
  geom_text(aes(label = round(MeanDecreaseGini, 1)),
            hjust = -0.1, size = 3.5, fontface = "bold") +
  labs(
    title    = "Random Forest: Variable Importance — Comorbidities",
    subtitle = sprintf("500 trees | AUC = %.3f | Simulated SAIL Databank Wales cohort",
                       auc(roc_rf_comorbid)),
    x        = "Mean Decrease Gini",
    y        = NULL,
    caption  = "Source: Adewole (2024) — Swansea University MSc Dissertation (Extended)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey40")
  )

ggsave("reports/13_rf_importance_comorbid.png", p_imp_comorbid,
       width = 9, height = 5, dpi = 180)
cat("✅ Saved: reports/13_rf_importance_comorbid.png\n")

# ─────────────────────────────────────────────────────────────────────────────
# MODEL COMPARISON: Logistic Regression vs Random Forest
# ─────────────────────────────────────────────────────────────────────────────
cat("\n--- Model Comparison (Risk Factors) ---\n")

# Logistic regression on same train/test split for fair comparison
logit_risk <- glm(
  severe_tag ~ age + gender_sex + genecode + obese_code +
               poisongas_code + respi_code + smkcode,
  data   = train_data %>% mutate(severe_tag = as.numeric(severe_tag) - 1),
  family = binomial
)
pred_logit_prob <- predict(logit_risk,
                            test_data[, features_risk],
                            type = "response")
roc_logit <- roc(as.numeric(test_data$severe_tag) - 1,
                  pred_logit_prob, quiet = TRUE)

comparison <- data.frame(
  Model = c("Logistic Regression", "Random Forest"),
  AUC   = c(round(auc(roc_logit), 3), round(auc(roc_rf_risk), 3))
)

cat("\nModel AUC Comparison (Risk Factors):\n")
print(comparison)

p_compare <- ggplot(comparison, aes(x = Model, y = AUC, fill = Model)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  geom_text(aes(label = sprintf("AUC = %.3f", AUC)),
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Logistic Regression" = "#85929E",
                                "Random Forest"       = "#C0392B")) +
  ylim(0, 1) +
  labs(
    title    = "Model Comparison: Logistic Regression vs Random Forest",
    subtitle = "Predicting Severe COPD — Risk Factors Model",
    x        = NULL,
    y        = "AUC (Area Under ROC Curve)",
    caption  = "Source: Adewole (2024) — Swansea University MSc Dissertation (Extended)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("reports/14_model_comparison_auc.png", p_compare,
       width = 7, height = 5, dpi = 180)
cat("✅ Saved: reports/14_model_comparison_auc.png\n")

cat("\n=== Random Forest analysis complete ===\n")
