# =============================================================================
# 02_logistic_regression.R
# Logistic Regression — Risk Factors & Comorbidities for Severe COPD
# Author: Oluwatosin Samson Adewole
# Mirrors methodology from MSc dissertation (Swansea University, 2024)
# =============================================================================

library(ggplot2)
library(dplyr)
library(broom)
library(pROC)

# --- Load data ---
copd <- read.csv("data/copd_simulated.csv")

cat("=== LOGISTIC REGRESSION ANALYSIS ===\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# MODEL 1: Risk Factors
# ─────────────────────────────────────────────────────────────────────────────
cat("--- Model 1: Risk Factors ---\n")

model_risk <- glm(
  severe_tag ~ age + gender_sex + genecode + obese_code +
               poisongas_code + respi_code + smkcode,
  data   = copd,
  family = binomial(link = "logit")
)

summary(model_risk)

# Tidy output with odds ratios
risk_tidy <- tidy(model_risk, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term      = recode(term,
      age            = "Age",
      gender_sex     = "Gender (Female)",
      genecode       = "Genetic (Alpha-1 AT)",
      obese_code     = "Obesity",
      poisongas_code = "Toxic Gas Exposure",
      respi_code     = "Respiratory Disease",
      smkcode        = "Smoking"
    ),
    significant = ifelse(p.value < 0.05, "Significant (p<0.05)", "Not significant")
  )

cat("\nOdds Ratios — Risk Factors:\n")
print(risk_tidy[, c("term", "estimate", "conf.low", "conf.high", "p.value")])

# Forest plot — Risk factors
p_risk <- ggplot(risk_tidy, aes(x = estimate, y = reorder(term, estimate),
                                 colour = significant)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.25, linewidth = 0.8) +
  scale_colour_manual(values = c("Significant (p<0.05)" = "#C0392B",
                                  "Not significant"      = "#85929E")) +
  scale_x_log10() +
  labs(
    title    = "Odds Ratios: Risk Factors for Severe COPD",
    subtitle = "Logistic Regression | Simulated SAIL Databank Wales cohort (n = 42,005)",
    x        = "Odds Ratio (log scale)",
    y        = NULL,
    colour   = NULL,
    caption  = "Source: Adewole (2024) — Swansea University MSc Dissertation"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey40"),
    legend.position = "bottom"
  )

ggsave("reports/04_logistic_risk_factors_OR.png", p_risk,
       width = 9, height = 5, dpi = 180)
cat("✅ Saved: reports/04_logistic_risk_factors_OR.png\n")

# ─────────────────────────────────────────────────────────────────────────────
# MODEL 2: Comorbidities
# ─────────────────────────────────────────────────────────────────────────────
cat("\n--- Model 2: Comorbidities ---\n")

model_comorbid <- glm(
  severe_tag ~ age + gender_sex + anxiety_code + asthma_code +
               diabetes_code + heart_code + hyper_code + osteo_code,
  data   = copd,
  family = binomial(link = "logit")
)

summary(model_comorbid)

comorbid_tidy <- tidy(model_comorbid, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = recode(term,
      age           = "Age",
      gender_sex    = "Gender (Female)",
      anxiety_code  = "Anxiety",
      asthma_code   = "Asthma",
      diabetes_code = "Diabetes",
      heart_code    = "Heart Disease",
      hyper_code    = "Hypertension",
      osteo_code    = "Osteoporosis"
    ),
    significant = ifelse(p.value < 0.05, "Significant (p<0.05)", "Not significant")
  )

cat("\nOdds Ratios — Comorbidities:\n")
print(comorbid_tidy[, c("term", "estimate", "conf.low", "conf.high", "p.value")])

p_comorbid <- ggplot(comorbid_tidy, aes(x = estimate, y = reorder(term, estimate),
                                         colour = significant)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.25, linewidth = 0.8) +
  scale_colour_manual(values = c("Significant (p<0.05)" = "#1A5276",
                                  "Not significant"      = "#85929E")) +
  scale_x_log10() +
  labs(
    title    = "Odds Ratios: Comorbidities Associated with Severe COPD",
    subtitle = "Logistic Regression | Simulated SAIL Databank Wales cohort (n = 42,005)",
    x        = "Odds Ratio (log scale)",
    y        = NULL,
    colour   = NULL,
    caption  = "Source: Adewole (2024) — Swansea University MSc Dissertation"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(colour = "grey40"),
    legend.position = "bottom"
  )

ggsave("reports/05_logistic_comorbidities_OR.png", p_comorbid,
       width = 9, height = 5, dpi = 180)
cat("✅ Saved: reports/05_logistic_comorbidities_OR.png\n")

# ─────────────────────────────────────────────────────────────────────────────
# ROC Curve & AUC Comparison
# ─────────────────────────────────────────────────────────────────────────────
cat("\n--- ROC / AUC ---\n")

roc_risk     <- roc(copd$severe_tag, fitted(model_risk),     quiet = TRUE)
roc_comorbid <- roc(copd$severe_tag, fitted(model_comorbid), quiet = TRUE)

cat(sprintf("AUC — Risk Factors   : %.3f\n", auc(roc_risk)))
cat(sprintf("AUC — Comorbidities  : %.3f\n", auc(roc_comorbid)))

# Build ROC plot data
roc_df <- bind_rows(
  data.frame(fpr = 1 - roc_risk$specificities,
             tpr = roc_risk$sensitivities,
             Model = sprintf("Risk Factors (AUC = %.3f)", auc(roc_risk))),
  data.frame(fpr = 1 - roc_comorbid$specificities,
             tpr = roc_comorbid$sensitivities,
             Model = sprintf("Comorbidities (AUC = %.3f)", auc(roc_comorbid)))
)

p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = Model)) +
  geom_line(linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  scale_colour_manual(values = c("#C0392B", "#1A5276")) +
  labs(
    title    = "ROC Curves: Logistic Regression Models for Severe COPD",
    subtitle = "Comparing predictive performance of Risk Factors vs Comorbidities",
    x        = "False Positive Rate (1 - Specificity)",
    y        = "True Positive Rate (Sensitivity)",
    colour   = "Model",
    caption  = "Source: Adewole (2024) — Swansea University MSc Dissertation"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("reports/06_roc_curves.png", p_roc,
       width = 8, height = 6, dpi = 180)
cat("✅ Saved: reports/06_roc_curves.png\n")
