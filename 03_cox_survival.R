# =============================================================================
# 03_cox_survival.R
# Cox Proportional Hazards Model + Kaplan-Meier Survival Analysis
# Author: Oluwatosin Samson Adewole
# Mirrors methodology from MSc dissertation (Swansea University, 2024)
# =============================================================================

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(broom)

# --- Load data ---
copd <- read.csv("copd_simulated.csv")

cat("=== SURVIVAL ANALYSIS ===\n\n")

# Create survival object
surv_obj <- Surv(time = copd$time_to_event, event = copd$event_status)

# ─────────────────────────────────────────────────────────────────────────────
# KAPLAN-MEIER CURVES
# ─────────────────────────────────────────────────────────────────────────────

# --- KM by Smoking Status ---
cat("--- Kaplan-Meier: Smoking Status ---\n")
km_smoking <- survfit(surv_obj ~ smkcode, data = copd)

p_km_smoking <- ggsurvplot(
  km_smoking,
  data         = copd,
  palette      = c("#2E86C1", "#C0392B"),
  legend.labs  = c("Non-Smoker", "Smoker"),
  legend.title = "Smoking Status",
  pval         = TRUE,
  risk.table   = TRUE,
  risk.table.height = 0.28,
  conf.int     = TRUE,
  xlab         = "Time (days)",
  ylab         = "Survival Probability",
  title        = "Kaplan-Meier Survival Curve by Smoking Status",
  subtitle     = "Simulated SAIL Databank Wales cohort (n = 42,005)",
  ggtheme      = theme_minimal(base_size = 12),
  fontsize     = 3.5
)

ggsave("reports/07_km_smoking.png",
       print(p_km_smoking), width = 10, height = 7, dpi = 180)
cat("✅ Saved: reports/07_km_smoking.png\n")

# --- KM by Age Group ---
cat("\n--- Kaplan-Meier: Age Groups ---\n")
copd <- copd %>%
  mutate(age_group = case_when(
    age < 50            ~ "< 50",
    age >= 50 & age < 65 ~ "50–64",
    age >= 65 & age < 75 ~ "65–74",
    age >= 75           ~ "75+"
  ))
copd$age_group <- factor(copd$age_group,
                          levels = c("< 50", "50–64", "65–74", "75+"))

km_age <- survfit(surv_obj ~ age_group, data = copd)

p_km_age <- ggsurvplot(
  km_age,
  data         = copd,
  palette      = c("#27AE60", "#F39C12", "#E67E22", "#C0392B"),
  legend.title = "Age Group",
  pval         = TRUE,
  conf.int     = FALSE,
  risk.table   = TRUE,
  risk.table.height = 0.30,
  xlab         = "Time (days)",
  ylab         = "Survival Probability",
  title        = "Kaplan-Meier Survival Curve by Age Group",
  subtitle     = "Older patients show markedly worse survival — consistent with thesis findings",
  ggtheme      = theme_minimal(base_size = 12),
  fontsize     = 3.5
)

ggsave("reports/08_km_age_groups.png",
       print(p_km_age), width = 10, height = 7, dpi = 180)
cat("✅ Saved: reports/08_km_age_groups.png\n")

# --- KM by Gender ---
cat("\n--- Kaplan-Meier: Gender ---\n")
km_gender <- survfit(surv_obj ~ gender_sex, data = copd)

p_km_gender <- ggsurvplot(
  km_gender,
  data         = copd,
  palette      = c("#2980B9", "#C0392B"),
  legend.labs  = c("Male", "Female"),
  legend.title = "Gender",
  pval         = TRUE,
  conf.int     = TRUE,
  risk.table   = TRUE,
  risk.table.height = 0.28,
  xlab         = "Time (days)",
  ylab         = "Survival Probability",
  title        = "Kaplan-Meier Survival Curve by Gender",
  subtitle     = "Male patients show higher hazard — reflected in Cox model HR",
  ggtheme      = theme_minimal(base_size = 12),
  fontsize      = 3.5
)

ggsave("reports/09_km_gender.png",
       print(p_km_gender), width = 10, height = 7, dpi = 180)
cat("✅ Saved: reports/09_km_gender.png\n")

# ─────────────────────────────────────────────────────────────────────────────
# COX PROPORTIONAL HAZARDS — Risk Factors
# ─────────────────────────────────────────────────────────────────────────────
cat("\n--- Cox Model: Risk Factors ---\n")

cox_risk <- coxph(
  surv_obj ~ age + gender_sex + genecode + obese_code +
             passive_code + respi_code + smkcode,
  data = copd
)

summary(cox_risk)

cox_risk_tidy <- tidy(cox_risk, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    term = recode(term,
      age            = "Age",
      gender_sex     = "Gender (Female)",
      genecode       = "Genetic (Alpha-1 AT)",
      obese_code     = "Obesity",
      passive_code   = "Passive Smoking",
      respi_code     = "Respiratory Disease",
      smkcode        = "Smoking"
    ),
    significant = ifelse(p.value < 0.05, "Yes", "No")
  )

cat("\nHazard Ratios — Risk Factors:\n")
print(cox_risk_tidy[, c("term", "estimate", "conf.low", "conf.high", "p.value")])

# Forest plot
p_cox_risk <- ggplot(cox_risk_tidy,
                      aes(x = estimate, y = reorder(term, estimate),
                          colour = significant)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.25, linewidth = 0.8) +
  scale_colour_manual(values = c("Yes" = "#C0392B", "No" = "#85929E"),
                      labels = c("Yes" = "Significant (p<0.05)",
                                 "No"  = "Not significant")) +
  scale_x_log10() +
  labs(
    title    = "Hazard Ratios: Cox Model — Risk Factors for COPD Severity",
    subtitle = sprintf("Concordance = %.3f | Simulated SAIL Databank Wales cohort",
                       summary(cox_risk)$concordance["C"]),
    x        = "Hazard Ratio (log scale)",
    y        = NULL,
    colour   = NULL,
    caption  = "Source: Adewole (2024) — Swansea University MSc Dissertation"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("reports/10_cox_risk_factors_HR.png", p_cox_risk,
       width = 9, height = 5, dpi = 180)
cat("✅ Saved: reports/10_cox_risk_factors_HR.png\n")

# ─────────────────────────────────────────────────────────────────────────────
# COX PROPORTIONAL HAZARDS — Comorbidities
# ─────────────────────────────────────────────────────────────────────────────
cat("\n--- Cox Model: Comorbidities ---\n")

cox_comorbid <- coxph(
  surv_obj ~ age + gender_sex + anxiety_code + asthma_code +
             diabetes_code + heart_code + hyper_code + osteo_code,
  data = copd
)

summary(cox_comorbid)

cox_comorbid_tidy <- tidy(cox_comorbid, exponentiate = TRUE, conf.int = TRUE) %>%
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
    significant = ifelse(p.value < 0.05, "Yes", "No")
  )

cat("\nHazard Ratios — Comorbidities:\n")
print(cox_comorbid_tidy[, c("term", "estimate", "conf.low", "conf.high", "p.value")])

p_cox_comorbid <- ggplot(cox_comorbid_tidy,
                          aes(x = estimate, y = reorder(term, estimate),
                              colour = significant)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.25, linewidth = 0.8) +
  scale_colour_manual(values = c("Yes" = "#1A5276", "No" = "#85929E"),
                      labels = c("Yes" = "Significant (p<0.05)",
                                 "No"  = "Not significant")) +
  scale_x_log10() +
  labs(
    title    = "Hazard Ratios: Cox Model — Comorbidities & COPD Severity",
    subtitle = sprintf("Concordance = %.3f | Simulated SAIL Databank Wales cohort",
                       summary(cox_comorbid)$concordance["C"]),
    x        = "Hazard Ratio (log scale)",
    y        = NULL,
    colour   = NULL,
    caption  = "Source: Adewole (2024) — Swansea University MSc Dissertation"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("reports/11_cox_comorbidities_HR.png", p_cox_comorbid,
       width = 9, height = 5, dpi = 180)
cat("✅ Saved: reports/11_cox_comorbidities_HR.png\n")

cat("\n=== Survival analysis complete ===\n")
