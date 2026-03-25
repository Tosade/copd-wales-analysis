# =============================================================================
# 01_descriptive_stats.R
# Exploratory Data Analysis — COPD Wales Simulated Dataset
# Author: Oluwatosin Samson Adewole
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# --- Load data ---
copd <- read.csv("data/copd_simulated.csv")

cat("=== COPD WALES — DESCRIPTIVE STATISTICS ===\n\n")

# --- 1. Overview ---
cat("--- Dataset Overview ---\n")
cat(sprintf("Total patients     : %d\n", nrow(copd)))
cat(sprintf("Severe COPD cases  : %d (%.1f%%)\n",
            sum(copd$severe_tag), mean(copd$severe_tag)*100))
cat(sprintf("Deaths             : %d (%.1f%%)\n",
            sum(copd$death_flag), mean(copd$death_flag)*100))
cat(sprintf("Hospital admissions: %d (%.1f%%)\n",
            sum(copd$admission_flag), mean(copd$admission_flag)*100))
cat(sprintf("Mean age           : %.1f years (SD: %.1f)\n",
            mean(copd$age), sd(copd$age)))
cat(sprintf("Male patients      : %.1f%%\n", mean(copd$gender_sex == 0)*100))

# --- 2. Risk Factor Prevalence Plot ---
risk_factors <- copd %>%
  summarise(
    Smoking          = mean(smkcode) * 100,
    Obesity          = mean(obese_code) * 100,
    `Passive Smoking` = mean(passive_code) * 100,
    `Resp. Disease`  = mean(respi_code) * 100,
    `Genetic (A1AT)` = mean(genecode) * 100,
    `Toxic Gas Exp.` = mean(poisongas_code) * 100
  ) %>%
  pivot_longer(everything(), names_to = "Factor", values_to = "Prevalence")

p1 <- ggplot(risk_factors, aes(x = reorder(Factor, Prevalence),
                                y = Prevalence, fill = Prevalence)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", Prevalence)),
            hjust = -0.1, size = 3.5, fontface = "bold") +
  scale_fill_gradient(low = "#AED6F1", high = "#1A5276") +
  coord_flip(ylim = c(0, 80)) +
  labs(
    title    = "Prevalence of COPD Risk Factors in Wales (2015–2022)",
    subtitle = "Simulated dataset based on SAIL Databank cohort (n = 42,005)",
    x        = NULL,
    y        = "Prevalence (%)",
    caption  = "Source: Adewole (2024) — Swansea University MSc Dissertation"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(colour = "grey40"),
    axis.text.y   = element_text(face = "bold")
  )

ggsave("reports/01_risk_factor_prevalence.png", p1,
       width = 9, height = 5, dpi = 180)
cat("\n✅ Saved: reports/01_risk_factor_prevalence.png\n")

# --- 3. Comorbidity Prevalence Plot ---
comorbidities <- copd %>%
  summarise(
    Asthma        = mean(asthma_code) * 100,
    Anxiety       = mean(anxiety_code) * 100,
    Hypertension  = mean(hyper_code) * 100,
    `Heart Disease` = mean(heart_code) * 100,
    Diabetes      = mean(diabetes_code) * 100,
    Osteoporosis  = mean(osteo_code) * 100,
    Tuberculosis  = mean(tb_code) * 100
  ) %>%
  pivot_longer(everything(), names_to = "Comorbidity", values_to = "Prevalence")

p2 <- ggplot(comorbidities, aes(x = reorder(Comorbidity, Prevalence),
                                 y = Prevalence, fill = Prevalence)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", Prevalence)),
            hjust = -0.1, size = 3.5, fontface = "bold") +
  scale_fill_gradient(low = "#A9DFBF", high = "#1E8449") +
  coord_flip(ylim = c(0, 55)) +
  labs(
    title    = "Prevalence of Comorbidities in COPD Patients in Wales",
    subtitle = "Simulated dataset based on SAIL Databank cohort (n = 42,005)",
    x        = NULL,
    y        = "Prevalence (%)",
    caption  = "Source: Adewole (2024) — Swansea University MSc Dissertation"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(colour = "grey40"),
    axis.text.y   = element_text(face = "bold")
  )

ggsave("reports/02_comorbidity_prevalence.png", p2,
       width = 9, height = 5, dpi = 180)
cat("✅ Saved: reports/02_comorbidity_prevalence.png\n")

# --- 4. Age distribution by severity ---
p3 <- ggplot(copd, aes(x = age, fill = factor(severe_tag))) +
  geom_histogram(bins = 35, alpha = 0.75, position = "identity") +
  scale_fill_manual(values = c("#85C1E9", "#C0392B"),
                    labels = c("Non-severe", "Severe / Death")) +
  labs(
    title    = "Age Distribution by COPD Severity",
    subtitle = "Severe cases concentrated in older age groups (45+)",
    x        = "Age (years)",
    y        = "Count",
    fill     = "Severity",
    caption  = "Source: Adewole (2024) — Swansea University MSc Dissertation"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("reports/03_age_distribution_by_severity.png", p3,
       width = 9, height = 5, dpi = 180)
cat("✅ Saved: reports/03_age_distribution_by_severity.png\n")

# --- 5. Severity rate by gender ---
gender_summary <- copd %>%
  group_by(gender_sex) %>%
  summarise(
    n          = n(),
    severe_pct = mean(severe_tag) * 100,
    death_pct  = mean(death_flag) * 100
  ) %>%
  mutate(gender_label = ifelse(gender_sex == 0, "Male", "Female"))

cat("\n--- Severity & Mortality by Gender ---\n")
print(gender_summary[, c("gender_label", "n", "severe_pct", "death_pct")])
