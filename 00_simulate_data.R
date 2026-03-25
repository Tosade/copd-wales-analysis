# =============================================================================
# 00_simulate_data.R
# Simulate COPD dataset based on SAIL Databank Wales study (2015-2022)
# Author: Oluwatosin Samson Adewole
# Note: Real data from SAIL Databank is confidential. This script generates
#       a synthetic dataset that mirrors the structure and distributions
#       observed in the original study for reproducibility and sharing.
# =============================================================================

set.seed(42)
n <- 42005  # Focus group size matching thesis cohort

# --- Demographics ---
age        <- round(rnorm(n, mean = 65, sd = 12))
age        <- pmax(pmin(age, 95), 35)   # Clamp to realistic COPD age range
gender_sex <- rbinom(n, 1, prob = 0.45) # 0 = Male, 1 = Female (males more prevalent)

# --- Risk Factors ---
# Probabilities calibrated to thesis findings
smkcode      <- rbinom(n, 1, prob = 0.62)  # Smoking: dominant risk factor
obese_code   <- rbinom(n, 1, prob = 0.28)  # Obesity
passive_code <- rbinom(n, 1, prob = 0.18)  # Passive smoking
respi_code   <- rbinom(n, 1, prob = 0.35)  # Other respiratory diseases
genecode     <- rbinom(n, 1, prob = 0.04)  # Alpha-1 antitrypsin deficiency (rare)
poisongas_code <- rbinom(n, 1, prob = 0.03) # Toxic gas exposure (rare)

# --- Comorbidities ---
asthma_code    <- rbinom(n, 1, prob = 0.40)
anxiety_code   <- rbinom(n, 1, prob = 0.38)
diabetes_code  <- rbinom(n, 1, prob = 0.22)
heart_code     <- rbinom(n, 1, prob = 0.25)
hyper_code     <- rbinom(n, 1, prob = 0.30)
osteo_code     <- rbinom(n, 1, prob = 0.15)
tb_code        <- rbinom(n, 1, prob = 0.02)

# --- Severity tag (outcome variable) ---
# Binary: 1 = severe COPD / death, 0 = non-severe
# Modelled using logistic function based on thesis coefficients
log_odds_severity <- (
  -5.33
  + 0.055  * age
  - 0.349  * gender_sex
  + 0.541  * genecode
  - 0.067  * obese_code
  + 1.003  * smkcode
  - 0.079  * respi_code
  - 0.417  * poisongas_code
  + 0.308  * diabetes_code
  + 0.244  * heart_code
  + 0.095  * anxiety_code
  - 2.648  * asthma_code
  + 0.035  * osteo_code
  - 0.064  * hyper_code
  + rnorm(n, 0, 0.5)  # Random noise
)
prob_severe  <- 1 / (1 + exp(-log_odds_severity))
severe_tag   <- rbinom(n, 1, prob = prob_severe)

# --- Survival / Time-to-event ---
# Study period: 2015-01-01 to 2022-03-02 (approx 2618 days)
study_start <- as.Date("2015-01-01")
study_end   <- as.Date("2022-03-02")
total_days  <- as.numeric(study_end - study_start)

# Base hazard influenced by key factors (smoking, age, heart disease)
base_lambda <- exp(
  -8
  + 0.047  * age
  + 0.763  * smkcode
  - 0.303  * passive_code
  + 0.154  * heart_code
  + 0.194  * diabetes_code
  + 0.065  * osteo_code
  - 0.239  * gender_sex
)

time_to_event <- rexp(n, rate = base_lambda)
time_to_event <- pmin(time_to_event, total_days)  # Censor at end of study
event_status  <- as.integer(time_to_event < total_days & severe_tag == 1)

# --- Admission & death flags ---
death_flag     <- rbinom(n, 1, prob = severe_tag * 0.32)
admission_flag <- rbinom(n, 1, prob = 0.45 + severe_tag * 0.3)

# --- Compile dataset ---
copd_data <- data.frame(
  patient_id     = 1:n,
  age            = age,
  gender_sex     = gender_sex,         # 0=Male, 1=Female
  smkcode        = smkcode,
  obese_code     = obese_code,
  passive_code   = passive_code,
  respi_code     = respi_code,
  genecode       = genecode,
  poisongas_code = poisongas_code,
  asthma_code    = asthma_code,
  anxiety_code   = anxiety_code,
  diabetes_code  = diabetes_code,
  heart_code     = heart_code,
  hyper_code     = hyper_code,
  osteo_code     = osteo_code,
  tb_code        = tb_code,
  severe_tag     = severe_tag,
  time_to_event  = round(time_to_event),
  event_status   = event_status,
  death_flag     = death_flag,
  admission_flag = admission_flag
)

# --- Save ---
write.csv(copd_data, file = "copd_simulated.csv", row.names = FALSE)

cat("✅ Simulated dataset created: copd_simulated.csv\n")
cat(sprintf("   Total patients : %d\n", n))
cat(sprintf("   Severe cases   : %d (%.1f%%)\n", sum(copd_data$severe_tag),
            mean(copd_data$severe_tag) * 100))
cat(sprintf("   Deaths         : %d (%.1f%%)\n", sum(copd_data$death_flag),
            mean(copd_data$death_flag) * 100))
cat(sprintf("   Events (Cox)   : %d (%.1f%%)\n", sum(copd_data$event_status),
            mean(copd_data$event_status) * 100))
