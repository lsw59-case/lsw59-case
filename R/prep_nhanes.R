# =============================================================================
# prep_nhanes.R
# Ingest, merge, and clean NHANES 2017-March 2020 pre-pandemic data
#   P_DEMO   Demographics
#   P_TRIGLY LDL cholesterol (lab; fasting subsample only)
#   P_BPXO   Oscillometric blood pressure (exam)
# =============================================================================

library(here)
library(nhanesA)
library(janitor)
library(tidyverse)   # tidyverse loaded last, per course convention

raw_dir <- here("data", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

get_nhanes <- function(table) {
  path <- file.path(raw_dir, paste0(table, ".rds"))
  if (!file.exists(path)) {
    nhanes(table, translated = FALSE) |>
      as_tibble() |>
      saveRDS(path)
  }
  readRDS(path)
}

demo_raw <- get_nhanes("P_DEMO")
trig_raw <- get_nhanes("P_TRIGLY")
bpxo_raw <- get_nhanes("P_BPXO")

dim(demo_raw); dim(trig_raw); dim(bpxo_raw)
names(trig_raw)
names(bpxo_raw)

demo <- demo_raw |>
  select(SEQN, RIDSTATR, RIAGENDR, RIDAGEYR, RIDRETH3)

trig <- trig_raw |>
  select(SEQN, LBDLDL)                      # LDL-C, Friedewald, mg/dL

bpxo <- bpxo_raw |>
  select(SEQN,
         BPXOSY1, BPXOSY2, BPXOSY3,         # systolic readings 1-3
         BPXODI1, BPXODI2, BPXODI3)         # diastolic readings 1-3

merged <- demo |>
  left_join(trig, by = "SEQN") |>
  left_join(bpxo, by = "SEQN") |>
  clean_names()

nrow(merged) == nrow(demo_raw)       # should be TRUE
n_distinct(merged$seqn) == nrow(merged)  # TRUE = no duplicate people

sbp_vars <- c("bpxosy1", "bpxosy2", "bpxosy3")
dbp_vars <- c("bpxodi1", "bpxodi2", "bpxodi3")

eligible <- merged |>
  filter(ridstatr == 2, between(ridageyr, 45, 79)) |>
  mutate(
    n_bp = rowSums(!is.na(pick(all_of(sbp_vars)))),
    sbp  = if_else(n_bp >= 2,
                   rowMeans(pick(all_of(sbp_vars)), na.rm = TRUE), NA_real_),
    dbp  = if_else(n_bp >= 2,
                   rowMeans(pick(all_of(dbp_vars)), na.rm = TRUE), NA_real_),
    sbp  = round(sbp, 1),
    dbp  = round(dbp, 1),
    
    age  = ridageyr,
    ldl  = lbdldl,
    
    sex = factor(riagendr, levels = c(1, 2), labels = c("Male", "Female")),
    
    race_eth = case_when(
      ridreth3 %in% c(1, 2) ~ "Hispanic",
      ridreth3 == 3         ~ "NH White",
      ridreth3 == 4         ~ "NH Black",
      ridreth3 == 6         ~ "NH Asian",
      ridreth3 == 7         ~ "Other/Multiracial"
    ),
    race_eth = factor(race_eth, levels = c("NH White", "NH Black", "Hispanic",
                                           "NH Asian", "Other/Multiracial"))
  )

count(eligible, sex)
count(eligible, race_eth)
summary(eligible$sbp)

nh <- eligible |>
  filter(!is.na(ldl), !is.na(sbp), !is.na(dbp),
         !is.na(sex), !is.na(race_eth), !is.na(age)) |>
  mutate(
    ldl_high = factor(if_else(ldl >= 130, "Yes", "No"), levels = c("Yes", "No")),
    
    ldl_cat = cut(ldl,
                  breaks = c(-Inf, 100, 130, 160, Inf), right = FALSE,
                  labels = c("<100", "100-129", "130-159", "160+")),
    
    bp_cat = case_when(
      sbp >= 140 | dbp >= 90 ~ "Stage 2",
      sbp >= 130 | dbp >= 80 ~ "Stage 1",
      TRUE                   ~ "Normal/Elevated"
    ),
    bp_cat = factor(bp_cat, levels = c("Normal/Elevated", "Stage 1", "Stage 2")),
    
    htn_range = factor(if_else(bp_cat == "Normal/Elevated", "No", "Yes"),
                       levels = c("Yes", "No")),
    
    ldl10 = ldl / 10
  ) |>
  select(seqn, age, sex, race_eth, ldl, ldl10, ldl_high, ldl_cat,
         sbp, dbp, bp_cat, htn_range)

nrow(nh)                        # must be between 500 and 7,500
count(nh, ldl_high, htn_range)  # every cell >= 30 (Analysis D)
count(nh, ldl_cat)              # every group >= 30 (Analysis C)
table(nh$race_eth, nh$bp_cat)   # every cell >= 15 (Analysis E)

stopifnot(
  "Duplicate SEQN values"         = !anyDuplicated(nh$seqn),
  "Sample outside 500-7,500 rule" = between(nrow(nh), 500, 7500),
  "LDL has < 15 unique values"    = n_distinct(nh$ldl) >= 15,
  "SBP has < 15 unique values"    = n_distinct(nh$sbp) >= 15,
  "Implausible LDL value"         = all(nh$ldl > 0 & nh$ldl < 500),
  "Implausible SBP value"         = all(between(nh$sbp, 60, 260)),
  "Implausible DBP value"         = all(between(nh$dbp, 20, 160))
)

# -----------------------------------------------------------------------------
# 5. Sample flow, codebook, and save
# -----------------------------------------------------------------------------

meta <- list(
  n_demo       = nrow(demo_raw),
  n_eligible   = nrow(eligible),
  n_with_ldl   = sum(!is.na(eligible$ldl)),
  n_complete   = nrow(nh),
  retrieved_on = as.Date(min(file.mtime(list.files(raw_dir, full.names = TRUE)))),
  prepared_on  = Sys.Date()
)

codebook <- tribble(
  ~Variable,   ~Type,    ~Description,                                                ~Original,
  "seqn",      "ID",     "Respondent sequence number",                                "SEQN",
  "age",       "Quant",  "Age at screening, years (45-79)",                           "RIDAGEYR",
  "sex",       "Binary", "Sex: Male, Female",                                         "RIAGENDR",
  "race_eth",  "5-cat",  "Race/ethnicity: NH White, NH Black, Hispanic, NH Asian, Other/Multiracial", "RIDRETH3",
  "ldl",       "Quant",  "LDL cholesterol, Friedewald equation (mg/dL)",              "LBDLDL",
  "ldl_high",  "Binary", "LDL of 130 mg/dL or higher: Yes, No",                       "from LBDLDL",
  "ldl_cat",   "4-cat",  "LDL category (mg/dL): <100, 100-129, 130-159, 160+",        "from LBDLDL",
  "sbp",       "Quant",  "Mean systolic BP of 2-3 oscillometric readings (mm Hg)",    "BPXOSY1-3",
  "dbp",       "Quant",  "Mean diastolic BP of 2-3 oscillometric readings (mm Hg)",   "BPXODI1-3",
  "bp_cat",    "3-cat",  "2017 ACC/AHA category: Normal/Elevated, Stage 1, Stage 2",  "from sbp, dbp",
  "htn_range", "Binary", "Measured BP of 130/80 or higher (Stage 1 or 2): Yes, No",   "from bp_cat"
)

saveRDS(nh,       here("data", "nh_ldl_bp.rds"))
saveRDS(meta,     here("data", "nh_meta.rds"))
saveRDS(codebook, here("data", "nh_codebook.rds"))


