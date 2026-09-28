# =============================================================================
# analyses.R
# All statistical results shown in the dashboard AND quoted in the blog.
# Requires: `nh` loaded, plus broom, Epi, and the tidyverse.
# All intervals use a 90% confidence level.
# =============================================================================

conf_level <- 0.90
alpha      <- 1 - conf_level

# Percentile bootstrap CI for a difference in means (group 1 minus group 2)
boot_mean_diff <- function(y, g, reps = 2000, conf = conf_level, seed = 431) {
  set.seed(seed)
  lv <- levels(droplevels(g))
  y1 <- y[g == lv[1]]
  y2 <- y[g == lv[2]]
  diffs <- replicate(reps, mean(sample(y1, replace = TRUE)) -
                       mean(sample(y2, replace = TRUE)))
  tibble(
    estimate  = mean(y1) - mean(y2),
    conf.low  = unname(quantile(diffs, (1 - conf) / 2)),
    conf.high = unname(quantile(diffs, 1 - (1 - conf) / 2)),
    p.value   = NA_real_
  )
}

# Proportion of a factor equal to a given level, as a percentage
pct <- function(x, level) 100 * mean(x == level)

kpi <- list(
  n          = nrow(nh),
  median_ldl = median(nh$ldl),
  pct_ldl130 = pct(nh$ldl_high, "Yes"),
  pct_htn    = pct(nh$htn_range, "Yes"),
  r_ldl_sbp  = cor(nh$ldl, nh$sbp)
)

b_summary <- nh |>
  group_by(ldl_high) |>
  summarise(n = n(), mean = mean(sbp), sd = sd(sbp), median = median(sbp),
            .groups = "drop")

b_results <- bind_rows(
  "Pooled t (equal variances)" =
    tidy(t.test(sbp ~ ldl_high, data = nh, var.equal = TRUE, conf.level = conf_level)),
  "Welch t (unequal variances)" =
    tidy(t.test(sbp ~ ldl_high, data = nh, conf.level = conf_level)),
  "Bootstrap (difference in means)" =
    boot_mean_diff(nh$sbp, nh$ldl_high),
  "Wilcoxon rank sum (location shift)" =
    tidy(wilcox.test(sbp ~ ldl_high, data = nh, conf.int = TRUE,
                     conf.level = conf_level, exact = FALSE)),
  .id = "Method"
) |>
  select(Method, estimate, conf.low, conf.high, p.value)

c_summary <- nh |>
  group_by(ldl_cat) |>
  summarise(n = n(), mean = mean(sbp), sd = sd(sbp), median = median(sbp),
            .groups = "drop")

c_aov     <- aov(sbp ~ ldl_cat, data = nh)
c_anova   <- tidy(c_aov)
c_r2      <- summary(lm(sbp ~ ldl_cat, data = nh))$r.squared
c_tukey   <- tidy(TukeyHSD(c_aov, conf.level = conf_level)) |>
  select(contrast, estimate, conf.low, conf.high, adj.p.value)
c_kruskal <- tidy(kruskal.test(sbp ~ ldl_cat, data = nh))

d_table <- table(`LDL >= 130` = nh$ldl_high, `BP >= 130/80` = nh$htn_range)
d_epi   <- twoby2(d_table, alpha = alpha, print = FALSE)

d_results <- tibble(
  Measure   = c("Risk (Pr BP >= 130/80) if LDL >= 130",
                "Risk (Pr BP >= 130/80) if LDL < 130",
                "Relative risk", "Odds ratio", "Risk difference"),
  estimate  = c(d_epi$table[1, 3], d_epi$table[2, 3], d_epi$measures[c(1, 2, 4), 1]),
  conf.low  = c(d_epi$table[1, 4], d_epi$table[2, 4], d_epi$measures[c(1, 2, 4), 2]),
  conf.high = c(d_epi$table[1, 5], d_epi$table[2, 5], d_epi$measures[c(1, 2, 4), 3])
)
d_pvalue <- d_epi$p.value[2]

e_table    <- table(nh$race_eth, nh$bp_cat)
e_min_cell <- min(e_table)
e_chisq    <- chisq.test(e_table)
e_tidy     <- tidy(e_chisq)
e_props    <- nh |>
  count(race_eth, bp_cat) |>
  group_by(race_eth) |>
  mutate(pct = 100 * n / sum(n)) |>
  ungroup()

m_main <- lm(sbp ~ ldl10 + age + sex + race_eth, data = nh)
m_int  <- lm(sbp ~ ldl10 * (sex + race_eth) + age, data = nh)

main_slope <- tidy(m_main, conf.int = TRUE, conf.level = conf_level) |>
  filter(term == "ldl10")

interaction_tests <- drop1(m_int, test = "F") |>
  tidy() |>
  filter(str_detect(term, ":"))

strata_slopes <- nh |>
  group_by(sex, race_eth) |>
  group_modify(~ tidy(lm(sbp ~ ldl10 + age, data = .x),
                      conf.int = TRUE, conf.level = conf_level) |>
                 mutate(n = nrow(.x))) |>
  ungroup() |>
  filter(term == "ldl10") |>
  select(sex, race_eth, n, estimate, conf.low, conf.high, p.value)

