# =============================================================================
# helpers.R
# Palettes, ggplot theme, number formatting, and plot builders shared by the
# dashboard and the blog.
# =============================================================================

# Color-blind-safe palettes (Okabe-Ito, plus a sequential blue for BP stage)
pal_sex  <- c(Male = "#0072B2", Female = "#D55E00")
pal_race <- c("NH White" = "#E69F00", "NH Black" = "#56B4E9",
              "Hispanic" = "#009E73", "NH Asian" = "#CC79A7",
              "Other/Multiracial" = "#555555")
pal_bp   <- c("Normal/Elevated" = "#C6DBEF", "Stage 1" = "#6BAED6",
              "Stage 2" = "#08519C")
pal_yesno <- c(Yes = "#08519C", No = "#C6DBEF")

theme_set(
  theme_minimal(base_size = 13) +
    theme(legend.position     = "bottom",
          panel.grid.minor    = element_blank(),
          plot.title.position = "plot")
)

lab_sbp <- "Mean systolic BP (mm Hg)"
lab_ldl <- "LDL cholesterol (mg/dL)"

fmt <- function(x, digits = 1) {
  formatC(x, format = "f", digits = digits, big.mark = ",")
}

fmt_ci <- function(est, lo, hi, digits = 1) {
  paste0(fmt(est, digits), " (90% CI: ", fmt(lo, digits), " to ", fmt(hi, digits), ")")
}

fmt_p <- function(p) if (p < 0.001) "p < 0.001" else paste0("p = ", fmt(p, 3))

ci_position <- function(lo, hi, null = 0) {
  if (lo > null) "lies entirely above" else if (hi < null) "lies entirely below" else "includes"
}



plot_ldl_sbp_scatter <- function(data) {
  ggplot(data, aes(x = ldl, y = sbp, color = sex)) +
    geom_point(alpha = 0.25, size = 1) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 1.2) +
    scale_color_manual(values = pal_sex) +
    labs(x = lab_ldl, y = lab_sbp, color = NULL)
}

plot_sbp_by_group <- function(data, group, xlab) {
  ggplot(data, aes(x = {{ group }}, y = sbp)) +
    geom_violin(fill = "grey92", color = NA) +
    geom_boxplot(width = 0.25, outlier.alpha = 0.3, fill = "white") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                 fill = "#D55E00", color = "black") +
    labs(x = xlab, y = lab_sbp,
         caption = "Diamond = group mean; box = median and IQR")
}

plot_stacked_pct <- function(data, x, fill, palette, xlab, filllab) {
  ggplot(data, aes(x = {{ x }}, fill = {{ fill }})) +
    geom_bar(position = "fill", color = "white") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = palette) +
    labs(x = xlab, y = "Percent of group", fill = filllab)
}

plot_sample_composition <- function(data) {
  ggplot(data, aes(y = fct_rev(race_eth), fill = sex)) +
    geom_bar(position = position_dodge(preserve = "single")) +
    scale_fill_manual(values = pal_sex) +
    labs(x = "Participants", y = NULL, fill = NULL)
}

plot_strata_forest <- function(slopes) {
  ggplot(slopes, aes(x = estimate, y = fct_rev(race_eth), color = sex)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high), linewidth = 0.9) +
    facet_wrap(~ sex) +
    scale_color_manual(values = pal_sex, guide = "none") +
    labs(x = "Change in mean SBP (mm Hg) per 10 mg/dL higher LDL,\nage-adjusted, with 90% CI",
         y = NULL)
}


