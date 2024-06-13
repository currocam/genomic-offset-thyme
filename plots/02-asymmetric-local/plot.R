source("plots/theme.R")
library(ggsci)

adjusted_r2 <- function(x, y){
  model <- lm(y ~ x)
  summary(model)$adj.r.squared
}

dat <- read_csv("results/local_adaptation_scenarios/asymmetric_fitness_variance.csv.gz")

dat2 <- dat |>
  group_by(model, seed, QTLs) |>
  summarise(
    "Causal" = adjusted_r2(-log(future_fitness), causal_offset),
    "Empirical" = adjusted_r2(-log(future_fitness), empirical),
    "All SNPs" = adjusted_r2(-log(future_fitness), all_snps_offset),
    "Incomplete (5+5 QTLs)" = adjusted_r2(-log(future_fitness), incomplete_offset),
    "Missing secondary trait (5+0 QTLs)" = adjusted_r2(-log(future_fitness), only_QTL1_offset),
    "Missing primary trait (0+5 QTLs)" = adjusted_r2(-log(future_fitness), only_QTL2_offset),
    "Environmental distance" =  adjusted_r2(-log(future_fitness), env_dist),
  )

p1 <- dat2 |>
  pivot_longer(-c(model, seed, QTLs), names_to = "R2adj") |>
  mutate(
    model = factor(model, levels = c("m3.3", "m3.4"), labels = c("Weak asymmetry", "Strong asymmetry")),
    R2adj = factor(R2adj, levels = colnames(dat2)[-c(1:3)]) 
    ) |>
  ggplot(aes(
    y = R2adj,
    x = value, colour = model)
  ) +
  geom_boxplot() +
  theme_classic() +
  scale_colour_frontiers(name = "") +
  theme(legend.position = "bottom")+
  xlab("Adjusted R squared") +
  ylab("")

ggsave(
  "plots/02-asymmetric-local/robustness.pdf",
  p1, device = cairo_pdf,
  width = fig.witdh, height = fig.height, units = "mm", dpi = "retina"
  )
