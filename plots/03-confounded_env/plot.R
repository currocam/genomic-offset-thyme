source("plots/theme.R")
library(ggsci)
adjusted_r2 <- function(x, y){
  model <- lm(y ~ x)
  summary(model)$adj.r.squared
}

dat <- read_csv("results/local_adaptation_scenarios/m4_offsets_uncorrelated.csv")
dat2 <- dat |>
  drop_na() |>
  pivot_longer(ends_with("_offset"), values_to = "offset", names_to = "type") |>
  group_by(model, seed, QTLs, type) |>
  summarise(
    R2adj = adjusted_r2(-log(future_fitness), offset)
  )

p1 <- dat2 |>
  filter(startsWith(type, "causal_random")) |>
  separate(type, into = c(NA, NA, "n", NA), sep = "_")|>
  mutate(n = as.numeric(n)) |>
  ggplot(aes(x = n, y = R2adj))+
  geom_jitter(aes(colour = as.factor(QTLs)), width=1.5, shape = 1) +
  geom_boxplot(aes(group = n), width=3.5, outlier.shape = NA, alpha = 0) +
  theme_classic() +
  scale_colour_npg(name = "Number of QTLs") +
  theme(legend.position = "bottom")+
  ylab("Adjusted R squared") +
  xlab("Number of added uncorrelated environmental factors")+
  ylim(c(0, 1))

dat <- read_csv("results/local_adaptation_scenarios/m4_offsets_confounded.csv.gz") |>
  drop_na()

dat2 <- dat |>
  group_by(model, seed, QTLs, mobility, fst) |>
  summarise(
    empirical = adjusted_r2(-log(future_fitness), empirical_offset),
    confounded = adjusted_r2(-log(future_fitness), empirical_confounded_offset)
  ) |>
  drop_na() |>
  mutate(diff = empirical - confounded)

p2 <- dat2 |>
  #filter(mobility > 0.10) |>
  rename("Not confounded" = empirical, "Confounded" = confounded) |>
  pivot_longer(c("Not confounded", Confounded), values_to = "val", names_to = "type") |>
  ggplot(aes(x = fst, y = val, colour = type, group = interaction(type, cut_interval(fst, n = 3)))) +
  geom_smooth(
    aes(x = fst, y = val, colour = type), se = TRUE,
    linetype = "dashed", method ="lm", inherit.aes = FALSE)+
  geom_point(aes(colour = type), shape = 1)+
  #geom_boxplot(outlier.shape = NA)+
  theme_classic() +
  scale_colour_jama(name = "Geometric genomic offset") +
  scale_fill_jama(name = "Geometric genomic offset") +
  theme(legend.position = "bottom") +
  ylab("Adjusted R squared")+
  xlab("Fst")

p2

ggsave(
  "plots/03-confounded_env/uncorrelated.pdf",
  p1, device = cairo_pdf,
  width = fig.witdh, height = fig.height, units = "mm", dpi = "retina"
  )

ggsave(
  "plots/03-confounded_env/confounded.pdf",
  p2, device = cairo_pdf,
  width = fig.witdh, height = fig.height, units = "mm", dpi = "retina"
  )
