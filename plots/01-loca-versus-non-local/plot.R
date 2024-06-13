source("plots/theme.R")
library(ggsci)

infile <- "results/local_adaptation_scenarios/quality_control.csv.gz"
dat <- read_csv(infile) |>
  mutate(
    scenario = as.factor(scenario),
    QTLs = str_remove_all(QTLs, ".Rds") |> as.integer(),
  ) |>
  rename(Method = method)

plot_method <- function(method){
  plot_data <- dat |>
    filter(seed =="s104", Method == method, QTLs == 100)
  ylim <- plot_data$offset_causal |> range() |> round(digits = 2)
  ylim2 <- -log(plot_data$shifted_fitness) |> range() |> round(digits = 2)
  p1 <- ggplot(plot_data)+
    scale_x_continuous(expand=c(0,0), limits=c(-1000, 1000)) +
    geom_point(aes(x = offset_empirical, y = offset_causal), shape = 1, colour = "#88CCEE")+
    geom_smooth(
      aes(x = offset_empirical, y = offset_causal),
      colour = "black",
      method = "lm",
      fullrange=TRUE,
      se = TRUE,
      linetype = "dashed",
      formula = "y ~ x"
    ) +
    theme_classic()+
    theme(legend.position = "bottom")+
    xlim(ylim)+
    coord_cartesian(xlim=ylim, ylim = ylim)+
    facet_wrap(~scenario, ncol = 1)+
    ylab("Causal geometric genomic offset")+
    xlab("Empirical geometric genomic offset")

  p2 <- dat |>
    filter(seed =="s104", Method == method, QTLs == 100) |>
    ggplot()+
    scale_x_continuous(expand=c(0,0), limits=c(-1000, 1000)) +
    geom_point(aes(x = offset_empirical, y = -log(shifted_fitness)), shape = 1, colour = "#CC6677")+
    geom_smooth(
      aes(x = offset_empirical, y = -log(shifted_fitness)),
      colour = "black",
      method = "lm",
      se = TRUE,
      linetype = "dashed",
      fullrange = TRUE,
      formula = "y ~ x"
    ) +
    theme_classic()+
    theme(legend.position = "bottom")+
    xlim(ylim)+
    coord_cartesian(xlim=ylim, ylim = ylim2)+
    facet_wrap(~scenario, ncol = 1)+
    xlab("Empirical geometric genomic offset")+
    ylab("Negative log shifted fitness")
    
  ggarrange(
      p1, p2,
      ncol=2
    )
}

ggsave(
  "plots/01-loca-versus-non-local/geometric_go.pdf",
  plot_method("Geometric GO"), device = cairo_pdf,
  width = fig.witdh, height = fig.height, units = "mm", dpi = "retina"
)

ggsave(
  "plots/01-loca-versus-non-local/RONA_go.pdf",
  plot_method("RONA"), device = cairo_pdf,
  width = fig.witdh, height = fig.height, units = "mm", dpi = "retina"
)

ggsave(
  "plots/01-loca-versus-non-local/RDA_go.pdf",
  plot_method("GO RDA"), device = cairo_pdf,
  width = fig.witdh, height = fig.height, units = "mm", dpi = "retina"
)

ggsave(
  "plots/01-loca-versus-non-local/GF_go.pdf",
  plot_method("GF GO"), device = cairo_pdf,
  width = fig.witdh, height = fig.height, units = "mm", dpi = "retina"
)

# Histogram

hist_dat <- dat |>
    filter(Method == "Geometric GO") |> 
    pivot_longer(starts_with("offset")) |>
    filter(QTLs == 100, seed == "s100") |>
    filter(name %in% c("offset_causal", "offset_empirical")) |>
    mutate(
      name = factor(name,
        levels = c("offset_causal", "offset_empirical"),
        labels = c("(A) Causal", "(B) Empirical")
      )) |>
    # Invert level order of scenario
    mutate(scenario = fct_rev(scenario))
hist_plot <- hist_dat |>
  ggplot(aes(colour=name, y=value, x=scenario)) + 
  geom_violin() +
  geom_jitter(alpha = 0.5, height = 0)+
  scale_fill_aaas(name = "Scenario")+
  scale_colour_aaas(name = "Scenario")+
  coord_flip()+
  labs(x = "", y ="Geometric genomic offset", colour = "Scenario")+
  theme_classic()+
  theme(legend.position = "bottom")+
  facet_wrap(~name)
hist_plot

ggsave(
  "plots/01-loca-versus-non-local/hist_offsets.pdf",
  hist_plot, device = cairo_pdf,
  width = fig.witdh, height = fig.height, units = "mm", dpi = "retina"
  )


adj.r.squared <- dat |>
  group_split(scenario, seed, QTLs, Method) |>
  map(\(data){
    model <- lm(formula = -log(shifted_fitness) ~ offset_empirical, data = data)
    summary(model)$adj.r.squared
  }) |>
  as.numeric()

res <- dat |>
  group_by(scenario, seed, QTLs, Method) |>
  summarize(n = n()) |>
  ungroup() |>
  mutate(adj.r.squared = adj.r.squared) |>
  mutate(Radj = Method |>
    str_remove_all("GO") |>
    str_squish() |>
    str_wrap(5),
    Radj = ifelse(Radj == "GF", "Gradient forest", Radj)
    )

ggsig <- res |>
  ggplot(aes(x = Radj, y = adj.r.squared, colour = scenario)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, shape = 1, aes(colour = scenario))+
  ylim(c(-0.02, 1))+
  xlab("Genomic offset")+
  ylab("Adjusted R squared")+
  scale_color_brewer(name = "Scenario", palette="Dark2")+
  theme_classic()+
  theme(legend.position = "bottom")
  
ggsave(
  "plots/01-loca-versus-non-local/ggsignificant.pdf",
  ggsig, device = cairo_pdf,
  width = fig.witdh, height = fig.height, units = "mm", dpi = "retina"
)





