source("plots/theme.R")

library(ggsci)
library(marginaleffects)
set.seed(555)
infile <- "results/msprime/thyme_past.csv.gz"
dat <- read_csv(infile) |>
  mutate(
    n = file |> str_extract("n[0-9]+") |> str_remove("n") |> as.numeric(),
    r2adj = replace_na(r2adj, 0),
    speed = case_match(
      n,
      100 ~ "Slow (50 pseudogenerations)",
      20 ~ "Intermediate (10 pseudogenerations)",
      5 ~ "Fast (2.5 pseudogenerations)",
    )
  )


# Fast scenario

plot_current <- function(dat){
  datplot <- dat |>
    filter(dataset == "test", category %in% c("past_env")) |>
    mutate(
      type = if_else(type == "causal", "Causal geometric genomic offset", "Empirical geometric genomic offset"),
      rowid = row_number(),
      n = n * fraction / 2
    )
  test <- datplot |>
    slice_sample(prop = 0.3)
  train <- datplot |>
    anti_join(test, by = "rowid")
  mod <- lm(r2adj ~ n*type, data = train)
  p <- predictions(mod, conf_level = .95) |> 
    inferences(
      R = 5,
      method = "conformal_cv+",
      conformal_test = test)
  print(mean(p$r2adj <= p$pred.high & p$r2adj >= p$pred.low) > .90)
  ggplot(p, aes(x = n, y = r2adj, colour = type)) +
    geom_ribbon(
      aes(n, ymin = pred.low, ymax = pred.high),linetype = "dashed", alpha = 0.0)+
    geom_smooth(aes(n, r2adj), method = "lm", se = TRUE) +
    geom_point(alpha = .3)+
    theme_classic()+
    theme(legend.position = "bottom")+
    scale_colour_cosmic(name = "")+
    xlab("Sampling timepoint (pseudogenerations since climate change started)")+
    ylab("Adjusted R squared")+
    ylim(c(-0.1, 0.2))+
    theme(plot.title = element_text(hjust = 0.5))
}

fast <- dat |>
    filter(n == 5) |>  
  plot_current() +
  ggtitle("Fast rate\n(2.5 pseudogenerations elapsed)")
middle <- dat |>
  filter(n == 20) |>  
  plot_current() +
  ggtitle("Intermediate rate\n(10 pseudogenerations elapsed)")
slow <- dat |>
  filter(n == 100) |>  
  plot_current() +
  ggtitle("Slow rate\n(50 pseudogenerations elapsed)")

library(patchwork)
final <- fast + middle + slow + plot_layout(guides = "collect", axis_titles = "collect") & theme_classic() &theme(legend.position = "bottom")
final

ggsave(
  "plots/06-thyme/current.pdf",
  final,
  device = cairo_pdf, scale = 1.5,
  width = fig.witdh, height = fig.height, units = "mm", dpi = "retina"
)


## Predict future fitness loss!!!

plot_future <- function(dat){
  datplot <- filter(dat, dataset == "test", category %in% c("past_env_future")) |>
    mutate(
      type = if_else(type == "causal", "Causal geometric genomic offset", "Empirical geometric genomic offset"),
      rowid = row_number(),
      n = n * fraction / 2
    )
  test <- datplot |>
    slice_sample(prop = 0.3)
  train <- datplot |>
    anti_join(test, by = "rowid")
  mod <- lm(r2adj ~ n*type, data = train)
  p <- predictions(mod, conf_level = .95) |> 
    inferences(
      R = 5,
      method = "conformal_cv+",
      conformal_test = test)
  print(mean(p$r2adj <= p$pred.high & p$r2adj >= p$pred.low) > .90)
  ggplot(p, aes(x = n, y = r2adj, colour = type)) +
    geom_point(alpha = .3)+
    geom_ribbon(
      aes(n, ymin = pred.low, ymax = pred.high),linetype = "dashed", alpha = 0)+
    geom_smooth(aes(n, r2adj),method = "lm", se = TRUE) +
    theme_classic()+
    theme(legend.position = "bottom")+
    xlab("Sampling timepoint (pseudogenerations since climate change started)")+
    ylab("Adjusted R squared")+
    scale_colour_cosmic(name = "")+
    ylim(c(-0.1, 0.3))+
    theme(plot.title = element_text(hjust = 0.5))
}

fast <- dat |>
  filter(n == 5) |>  
  plot_future() +
  ggtitle("Fast rate\n(2.5 pseudogenerations elapsed)")
middle <- dat |>
  filter(n == 20) |>  
  plot_future() +
  ggtitle("Intermediate rate\n(10 pseudogenerations elapsed)")

slow <- dat |>
  filter(n == 100) |>  
  plot_future() +
  ggtitle("Slow rate\n(50 pseudogenerations elapsed)")

final <- fast + middle + slow + plot_layout(guides = "collect" , axis_titles = "collect")& theme_classic() & theme(legend.position = "bottom")
final

ggsave(
  "plots/06-thyme/future.pdf",
  final,
  device = cairo_pdf, scale = 1.5,
  width = fig.witdh, height = fig.height, units = "mm", dpi = "retina"
)

dat |>
  filter(
    fraction == 0.0, category == "past_env", type == "empirical", 
    dataset %in% c("test_ecotypeA", "test_ecotypeB")
      ) |>
  mutate(
    n = n * fraction / 2,
    dataset = case_match(
      dataset,
      "test_ecotypeA" ~ "Drough-tolerant",
      "test_ecotypeB" ~ "Freezing-tolerant",
    ),
      type = if_else(type == "causal", "Causal geometric genomic offset", "Empirical geometric genomic offset")
  ) |>
  ggplot(aes(x = r2adj, fill = dataset)) +
  theme_classic()+
  geom_histogram(colour = "black")+
  theme(legend.position = "bottom")+
  scale_fill_manual(name = "Ecotype", values = c("firebrick", "lightblue"))+
  xlab("Adjusted R squared")+
  ylab("Frequency")

ggsave(
  "plots/06-thyme/ecotypes.pdf",
  last_plot(),
  device = cairo_pdf,
  width = fig.witdh, height = fig.height, units = "mm", dpi = "retina"
)







