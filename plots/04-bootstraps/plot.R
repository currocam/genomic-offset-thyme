source("plots/theme.R")
library(ggsci)

infile <-"results/local_adaptation_scenarios/quality_control_bootstrapping.csv.gz"
data <- read_csv(infile) |>
  filter(FDR == 0.05) |>
  mutate(
    method = case_when(
      method == "GeometricGO" ~ "Geometric",
      method == "GradientForestGO" ~ "Gradient Forest",
      method == "RDAGO" ~ "RDA",
      method == "RONA" ~ "RONA",
    )
  )

p1 <- data |>
  mutate(file = fct_reorder(file, causal)) |>
  ggplot()+
  geom_errorbar(aes(y = file, xmin = minconf95,xmax = maxconf95))+
  geom_point(aes(y = file, x = causal,colour = "Causal"), size = 1.5)+
  geom_point(aes(y = file, x = empirical,colour = "Empirical"), size = 1.5)+
  theme_classic()+
  theme(
    legend.position = "bottom",axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
        )+
  scale_color_aaas(name = "Genomic offset")+
  xlab("Spearman correlation 95% confidence interval ")+
  ylab("Random seeds")+
  xlim(c(0, 1))+
  facet_wrap(~method)
p1

ggsave(
  "plots/04-bootstraps/error_bars.pdf",
  p1, device = cairo_pdf,
  width = fig.witdh, height = fig.height, units = "mm", dpi = "retina"
)

data |>
  ggplot(aes(y = method, x = bootruntime))+
  geom_boxplot()+
  xlab("100 bootstrap iterations runtime (seconds) ")+
  ylab("")+
  scale_x_log10()+
  theme_classic()

ggsave(
  "plots/04-bootstraps/time.pdf",
  last_plot(), device = cairo_pdf,
  width = fig.witdh, height = fig.height, units = "mm", dpi = "retina"
)

infile <- "results/local_adaptation_scenarios/regression_bootstraps.csv.gz"
data <- read_csv(infile) |>
  filter(file == "steps/slim/m3.2_s102_nQTL1s100_nQTL2s_100.Rds")

# First, we do kmeans with x and y to create populations
set.seed(123)
prepare_data <- function(data){
  km.res <- kmeans(data[, c("x", "y")], 11^2, nstart = 25)
  data |>
    mutate(
      cluster = km.res$cluster,
      method = case_when(
        method == "GeometricGO" ~ "Geometric",
        method == "GradientForestGO" ~ "Gradient Forest",
        method == "RDAGO" ~ "RDA",
        method == "RONA" ~ "RONA",
      ) 
    ) |>
    mutate(across(starts_with("boot"), rank)) |>
    group_by(cluster, method) |>
    summarise(
      across(c(
        starts_with("boot"), empirical, causal,
        negative_log_fitness, ), mean)
    ) |>
    ungroup() |>
    mutate(
      across(starts_with("boot"), rank),
      empirical = rank(empirical),
      causal = rank(causal),
      negative_log_fitness = rank(negative_log_fitness),
    ) |>
    pivot_longer(starts_with("boot"), names_to = "boot_number", values_to = "rank") |>
    group_by(method, empirical, causal, negative_log_fitness, cluster) |>
    summarise(
      minconf95 = quantile(rank, 0.025),
      maxconf95 = quantile(rank, 0.975),
      median = median(rank)
    )
}

datplot <- data |>
  group_split(method) |>
  map_dfr(prepare_data)



p2 <- datplot|>
    ggplot() +
    geom_abline(linetype = "dashed")+
    geom_point(aes(y = causal, x = empirical, colour = method))+
    geom_errorbar(aes(y = causal, xmin = minconf95, xmax = maxconf95, colour = method),width = 1)+
    xlab("95% confidence interval ranked \npopulation mean genomic offset")+
    ylab("Rank of population mean \n causal genomic offset")+
    theme(legend.position = "bottom")+
    scale_color_brewer(name = "Method", palette="Dark2")+
    facet_wrap(~method)
  
p3 <- datplot|>
    ggplot() +
    geom_abline(linetype = "dashed")+
    geom_point(aes(y = negative_log_fitness, x = empirical, colour = method))+
    geom_errorbar(aes(y = negative_log_fitness, xmin = minconf95, xmax = maxconf95,  colour = method), width = 1)+
    xlab("95% confidence interval ranked \npopulation mean genomic offset")+
    ylab("Rank of population mean \nnegative logarithm of shifted fitness")+
    theme(legend.position = "none")+
    scale_color_brewer(name = "Method", palette="Dark2")+
    facet_wrap(~method)

# Combine them with shared legend in the bottom
library(patchwork)
p2 + p3 + plot_layout(guides='collect') &
  theme_classic() & theme(legend.position = "none")
ggsave(
  "plots/04-bootstraps/rank_custom.pdf",
  last_plot(), device = cairo_pdf,
  width = fig.witdh, height = fig.height, units = "mm", dpi = "retina"
)
