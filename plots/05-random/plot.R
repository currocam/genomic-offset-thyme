source("plots/theme.R")

res <- read_csv("results/local_adaptation_scenarios/spurious_genomic_offsets.csv.gz")

library(ggsci)
res |>
  filter(!is.na(seed)) |>
  pivot_longer(-c(seed, L, n_uncorrelated)) |>
  mutate(
    label = case_when(
      name == "with_testing" & is.na(value) ~ "No candidate found after testing",
      name == "with_testing" & !is.na(value) ~ "Spurious estimates with\nputatively adaptative loci",
      name == "without_testing" ~ "Spurious estimates with\n10% randomly selected loci"
    ),
    value = replace_na(value, 0)
    ) |>
  ggplot(aes(x = value, fill = label))+
  geom_histogram(colour = "black", bins = 30)+
  theme_classic()+
  geom_vline(xintercept = 0, linetype = "dashed")+
  scale_fill_tron(name = "Genomic offset")+
  theme(legend.position = "bottom",axis.ticks.y = element_blank())+
  xlab("Adjusted R squared")+
  ylab("Frequency")

ggsave(
  "plots/05-random/random.pdf",
  device = cairo_pdf,
  width = fig.witdh, height = fig.height, units = "mm", dpi = "retina"
)

