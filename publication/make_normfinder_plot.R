#!/usr/bin/env Rscript

library(ggplot2)

data <- read.csv("stabilities.csv")

data[["gene_id"]] <- factor(data[["gene_id"]], levels = unique(data[["gene_id"]]))
data <- data[order(data[["normfinder_stability_value"]], decreasing = FALSE), ]

ggplot(data, aes(x = gene_id, y = normfinder_stability_value)) +
  geom_bar(stat = "identity", fill = "deeppink4") +
  labs(x = "Gene", y = "NormFinder stability value") +
  theme_minimal() +
  theme(
    text=element_text(size=16,  family="Liberation Sans"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("normfinder_stability_values.png", width = 8, height = 6, dpi = 300) 
