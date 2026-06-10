#!/usr/bin/env Rscript

library(ggplot2)

data <- read.csv("m_measures.csv")

data[["gene_id"]] <- factor(data[["gene_id"]], levels = unique(data[["gene_id"]]))
data <- data[order(data[["m_measure"]], decreasing = FALSE), ]
print(data)
ggplot(data, aes(x = gene_id, y = m_measure)) +
  geom_bar(stat = "identity", fill = "darkcyan") +
  labs(x = "Gene", y = "GeNorm M-measure") +
  theme_minimal() +
  theme(
    text=element_text(size=16,  family="Liberation Sans"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("genorm_m_measures.png", width = 8, height = 6, dpi = 300) 
