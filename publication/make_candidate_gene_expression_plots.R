library(ggplot2)
library(dplyr)
library(arrow)
library(tidyr)
library(tibble)

counts <- arrow::read_parquet("all_counts.missing_values_imputed.parquet")

candidate_data <- read.csv("candidate_genes.csv")
candidate_gene_ids <- candidate_data$gene_id

# keep only candidate gene counts
counts <- counts |> filter(gene_id %in% candidate_gene_ids)

# recode gene_id to gene_name
gene_names <- candidate_data %>%
  select(gene_id, gene_name) %>%
  deframe()

counts <- counts %>%
  mutate(gene_id = recode(gene_id, !!!gene_names))
  

# pivot to long format and sort gene_id by desired order
desired_order <- candidate_data[order(candidate_data$section, candidate_data$rank), ]$gene_name
data_long <- counts %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "value") %>%
  mutate(gene_id = factor(gene_id, levels = rev(desired_order)))
    
print(data_long)

g <- ggplot(data_long, aes(x = value, y = gene_id)) +
  geom_boxplot() +
  theme_bw() + 
  labs(y = NULL, x = "Normalized gene expression") +
  theme(
    axis.text.y = element_text(size=12)
  )

g


