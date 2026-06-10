library(flextable)
library(dplyr)

set_flextable_defaults(font.family = "Liberation Sans")

data <- read.csv("candidate_genes.csv")

usual_genes <- c("Act5C", "Gapdh1", "RpL32")

data <- data[
  order(data$stability_score, decreasing = FALSE),
  c("gene_id", "gene_name", "section", "rank", "stability_score", "normfinder_score", "genorm_score", "cv", "rcvm")
] 

print(data)

data <- data %>%
  mutate(
    stability_score = round(stability_score, 3),
    normfinder_score = round(normfinder_score, 3),
    genorm_score = round(genorm_score, 3),
    cv = round(cv, 3),
    rcvm = round(rcvm, 3)
  ) %>%
  rename(
    "Stability score" = stability_score,
    "Gene ID" = gene_id,
    "Gene name" = gene_name,
    "Section" = section,
    "Rank" = rank,
    "NormFinder score" = normfinder_score,
    "GeNorm score" = genorm_score,
    "CV" = cv,
    "RCVm" = rcvm
  )


ft <- flextable(data) %>%
  theme_booktabs() %>%
  bold(i = ~ `Gene name` %in% usual_genes)

ft

save_as_image(ft, path = "candidate_gene_table.png")