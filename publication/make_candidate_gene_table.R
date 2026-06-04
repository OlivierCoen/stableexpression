library(flextable)

data <- 

ft <- flextable(data) %>%
  set_caption("Gene Expression Analysis") %>%
  bold(i = 1) %>%  # Bold header
  color(i = ~ p_value < 0.05, j = 3, color = "red") %>%  # Color p-values
  add_header_row(values = c("", "Values"), colwidths = c(1, 2)) %>%
  theme_booktabs()  # Professional styling

ft