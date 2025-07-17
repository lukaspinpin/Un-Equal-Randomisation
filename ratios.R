library(ggplot2)
library(Rlab) # This library is not used in the provided code, but kept for completeness
library(patchwork)

p0 <- c(0.01, 0.1, 0.3, 0.5)
p1 <- seq(0, 1, 0.005)
plot_list <- list()

for (i in 1:length(p0)) {
  p_CR <- rep(0.5, length(p1))
  p_N <- c()
  p_R <- c()
  p_RG <- c()
  
  for (j in 1:length(p1)) {
    p_N[j] <- sqrt(p1[j] * (1 - p1[j])) /
      (sqrt(p0[i] * (1 - p0[i])) + sqrt(p1[j] * (1 - p1[j])))
    p_R[j] <- sqrt(p1[j]) / (sqrt(p1[j]) + sqrt(p0[i]))
    p_RG[j] <- p1[j] / (p1[j] + p0[i])
  }
  
  data <- data.frame(p1 = p1, p_CR = p_CR, p_N = p_N, p_R = p_R, p_RG = p_RG)
  
  # Title using bquote and curly brace simulation
  if (i < 4) {
    plot_title <- bquote(p[0] == "{" * .(p0[i]) * "," ~ .(1 - p0[i]) * "}")
  } else {
    plot_title <- bquote(p[0] == .(p0[i]))
  }
  
  plot <- ggplot(data, aes(x = p1)) +
    geom_line(aes(y = p_CR, color = "p_CR"), linetype = "solid", size = 1.75) +
    geom_line(aes(y = p_N, color = "p_N"), linetype = "dashed", size = 1.75) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      x = expression(p[1]),
      y = expression(rho),
      title = plot_title
    ) +
    scale_color_manual(
      name = "Allocation",
      values = c("p_CR" = "lightgreen", "p_N" = "lightblue"),
      labels = c(expression(rho[ER]), expression(rho[N[1]]))
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      text = element_text(size = 20),
      legend.position = "bottom",
      legend.text.align = 0,
      legend.key.width = unit(2, "cm"),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 20),
      # Add this line to control legend text size
      legend.text = element_text(size = 26) # Adjust size as needed
    )
  
  plot_list[[i]] <- plot
}

final_plot <- (plot_list[[1]] + plot_list[[2]]) / (plot_list[[3]] + plot_list[[4]]) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom",
                                          legend.text = element_text(size = 18)) # Also add here for consistency

ggsave("combined_figure_binary.pdf", final_plot, width = 12, height = 12, dpi = 300)