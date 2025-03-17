# Load necessary libraries
library(ggplot2)
library(dplyr)
library(rstatix)

create_significance_plot <- function(data, y_variable, y_axis_label) {
  # Perform ANOVA for each treat_day
  anova_results <- data %>%
    group_by(treat_day) %>%
    anova_test(as.formula(paste(y_variable, "~ treatement"))) %>%
    ungroup()
  
  # Perform Tukey's post-hoc test (which follows an ANOVA internally)
  posthoc_results <- data %>%
    group_by(treat_day) %>%
    tukey_hsd(as.formula(paste(y_variable, "~ treatement"))) %>%
    ungroup() %>%
    filter(p.adj < 0.05)  # Keep only significant comparisons
  
  # Calculate y.positions for plotting the significance bars if any significant comparisons exist
  if (nrow(posthoc_results) > 0) {
    posthoc_results <- posthoc_results %>%
      group_by(treat_day) %>%
      mutate(y.position = max(data[[y_variable]], na.rm = TRUE) +
               seq(0.1, by = 0.1, length.out = n())) %>%
      ungroup()
  }
  
  # Create the base plot
  bxp <- ggplot(data, aes(x = treatement, y = .data[[y_variable]], group = treatement, color = treatement)) +
    stat_summary(fun = "mean", geom = "pointrange",
                 fun.max = function(x) mean(x) + sd(x),
                 fun.min = function(x) mean(x) - sd(x),
                 fatten = 1.8) +
    geom_point(shape = 2, position = position_jitter(width = 0.07, seed = 84)) +
    scale_y_continuous(name = y_axis_label, breaks = scales::breaks_pretty(10)) +
    scale_colour_brewer(palette = "Dark2") +
    theme_classic() +
    theme(axis.text = element_text(size = 9),
          legend.position = "none",
          panel.border = element_rect(colour = "black", fill = NA)) +
    facet_wrap(~treat_day)
  
  # Add significance bars if any significant post-hoc comparisons exist
  if (nrow(posthoc_results) > 0) {
    bxp <- bxp +
      stat_pvalue_manual(posthoc_results, label = "p.adj.signif", tip.length = 0.01)
  } else {
    message("No significant differences found between groups.")
  }
  
  # Return a list containing the plot, the ANOVA table, and the post-hoc results table
  return(list(
    plot = bxp,
    anova = anova_results,
    posthoc = posthoc_results
  ))
}
