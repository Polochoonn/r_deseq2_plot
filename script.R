data_mutations = read.csv("./comp_by_mutation.csv")

install.packages("ggpubr")

library("ggpubr")
library("dplyr")
library("ggplot2")
library("rstatix")
library(purrr)
library(tidyr)
library(emmeans)
library(tibble)





# Specify the mutations for comparison
mutations_to_compare <- c("EGFR", "KRAS", "none", "STK11", "TP53")

# Create an empty list to store results for each comparison
results_list <- list()

for (comp in unique(data_mutations$compartment)) {
  compartment_data <- data_mutations %>% filter(compartment == comp)
  
  # Loop through each mutation for comparison
  for (mut in mutations_to_compare) {
    # Create a subset for the mutation of interest
    mut_subset <- compartment_data %>% filter(single_mutation == mut)
    
    # Create a subset for the other mutations
    other_mutations_subset <- compartment_data %>% filter(single_mutation != mut)
    
    # Perform a t-test
    t_test_result <- t.test(mut_subset$fraction, other_mutations_subset$fraction)
    
    # Store the results in the list
    results_list[[paste(comp, mut, sep = "_")]] <- t_test_result
  }
}

# Create a data frame to organize the results
results_df <- tibble::tibble(
  compartment = character(),
  single_mutation = character(),
  group2 = character(),
  comparison = character(),
  statistic = double(),
  p_value = double(),
  parameter = double(),
  method = character(),
  alternative = character(),
  conf_int_low = double(),
  conf_int_high = double(),
  p.value.signif = character()  # New column for significance annotation
)

# Fill the data frame with t-test results
for (comp in unique(data_mutations$compartment)) {
  for (mut in mutations_to_compare) {
    result_name <- paste(comp, mut, sep = "_")
    result <- results_list[[result_name]]
    
    # Extract relevant information from t-test result
    t_test_info <- tibble::tibble(
      compartment = comp,
      single_mutation = mut,
      group2 = "other",
      comparison = result_name,
      statistic = result$statistic,
      p_value = result$p.value,
      parameter = result$parameter,
      method = result$method,
      alternative = result$alternative,
      conf_int_low = result$conf.int[1],
      conf_int_high = result$conf.int[2],
      p.value.signif = ifelse(result$p.value < 0.001, "***",
                              ifelse(result$p.value < 0.01, "**",
                                     ifelse(result$p.value < 0.05, "*", "")))
    )
    
    # Append the information to the main results data frame
    results_df <- dplyr::bind_rows(results_df, t_test_info)
  }
}

# Print the results table
print(results_df)




############ Plot avec statistiquen en Y  ##########################

p_combined <- ggplot(results_df, aes(x = compartment, y = statistic, fill = single_mutation)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), alpha = 0.7) +
  geom_errorbar(aes(ymin = conf_int_low, ymax = conf_int_high),
                position = position_dodge(0.8), width = 0.25) +
  geom_text(aes(label = p.value.signif),
            position = position_dodge(width = 0.8), vjust = -0.5) +
  labs(title = "Bar Plot of Fraction by Compartment and Mutation",
       x = "Compartment", y = "Statistic") +
  theme_minimal()

print(p_combined)




############ Plot avec fraction en Y comme sur le tien mais il est degueu ##########################


summary_df <- merged_df %>%
  group_by(compartment, single_mutation) %>%
  summarise(
    p.value.signif = first(p.value.signif),
    fraction = first(fraction),  # Add this line to include the "fraction" variable
    .groups = "drop"
  )




# Bar plot with unique significance labels
p_combined <- ggplot(merged_df, aes(x = compartment, y = fraction, fill = single_mutation)) +
  geom_bar(stat = "identity", position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = conf_int_low, ymax = conf_int_high),
                position = position_dodge(0.8), width = 0.25) +
  geom_text(data = summary_df,
            aes(label = p.value.signif),
            position = position_dodge(width = 0.8), vjust = -0.5) +
  labs(title = "Bar Plot of Fraction by Compartment and Mutation",
       x = "compartment", y = "fraction") +
  theme_minimal()

print(p_combined)






############ Ton PLot ##########################

# Add a column indicating significance

# Bar plot with all groups






stat<- data_mutations %>%
  group_by(compartment) %>%
  t_test(fraction ~ single_mutation) %>%
  add_xy_position(fun = "mean_sd", x = "compartment", dodge = 0.8) %>% 
  filter(p.adj.signif != "ns")


significant_comparisons_stat <- stat %>%
  filter(p.adj < 0.05)  # You can adjust the significance level as needed





p_combined <- (ggbarplot(data_mutations, x = "compartment", y = "fraction", add = "mean_sd",
                         fill = "single_mutation",
                         position = position_dodge(0.8)
)   # Set the y-axis limits       
# + stat_compare_means(aes(group = compartment), label = "p.signif", label.y = 1)
+ stat_pvalue_manual( stat,  label = "p.adj.signif")  

# + stat_compare_means(comparisons=fraction ~ single_mutation)
)
p_combined



