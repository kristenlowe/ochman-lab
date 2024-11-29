library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(gridExtra)
library(grid)

setwd("/stor/scratch/Ochman/kristen/pangenome")

# import data
lifestyles <- read_csv("57_pangenome_lifestyles.csv")
species_pan <- read_csv("bac_57_difference_levels.csv")

life_df <- inner_join(x = species_pan, y = lifestyles, 
                 by = c("species" = "Species"))

properties <- c("length", "GC", "GC_3rd", "polarAA", "hydrophobicAA", "metabol", 
              "cai", "instability", "disorder", "aa_comp_bias", "transmembrane", "tm_coverage")

# Reshape the data
life_df_long <- life_df %>%
  filter(property == "cai") %>%
  pivot_longer(cols = c(core_value, cloud_value), 
               names_to = "Group", 
               values_to = "Value") 

# Create the boxplot
ggplot(life_df_long, aes(x = Intra_or_extracellular, y = Value, fill = Group)) +
  geom_boxplot() +
  labs(x = "Intra or Extracellular", y = "Value", fill = "Group", title = "cai") +
  #ylim(0, 1000) +
  theme_minimal()




# Define the properties
properties <- c("length", "GC", "GC_3rd", "polarAA", "hydrophobicAA", "metabol", 
                "cai", "instability", "disorder", "aa_comp_bias", "transmembrane", "tm_coverage")

# Function to filter out outliers
filter_outliers <- function(data, value_col) {
  lower <- quantile(data[[value_col]], 0.05, na.rm = TRUE)
  upper <- quantile(data[[value_col]], 0.95, na.rm = TRUE)
  data %>% filter(!!sym(value_col) >= lower, !!sym(value_col) <= upper)
}

# HOST LOCATION
intra_or_extracellular_plots <- map(properties, function(prop) {
  # Reshape and filter data for the current property
  life_df_long <- life_df %>%
    filter(property == prop) %>%
    pivot_longer(cols = c(core_value, cloud_value), 
                 names_to = "Group", 
                 values_to = "Value") %>%
    group_by(Group) %>%
    group_modify(~ filter_outliers(.x, "Value")) %>% # Remove outliers
    
    ungroup()
  
  # Create the boxplot for the current property
  ggplot(life_df_long, aes(x = Intra_or_extracellular, y = Value, fill = Group)) +
    geom_boxplot(outlier.shape = NA) + # Suppress outliers from being plotted
    labs(x = "", 
         y = "Value", 
         fill = "Group", 
         title = prop) +
    theme_minimal()
})

# Combine all plots into a single image
grid.arrange(
  arrangeGrob(grobs = intra_or_extracellular_plots, ncol = 3), # The grid of plots
  top = textGrob("Host Location (Intra or Extracellular)", gp = gpar(fontsize = 16, fontface = "bold"))
)

#combined_plot <- grid.arrange(grobs = intra_or_extracellular_plots, ncol = 3)

# Save the combined plot
#ggsave("combined_boxplots.png", combined_plot, width = 15, height = 10)


# HOST OR FREE LIVING
# Create a list of plots
host_or_free_living_plots <- map(properties, function(prop) {
  # Reshape and filter data for the current property
  life_df_long <- life_df %>%
    filter(property == prop) %>%
    pivot_longer(cols = c(core_value, cloud_value), 
                 names_to = "Group", 
                 values_to = "Value") %>%
    group_by(Group) %>%
    group_modify(~ filter_outliers(.x, "Value")) %>% # Remove outliers
    
    ungroup()
  
  # Create the boxplot for the current property
  ggplot(life_df_long, aes(x = Host_or_free, y = Value, fill = Group)) +
    geom_boxplot(outlier.shape = NA) + # Suppress outliers from being plotted
    labs(x = "", 
         y = "Value", 
         fill = "Group", 
         title = prop) +
    theme_minimal()
})

# Combine all plots into a single image
grid.arrange(
  arrangeGrob(grobs = host_or_free_living_plots, ncol = 3), # The grid of plots
  top = textGrob("Host or Free Living", gp = gpar(fontsize = 16, fontface = "bold"))
)


# HOST RELIANCE
# Create a list of plots
host_reliance_plots <- map(properties, function(prop) {
  # Reshape and filter data for the current property
  life_df_long <- life_df %>%
    filter(property == prop) %>%
    pivot_longer(cols = c(core_value, cloud_value), 
                 names_to = "Group", 
                 values_to = "Value") %>%
    group_by(Group) %>%
    group_modify(~ filter_outliers(.x, "Value")) %>% # Remove outliers
    
    ungroup()
  
  # Create the boxplot for the current property
  ggplot(life_df_long, aes(x = Obligate_facultative, y = Value, fill = Group)) +
    geom_boxplot(outlier.shape = NA) + # Suppress outliers from being plotted
    labs(x = "", 
         y = "Value", 
         fill = "Group", 
         title = prop) +
    theme_minimal()
})

# Combine all plots into a single image
grid.arrange(
  arrangeGrob(grobs = host_reliance_plots, ncol = 3), # The grid of plots
  top = textGrob("Host Reliance", gp = gpar(fontsize = 16, fontface = "bold"))
)



# MOTILITY
# Create a list of plots
motility_plots <- map(properties, function(prop) {
  # Reshape and filter data for the current property
  life_df_long <- life_df %>%
    filter(property == prop) %>%
    pivot_longer(cols = c(core_value, cloud_value), 
                 names_to = "Group", 
                 values_to = "Value") %>%
    group_by(Group) %>%
    group_modify(~ filter_outliers(.x, "Value")) %>% # Remove outliers
    
    ungroup()
  
  # Create the boxplot for the current property
  ggplot(life_df_long, aes(x = Motility, y = Value, fill = Group)) +
    geom_boxplot(outlier.shape = NA) + # Suppress outliers from being plotted
    labs(x = "", 
         y = "Value", 
         fill = "Group", 
         title = prop) +
    theme_minimal()
})

# Combine all plots into a single image
grid.arrange(
  arrangeGrob(grobs = motility_plots, ncol = 3), # The grid of plots
  top = textGrob("Motility", gp = gpar(fontsize = 16, fontface = "bold"))
)



# EFFECT ON HOST
# Create a list of plots
effect_on_host_plots <- map(properties, function(prop) {
  # Reshape and filter data for the current property
  life_df_long <- life_df %>%
    filter(property == prop) %>%
    pivot_longer(cols = c(core_value, cloud_value), 
                 names_to = "Group", 
                 values_to = "Value") %>%
    group_by(Group) %>%
    group_modify(~ filter_outliers(.x, "Value")) %>% # Remove outliers
    
    ungroup()
  
  # Create the boxplot for the current property
  ggplot(life_df_long, aes(x = Effect_on_host, y = Value, fill = Group)) +
    geom_boxplot(outlier.shape = NA) + # Suppress outliers from being plotted
    labs(x = "", 
         y = "Value", 
         fill = "Group", 
         title = prop) +
    theme_minimal()
})

# Combine all plots into a single image
grid.arrange(
  arrangeGrob(grobs = effect_on_host_plots, ncol = 3), # The grid of plots
  top = textGrob("Effect on Host", gp = gpar(fontsize = 16, fontface = "bold"))
)


# HOST TYPE
# Create a list of plots
host_type_plots <- map(properties, function(prop) {
  # Reshape and filter data for the current property
  life_df_long <- life_df %>%
    filter(property == prop) %>%
    pivot_longer(cols = c(core_value, cloud_value), 
                 names_to = "Group", 
                 values_to = "Value") %>%
    group_by(Group) %>%
    group_modify(~ filter_outliers(.x, "Value")) %>% # Remove outliers
    
    ungroup()
  
  # Create the boxplot for the current property
  ggplot(life_df_long, aes(x = Host_type, y = Value, fill = Group)) +
    geom_boxplot(outlier.shape = NA) + # Suppress outliers from being plotted
    labs(x = "", 
         y = "Value", 
         fill = "Group", 
         title = prop) +
    theme_minimal()
})

# Combine all plots into a single image
grid.arrange(
  arrangeGrob(grobs = host_type_plots, ncol = 3), # The grid of plots
  top = textGrob("Host Type", gp = gpar(fontsize = 16, fontface = "bold"))
)
