library(ggplot2)
library(tidyr)
alphadiv2 <- read.csv("C:/Users/Lindsay Alma/Dropbox/Eelgrass/West Coast TagSeq primers/FullDataset/alphadiv_823.csv", header = TRUE)
head(alphadiv2)

####################################STATE################################################
# Custom color palette for each state
box_colors <- c(    
  "AK" = "orange",
  "BB" = "gold",
  "BC" = "darkgreen",
  "OR" = "blue",
  "SD" = "lightblue",
  "WA" = "purple"
)

# Reorder the State factor levels
alphadiv2$State <- factor(alphadiv2$State, levels = c("AK", "BC", "WA", "OR", "BB", "SD"))



# Reshape the data to long format for easier plotting
alphadiv_long <- alphadiv2 %>%
  pivot_longer(cols = c(observed, shannon), names_to = "Diversity_Metric", values_to = "Value")

# Plot boxplots with custom colors, reordered states, and updated title
ggplot(alphadiv_long, aes(x = State, y = Value, fill = State)) +
  geom_boxplot() +
  facet_wrap(~Diversity_Metric, scales = "free_y", 
             labeller = as_labeller(c("observed" = "Observed", "shannon" = "Shannon"))) +  # Capitalize facet labels
  scale_fill_manual(values = box_colors) +  # Apply custom colors
  theme_classic() +  # Removes grid
  labs(y = "Diversity Value", x = "State") +
  theme(legend.position = "none",
        strip.text = element_text(size = 14),
        strip.background = element_blank(),
        axis.line = element_line(color = "black"))  # Adds axis lines




###############################Sample_Type############################

# Custom color palette for Sample_Type
sample_colors <- c("L" = "orange4", "G" = "green")


# Plot boxplots with Sample_Type instead of State
ggplot(alphadiv_long, aes(x = Sample_Type, y = Value, fill = Sample_Type)) +
  geom_boxplot() +
  facet_wrap(~Diversity_Metric, scales = "free_y", 
             labeller = as_labeller(c("observed" = "Observed", "shannon" = "Shannon"))) +  
  scale_fill_manual(values = c("G" = "green", "L" = "orange4")) +  
  scale_x_discrete(labels = c("G" = "Asymptomatic", "L" = "Symptomatic")) +  # Change x-axis labels
  theme_classic() +  
  labs(y = "Diversity Value", x = "Sample Type") +
  theme(legend.position = "none",
        strip.text = element_text(size = 14),
        strip.background = element_blank(),
        axis.line = element_line(color = "black"))



#############################LEsion_type#############
# Filter out "Likely Crescent" and "Possibly Crescent" from the dataset
alphadiv_long_filtered <- alphadiv_long %>%
  filter(!Lesion_Type %in% c("Likely Crescent", "Possibly Crescent"))

# Custom color palette for Lesion_Type
lesion_colors <- c("Crescent" = "cornsilk4", "Not crescent" = "forestgreen")

# Plot boxplots with Lesion_Type
ggplot(alphadiv_long_filtered, aes(x = Lesion_Type, y = Value, fill = Lesion_Type)) +
  geom_boxplot() +
  facet_wrap(~Diversity_Metric, scales = "free_y", 
             labeller = as_labeller(c("observed" = "Observed", "shannon" = "Shannon"))) +  
  scale_fill_manual(values = c("Crescent" = "cornsilk4", "Not crescent" = "forestgreen")) +  
  scale_x_discrete(labels = c("Crescent" = "Crescent", "Not crescent" = "Not Crescent")) +  # Change x-axis labels
  theme_classic() +  
  labs(y = "Diversity Value", x = "Lesion Type") +
  theme(legend.position = "none",
        strip.text = element_text(size = 14),
        strip.background = element_blank(),
        axis.line = element_line(color = "black"))
