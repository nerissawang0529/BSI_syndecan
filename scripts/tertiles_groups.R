library(ggplot2)
library(dplyr)

merged_data <- read.csv("original_data/merged_data.csv")

# Example dataset
set.seed(123)
data <- data.frame(value = merged_data$syndecan_1)

# Create tertiles
data <- data %>%
  mutate(tertile = ntile(value, 3))  # ntile() divides into 3 groups

# Plot histogram colored by tertile
ggplot(data, aes(x = value, fill = factor(tertile))) +
  
  geom_histogram( position = "identity", alpha = 0.6) +
  scale_x_continuous(breaks = c(seq(6.5, 9.169963, by = 0.1),
                                seq(9.169963, 9.718656, by = 0.1),
                                seq( 9.718656,12, by = 0.1))) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) + # Custom colors
  labs(fill = "Tertile") +
  theme_minimal()


head(data)

max(data$value[data$tertile==1])
min(data$value[data$tertile==2])

tertile_cuts <- unname(quantile(merged_data$syndecan_1, c(0.333, 0.667)))

ggplot(data, aes(x = value, fill = factor(tertile))) +
  geom_histogram(breaks = c(rev(9.169963 - (1:12 * 0.2743468)),
                            tertile_cuts[1], sum(tertile_cuts)/2 , tertile_cuts[2],
                            (1:10 * 0.2743468) + 9.718656), alpha = 0.9,
                            colour="black") + 
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  labs(fill = "Tertile", x = "Value", y = "Count") +
  theme_minimal()

