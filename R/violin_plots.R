library(ggplot2)

# Generate sample data
set.seed(123)
data <- data.frame(
  condition = rep(c("enfermo", "sano"), each = 40),
  cell_type = rep(c("tipo1", "tipo2"), times = 40),
  value = c(rnorm(40, 5, 1), rnorm(40, 7, 1.5), rnorm(40, 4, 1), rnorm(40, 6, 1.5))
)

# Plot violin
ggplot(data, aes(x = condition, y = value, fill = cell_type)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = c("orange", "lightblue")) +
  facet_wrap(~ condition, scales = "free_x") +
  theme_bw()


library(ggplot2)

# Example data
set.seed(123)
df <- data.frame(
  condition = factor(rep(c("healthy", "diseased"), each = 200)),
  cell_type = factor(rep(c("cell type 1", "cell type 2"), each = 100)),
  value = c(rnorm(100, 5, 2), rnorm(100, 10, 2), rnorm(100, 7, 2), rnorm(100, 12, 2))
)

# Create a summary data frame to calculate the statistics for each group
summary_df <- aggregate(value ~ condition + cell_type, data = df,
                        FUN = function(x) c(mean = mean(x), sd = sd(x)))

# Create the violin plot
ggplot(df, aes(x = condition, y = value, fill = cell_type)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_point(data = summary_df, aes(x = condition, y = value[, "mean"], group = cell_type),
             color = "white", size = 3, shape = 21, fill = "black") +
  geom_errorbar(data = summary_df, aes(x = condition, y = value[, "mean"], ymin = value[, "mean"] - value[, "sd"], ymax = value[, "mean"] + value[, "sd"], group = cell_type),
                color = "black", width = 0.2, size = 1) +
  scale_fill_manual(values = c("#FFC0CB", "#1E90FF"), name = "Cell type") +
  xlab("Condition") +
  ylab("Value") +
  ggtitle("Violin plot by condition and cell type")



library(ggplot2)

# Example data with random values
set.seed(123)
data <- data.frame(
  value = rnorm(1000),
  condition = rep(c("sick", "healthy"), each = 500),
  cell_type = rep(c("type1", "type2"), times = 500)
)

# Create violin plot
ggplot(data, aes(x = condition, y = value, fill = cell_type)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4")) +
  labs(x = "Condition", y = "Value", fill = "Cell Type")


library(ggplot2)

# Generate example data
set.seed(123)
sick_data <- data.frame(value = rnorm(100, mean = 10, sd = 2), group = "sick")
healthy_data <- data.frame(value = rnorm(100, mean = 8, sd = 1.5), group = "healthy")
data <- rbind(sick_data, healthy_data)

# Create violin plot
ggplot(data, aes(x = group, y = value)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_jitter(width = 0.1) +
  labs(x = "Group", y = "Value") +
  ggtitle("Distribution of Value by Group") +
  theme_bw()

