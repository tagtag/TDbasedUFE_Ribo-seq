#MOFA2
dim(Z)<- c(dim(Z)[1],18,3) #not scaled

mRNA <- Z[,,1];rownames(mRNA) <- x[,1]
trans <- Z[,,2]; rownames(trans) <- x[,1]
protein <- Z[,,3] ; rownames(protein) <- x[,1]
protein <- protein[apply(protein, 1, sd) > 0, ]
require(MOFA2)
# Compile the three omics datasets into a list
data_list <- list(
  mRNA = mRNA,       # Transcriptome
  trans = trans,   # translatome (Ribo-seq)
  protein = protein #  proteome
)

#Creating a MOFA object (the same applies to three-tier systems)
MOFAobject <- create_mofa(data_list)
# Retrieving Model Options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10

# Training Option: Relax ARD (Feature Selection)
train_opts <- get_default_training_options(MOFAobject)
# Increase the number of training iterations so that even very small changes are not overlooked
train_opts$maxiter <- 2000
train_opts$convergence_mode <- "slow"
train_opts$drop_factor_threshold <- 0

# [Important] Execute this function separately before running `run_mofa`
MOFAobject <- prepare_mofa(MOFAobject, 
                           model_options = model_opts, 
                           training_options = train_opts)

# Run Training
MOFAobject <- run_mofa(MOFAobject, use_basilisk = FALSE)


plot_variance_explained(MOFAobject)
 get_dimensions(MOFAobject)$K


library(MOFA2)
library(ggplot2)
library(dplyr)
library(patchwork) # For arranging plots

# 1. Extraction and Aggregation of Weight Data (Factors 1–3, All Views)
# Here, we use the average to determine the "overall trend."
summary_weights <- get_weights(MOFAobject, factors = 1:4, as.data.frame = TRUE) %>%
  group_by(factor, view) %>%
  summarize(value = mean(value), .groups = 'drop') %>%
  # Convert the view name to a number (k) (e.g., k=1, 2, 3, 4)
  mutate(k = as.numeric(as.factor(view)))

# 2. A function to create individual plots for each factor
create_stem_plot <- function(f_name, y_label) {
  df_sub <- filter(summary_weights, factor == f_name)
  
  ggplot(df_sub, aes(x = k, y = value)) +
    geom_segment(aes(xend = k, yend = 0), size = 1) + # 0から伸びる垂直線
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") + # 基準線
    scale_x_continuous(breaks = c(1, 2, 3), limits = c(0.8, 3.2)) +
    labs(y = y_label, x = "k") +
    theme_bw() +
    theme(panel.grid = element_blank())
}

# 3. Create three plots and arrange them vertically
p1 <- create_stem_plot("Factor1", "U1k")
p2 <- create_stem_plot("Factor2", "U2k")
p3 <- create_stem_plot("Factor3", "U3k")
p4 <- create_stem_plot("Factor1", "U4k")

# Show All
p1 / p2 / p3 / p4
