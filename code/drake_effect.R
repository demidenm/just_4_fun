library(MASS)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(png)
library(grid)


# set output to write figures to
output = "~/Downloads/animated/rcorr/"

drake_yay = readPNG("~/Downloads/animated/drake_sig.png")
drake_yay_grob <- rasterGrob(drake_yay, interpolate = TRUE,width = .5, height = .5, just = c(1,1))
drake_boo = readPNG("~/Downloads/animated/drake_nonsig.png")
drake_boo_grob <- rasterGrob(drake_boo, interpolate = TRUE,width = .5, height = .5, just = c(1,1))

# set means and true are for empirical correlation
means = c(5,5)
true_r = .30
r = matrix(c(1,true_r,
             true_r,1), 
  nrow = 2)
# create data.frame for above
df_r1 = as.data.frame(
  mvrnorm(n = 10000,
         mu = means,
         Sigma = r,
         empirical = TRUE)
)

# set parameters for resampled correlatoins, the n number to iteration, subsample list and dataframe to save to
n_iterations <- 200
subsamp_list <- list()
sub_n = 30
results <- data.frame(iteration = integer(),
                      estimate = numeric(),
                      pvalue = numeric())


# loop over iterations to create the separate datasets
for (i in 1:n_iterations) {
  # Sample w/o replacement n = 50
  subsamp <- df_r1 %>% sample_n(size = sub_n, replace = FALSE)
  
  # Store the subsample in the list
  subsamp_list[[i]] <- subsamp
  
  # Calculate the correlation between V1 and V2
  corr <- cor.test(subsamp$V1, subsamp$V2)
  
  # Extract the estimate and p-value
  r_val <- corr$estimate
  r_pval <- round(corr$p.value,3)
  
  # Append the results to the dataframe
  results <- results %>% add_row(iteration = i, estimate = r_val, pvalue = r_pval)
}

# select slice of first correlation that is at least .05. Most convienient way to get at which r <= .05
sig_r = results %>% select(estimate, pvalue) %>% 
  filter(pvalue <= .05) %>% 
  arrange(desc(pvalue)) %>% 
  slice(1)

# Population scatter plot to reuse
population <- df_r1 %>% 
  ggplot(aes(x = V1, y = V2)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y ~ x", se = FALSE, color = "red") +
  ggtitle("Population Correlation") +
  labs(x = "V1", y = "V2") +
  theme_minimal()

# loop over iterations to subsample our list to create custom distirbution of r and subsample plot
# save images to folder
for (df in 1:n_iterations) {
  subsamp_plot = subsamp_list[[df]] %>% 
    ggplot(aes(x = V1, y = V2)) +
    geom_point() +
    geom_smooth(method = "lm", formula = "y ~ x", se = FALSE, color = "red") +
    ggtitle(paste("Correlation for Subsample:", df)) +
    annotate(geom = "text",
             x = Inf, y = Inf, label = as.character(round(results[df,2],2)), 
             hjust = 1, vjust = 1)+
    labs(x = "V1", y = "V2") +
    theme_minimal()
  # if p < .05 or > .05
  if (results[df,3] < .05) {
    subsamp_plot = subsamp_plot + annotation_custom(drake_yay_grob, xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf)
  } else  {
    subsamp_plot = subsamp_plot + annotation_custom(drake_boo_grob, xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf)
  }
  
  distribution_r <- results %>% 
    ggplot(aes(x = estimate)) +
    geom_histogram(bins = 35, color = "black", fill = "white") +
    geom_vline(xintercept = true_r, color = "blue", linetype = "dashed", size = 1.3) +
    geom_vline(xintercept = sig_r[1,1], color = "black", linetype = "solid", size = 1.3) +
    geom_vline(xintercept = as.numeric(round(results[df,2],2)), linetype = "dashed", color = "red",size = 1)+
    labs(title = paste("Distribution of r values for samples N =", sub_n), 
         subtitle = "BLUE: population r; BLACK: p < .05 boundary; RED: Subsample r") +
    theme_minimal()
  
  plt = population + distribution_r + subsamp_plot
  file_name = paste0(output,df,".png")
  ggsave(filename = file_name, plot = plt, width = 15, height = 5, dpi = 200)
  
}
