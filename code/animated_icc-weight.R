library(ggplot2)
library(tidyverse)
library(magick)
library(psych)
library(MASS)
library(grid)
library(patchwork)

create_gif <- function(path_files, delay, gif_out) {
  # get files and sort oldest to newest
  files_list <- list.files(path_files, full.names = TRUE)
  file_info <- file.info(files_list)
  file_list_old2new <- files_list[order(file_info$mtime)]
  # read in images from list, then convert to gif with fps-delay
  images <- lapply(file_list_old2new, image_read) 
  gif_img <- image_animate(image_join(images), delay = delay)
  # write gif out
  image_write(image = gif_img, path = gif_out)
}

pal<- c("#1D91C0","#67001F","#CB181D")

# inp dir
code_dir = "/Users/michaeldemidenko/Desktop/Academia/Stanford/6_Presentations/r_code"

### Observe= True + Error
set.seed(111)
# True weight
true_wght <- rep(171,100)
n_weights <- 100
obs_weight <- true_wght + rnorm(n_weights, mean = 0, sd = 1.5)

# Create a directory to store intermediate histogram images
dir.create(paste(code_dir,"/histograms"), showWarnings = FALSE)

# Generate histograms for each sample size and save them as images
for (i in 1:n_weights) {
  # Subset weights
  sample_weights <- obs_weight[1:i]
  
  # Create histogram
  p <- ggplot(data.frame(weights = sample_weights), aes(x = weights)) +
    geom_histogram(binwidth = 1, fill = "white", color = "black") +
    ggtitle(paste("Sample Size:", i)) +
    geom_vline(xintercept = true_wght, linetype = "dashed",alpha = .5) +
    ylim(0,50) +
    xlim(min(obs_weight),max(obs_weight)) +
    labs(x = "Weight", y = "Frequency") +
    theme_minimal()
  
  # Save the plot as an image
  ggsave(paste0("~/Downloads/histograms/histogram_", sprintf("%03d", i), ".jpeg"), p, width = 6, height = 4)
}

create_gif(path_files = paste(code_dir,"/histograms"), delay = 15, 
           gif_out = paste(code_dir,"/histograms/histograms.gif"))


### ICC variations
icc_out <- paste0(code_dir,"/icc_var/")
dir.create(icc_out, showWarnings = FALSE)

# set N & correlation size
n = 150
cor_r = .7

# set range and interval in means & variances
range_means = seq(0,4.95, by = .05)
range_vars = seq(1,10.9, by = .1)


for (mean in range_means) {

  # assign mean, var and estimate cov for specified r. Will vary mean 2
  mean1 <- 0
  mean2 <- mean
  mean3 <- mean - mean*.1
  var1 <- 3
  var2 <- 3
  var3 <- 3
  cov12 <- cor_r * sqrt(var1 * var2)
  cov13 <- cor_r * sqrt(var1 * var3)
  cov23 <- cor_r * sqrt(var2 * var3)
  
  cov_matrix <- matrix(
    c(var1, cov12, cov13,
      cov12, var2, cov23,
      cov13, cov23, var3),
    nrow = 3, ncol = 3,
    byrow = TRUE
  )
  
  # simulate data
  set.seed(123)
  dat <- data.frame(mvrnorm(n = n, 
                            mu = c(mean1, mean2,mean3),
                            Sigma = cov_matrix,
                            empirical = TRUE))
  colnames(dat) <- c("T1","T2","T3")
  
  diff_means = mean1 - mean2
  #calc correlation
  cor_r = round(cor(dat$T1, dat$T2),2)
  
  # calculate ICCs
  icc_1 <- ICC(dat, lmer = FALSE)$results$ICC[1]
  icc_2 <- ICC(dat, lmer = FALSE)$results$ICC[2]
  icc_3 <- ICC(dat, lmer = FALSE)$results$ICC[3]
  # create ICC df
  icc_df <- data.frame(
    icc = c("ICC1", "ICC2", "ICC3"),
    value = c(icc_1, icc_2, icc_3)
  )
  
  # plot population corr, means/vars boxplot and ICCs
  corr_plot <- dat %>% 
    ggplot(aes(x = T1, y = T2)) +
    geom_point(position = position_jitter()) +
    geom_smooth(method = "lm", formula = "y ~ x", colour = "red") + 
    labs(title = "Correlation between T1 & T2") +
    annotate(geom = "text",
             x = Inf, y = Inf, label = as.character(cor_r), 
             hjust = 1, vjust = 1)+
    theme_minimal()
  
  mean_plt <- dat %>% 
    gather(key = "variables",value = "scores", T1:T3) %>% 
    ggplot(aes(x = variables, y = scores, fill = variables)) +
    geom_boxplot() +
    labs(title = ("Scores Across Time"), 
         subtitle = paste0("Mu1: ",mean1, "; Mu2: ", mean2, "; Mu3: ", mean3)) +
    xlab("Variables") +
    scale_x_discrete()+
    scale_fill_manual(values = pal) +
    theme_minimal()
  
  icc_plt <- icc_df %>% 
    ggplot(aes(x = icc, y = value, fill = icc)) +
    geom_bar(stat = "identity") +
    labs(title = "ICC(1), ICC(2,1) and ICC(3,1) across Times")+
    scale_fill_manual(values = pal) +
    ylim(c(-1,1))+
    ylab("") +
    xlab("")+
    theme_minimal()
  
  # combined plots using patchwork
  plt = mean_plt + corr_plot + icc_plt
  plt = plt + plot_annotation(
    title = "Change in mean impact on ICCs?",
    caption = paste0("ICC extracted from Psych package ICC().\n","N:", n, " Empirical Corr: ", cor_r)
  )
  
  # create file name and save out png as specified size
  file_name = paste0(icc_out,"means/meansdiff",mean,".png")
  ggsave(filename = file_name, plot = plt, width = 15, height = 5, dpi = 200)
}


for (vars in range_vars) {
  # assign mean, var and estimate cov for specified r
  mean1 <- 0
  mean2 <- 0
  mean3 <- 0
  
  var1 <- 3
  var2 <- 3 + vars
  var3 <- 3
  
  cov12 <- cor_r * sqrt(var1 * var2)
  cov13 <- cor_r * sqrt(var1 * var3)
  cov23 <- cor_r * sqrt(var2 * var3)
  
  cov_matrix <- matrix(
    c(var1, cov12, cov13,
      cov12, var2, cov23,
      cov13, cov23, var3),
    nrow = 3, ncol = 3,
    byrow = TRUE
  )
  
  # simulate data
  set.seed(123)
  dat <- data.frame(mvrnorm(n = n, 
                            mu = c(mean1, mean2, mean3),
                            Sigma = cov_matrix,
                            empirical = TRUE))
  colnames(dat) <- c("T1","T2","T3")
  
  diff_means = mean1 - mean2
  #calc correlation
  cor_r = round(cor(dat$T1, dat$T2),2)
  
  # calculate ICCs
  icc_1 <- ICC(dat, lmer = FALSE)$results$ICC[1]
  icc_2 <- ICC(dat, lmer = FALSE)$results$ICC[2]
  icc_3 <- ICC(dat, lmer = FALSE)$results$ICC[3]
  # create ICC df
  icc_df <- data.frame(
    icc = c("ICC1", "ICC2", "ICC3"),
    value = c(icc_1, icc_2, icc_3)
  )
  
  # plot population corr, means/vars boxplot and ICCs
  corr_plot <- dat %>% 
    ggplot(aes(x = T1, y = T2)) +
    geom_point(position = position_jitter()) +
    geom_smooth(method = "lm", formula = "y ~ x", colour = "red") + 
    labs(title = "Correlation between T1 & T2") +
    annotate(geom = "text",
             x = Inf, y = Inf, label = as.character(cor_r), 
             hjust = 1, vjust = 1)+
    theme_minimal()
  
  mean_plt <- dat %>% 
    gather(key = "variables",value = "scores", T1:T3) %>% 
    ggplot(aes(x = variables, y = scores, fill = variables)) +
    geom_boxplot() +
    labs(title = ("Scores across Time"), 
         subtitle = paste0("Var1: ",var1, "; Var2: ", var2, "; Var3: ", var3)) +
    xlab("Variables") +
    scale_x_discrete()+
    scale_fill_manual(values = pal) +
    theme_minimal()
  
  icc_plt <- icc_df %>% 
    ggplot(aes(x = icc, y = value, fill = icc)) +
    geom_bar(stat = "identity") +
    labs(title = "ICC(1), ICC(2,1) and ICC(3,1) across Measurement Occassions")+
    scale_fill_manual(values = pal) +
    ylim(c(-1,1))+
    ylab("") +
    xlab("")+
    theme_minimal()
  
  # combined plots using patchwork
  plt = mean_plt + corr_plot + icc_plt
  plt = plt + plot_annotation(
    title = "Change in variances impact on ICCs?",
    caption = paste0("ICC extracted from Psych package ICC().\n","N:", n, " Empirical Corr: ", cor_r)
  )
  
  # create file name and save out png as specified size
  file_name = paste0(icc_out,"variances/varsdiff",vars,".png")
  ggsave(filename = file_name, plot = plt, width = 15, height = 5, dpi = 200)
  
}

create_gif(path_files = paste0(icc_out,"/means/"), delay = 15, 
           gif_out = paste0(icc_out,"mean_impact-icc.gif"))

create_gif(path_files = paste0(icc_out,"/variances/"), delay = 15, 
           gif_out = paste0(icc_out,"variances_impact-icc.gif"))

