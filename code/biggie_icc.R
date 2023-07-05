library(psych)
library(MASS)
library(png)
library(grid)
library(patchwork)
library(tidyverse)
library(viridis)

output <- "./out/icc/"
big_img <- "../images/"

# set N & correlation size
n = 100
cor_r = .8

# set range and interval in means & variances
range_means = seq(0,4.95, by = .05)
range_vars = seq(1,10.9, by = .1)

count_means <- 1
for (mean in range_means) {
  
  if (count_means < 15) {
    count_means <- 1
  } else {
    count_means <- count_means + 1
  }
  
  # load image & rasterize it, reducing size and adjusting to bottom/left corner
  big = readPNG(paste0("../images/dream",count_means,".png"))
  big_grob <- rasterGrob(big, interpolate = TRUE,width = .35, height = .35, just = c(1.3,1.3))
  
  # assign mean, var and estimate cov for specified r. Will vary mean 2
  mean1 <- 0
  mean2 <- mean
  var1 <- 3
  var2 <- 3
  cov <- cor_r *sqrt(var1*var2)
  
  cov_matrix <- matrix(
    c(var1, cov,
      cov,var2), 
    nrow = 2, ncol = 2, 
    byrow = TRUE
  )
  
  # simulate data
  set.seed(111)
  dat <- data.frame(mvrnorm(n = n, 
                            mu = c(mean1, mean2),
                            Sigma = cov_matrix,
                            empirical = TRUE))
  colnames(dat) <- c("T1","T2")
  
  diff_means = mean1 - mean2
  #calc correlation
  cor_r = cor(dat$T1, dat$T2)
  
  # calculate ICCs
  icc_1 <- ICC(dat)$results$ICC[1]
  icc_2 <- ICC(dat)$results$ICC[2]
  icc_3 <- ICC(dat)$results$ICC[3]
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
             x = Inf, y = Inf, label = as.character(round(cor_r,2)), 
             hjust = 1, vjust = 1)+
    theme_minimal()
  
  mean_plt <- dat %>% 
    gather(key = "variables",value = "scores", T1:T2) %>% 
    ggplot(aes(x = variables, y = scores, fill = variables)) +
    geom_boxplot() +
    labs(title = ("Boxplots of T1 & T2"), 
         subtitle = paste0("Mu1: ",mean1, "; Mu2: ", mean2)) +
    xlab("Variables") +
    scale_x_discrete()+
    scale_fill_brewer(palette = "Set1")+
    theme_minimal()
  
  icc_plt <- icc_df %>% 
    ggplot(aes(x = icc, y = value, fill = icc)) +
    geom_bar(stat = "identity") +
    labs(title = "ICC(1), ICC(2,1) and ICC(3,1) across Time1 & Time2")+
    scale_fill_brewer(palette = "Set1")+
    ylim(c(-1,1))+
    ylab("") +
    xlab("")+
    theme_minimal()
  icc_plt <- icc_plt + annotation_custom(big_grob, xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf)
  
  # combined plots using patchwork
  plt = mean_plt + corr_plot + icc_plt
  plt = plt + plot_annotation(
    title = "T2 change in mean impact on ICCs?",
    caption = "giggles, no?"
  )
  
  # create file name and save out png as specified size
  file_name = paste0(output,"meansdiff",mean,".png")
  ggsave(filename = file_name, plot = plt, width = 15, height = 5, dpi = 200)
}


count_var <- 1
for (vars in range_vars) {
  
  if (count_var < 15) {
    count_var <- 1
  } else {
    count_var <- count_var + 1
  }
  
  # load image
  big = readPNG(paste0("../images/dream",count_var,".png"))
  big_grob <- rasterGrob(big, interpolate = TRUE,width = .35, height = .35, just = c(1.3,1.3))
  
  # assign mean, var and estimate cov for specified r
  mean1 <- 2
  mean2 <- 2.5
  var1 <- 2
  var2 <- 2 + vars
  cov <- cor_r *sqrt(var1*var2)
  
  cov_matrix <- matrix(
    c(var1, cov,
      cov,var2), 
    nrow = 2, ncol = 2, 
    byrow = TRUE
  )
  
  # simulate data
  set.seed(111)
  dat <- data.frame(mvrnorm(n = n, 
                            mu = c(mean1, mean2),
                            Sigma = cov_matrix,
                            empirical = TRUE))
  colnames(dat) <- c("T1","T2")
  
  diff_means = mean1 - mean2
  #calc correlation
  cor_r = cor(dat$T1, dat$T2)
  
  # calculate ICCs
  icc_1 <- ICC(dat)$results$ICC[1]
  icc_2 <- ICC(dat)$results$ICC[2]
  icc_3 <- ICC(dat)$results$ICC[3]
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
             x = Inf, y = Inf, label = as.character(round(cor_r,2)), 
             hjust = 1, vjust = 1)+
    theme_minimal()
  
  mean_plt <- dat %>% 
    gather(key = "variables",value = "scores", T1:T2) %>% 
    ggplot(aes(x = variables, y = scores, fill = variables)) +
    geom_boxplot() +
    labs(title = ("Boxplots of T1 & T2"), 
         subtitle = paste0("Var1: ",var1, "; Var2: ", var2)) +
    xlab("Variables") +
    scale_x_discrete()+
    scale_fill_brewer(palette = "Set1")+
    theme_minimal()
  
  icc_plt <- icc_df %>% 
    ggplot(aes(x = icc, y = value, fill = icc)) +
    geom_bar(stat = "identity") +
    labs(title = "ICC(1), ICC(2,1) and ICC(3,1) across Time1 & Time2")+
    scale_fill_brewer(palette = "Set1")+
    ylim(c(-1,1))+
    ylab("") +
    xlab("")+
    theme_minimal()
  icc_plt <- icc_plt + annotation_custom(big_grob, xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf)
  
  # combined plots using patchwork
  plt = mean_plt + corr_plot + icc_plt
  plt = plt + plot_annotation(
    title = "T2 change in variance impact on ICCs?",
    caption = "giggles, no?"
  )
  
  # create file name and save out png as specified size
  file_name = paste0(output,"varsdiff",vars,".png")
  ggsave(filename = file_name, plot = plt, width = 15, height = 5, dpi = 200)
}
