# CSRT Threshold validation: Validation iteration
# October 9, 2025

# Data contains nearest neighbor climate distance estimates (Edist) for 500m elevation transects of a set of randomized validation points
# The original dataset is 1GB so the sample file here has been truncated, but the original results output has been provided at the end of the script.

library(tidyverse)
library(skimr)
library(yaImpute)
library(gower)
library(lmtest) #dwtest
library(fBasics) #JarqueBeraTest
library(broom) #lm summary


# Load output from nearest neighbor Edist elevation transect
# (file size truncated for github)
output <- read.csv("data/git_makepts_edist.csv",header=TRUE)

# Count original locations used for elevation gradient
loc <- output %>% 
  distinct(loc)
nrow(loc)

# Bootstrap
results_list <- vector("list", 100)

for (i in 1:100) {
  
  # Subsample 90%
  set.seed(i)
  output_sub <- output[sample(nrow(output), size = 0.9 * nrow(output)), ]
  
  # Process data
  q <- output_sub %>%
    mutate(diff = abs(Nelev - elev)) %>%
    filter(Edist > 0) %>%
    select(Edist, diff) %>%
    group_by(diff) %>%
    summarise(meanE = mean(Edist),
              sdE = sd(Edist), .groups = "drop")
  
  # Run linear model
  qreg <- lm(meanE ~ diff, data = q)
  
  # Thresholds (with NA safety if diff values not found)
  thresh.300 <- q %>% filter(diff == 300) %>% pull(meanE)
  thresh.200 <- q %>% filter(diff == 200) %>% pull(meanE)
  mid <- if (length(thresh.200) == 1 && length(thresh.300) == 1) {
    thresh.200 + 0.5 * (thresh.300 - thresh.200)
  } else {
    NA
  }
  
  # Assumption checks with robust error handling
  jb_test <- tryCatch(jarqueberaTest(qreg$residuals), error = function(e) NULL)
  dw_test <- tryCatch(dwtest(qreg), error = function(e) NULL)
  
  # Extract p-values only if tests succeeded
  jb_pval <- if (!is.null(jb_test) && inherits(jb_test, "fHTEST")) jb_test@test$p.value else NA
  dw_pval <- if (!is.null(dw_test) && inherits(dw_test, "htest")) dw_test$p.value else NA
  
  # Model stats
  glance_stats <- glance(qreg)
  adj_r2   <- glance_stats$adj.r.squared
  model_p  <- glance_stats$p.value
  df_resid <- glance_stats$df.residual
  
  # Combine all results
  results_list[[i]] <- data.frame(
    iteration = i,
    thresh.200 = ifelse(length(thresh.200) == 0, NA, thresh.200),
    mid = mid,
    thresh.300 = ifelse(length(thresh.300) == 0, NA, thresh.300),
    adj_r2 = adj_r2,
    model_p = model_p,
    df_resid = df_resid,
    jb_pval = jb_pval,
    dw_pval = dw_pval
  )
}

# Combine into final table
results_df <- do.call(rbind, results_list)

# Save
# write.csv(results_df, "output/thresh_iteration_stats.csv", row.names = FALSE)
# Original results file provided output/thresh_iteration_stats_09_100.csv

# Optional: View result
print(results_df)
mean(results_df$thresh.200)
mean(results_df$thresh.300)
sd(results_df$thresh.200)
sd(results_df$thresh.300)
