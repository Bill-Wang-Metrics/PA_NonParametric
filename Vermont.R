########## EBT rollout timing differences
# VIM
a <- c(2, 4, 4, 9, 11)
b <- c(1, 2, 6, 7, 7, 9, 11)
wilcox.test(a, b, alternative = "two.sided")  


# ARM
a <- c(2, 6, 9, 11)
b <- c(1, 2, 4, 4, 6, 7, 7, 9, 11)
wilcox.test(a, b, alternative = "two.sided") 

####################### VIM

# WIC Retention Data (%)
sites <- c("Burlington", "Bennington*", "White River", "Rutland*", "Springfield",
           "Newport*", "St. Johnsbury*", "Barre", "Brattleboro", "St. Albans*",
           "Middlebury*", "Morrisville*")
year_2015 <- c(68.1, 74.8, 70.8, 69.9, 66.9, 81.5, 79.9, 66.8, 72.6, 73.2, 84.0, 83.8)
year_2016 <- c(62.6, 73.7, 63.3, 65.9, 62.5, 73.7, 77.3, 61.4, 70.3, 70.4, 79.7, 84.3)
year_2017 <- c(62.3, 63.3, 66.3, 59.2, 64.7, 59.4, 70.0, 61.9, 67.6, 59.0, 72.4, 77.9)

# Create treatment indicator (0 = control/*, 1 = treatment)
D <- c( 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0)


library(combinat) # For combn() function

par(mar = c(5, 6, 4, 2) + 0.1)

# Data from your table (2015 vs 2017 changes)
changes <- year_2017 - year_2015
treatment_indices <- which(D == 1)

# Function to calculate DiD for given treatment assignment
calc_did <- function(treated_indices) {
  mean(changes[treated_indices]) - mean(changes[-treated_indices])
}


# Generate all possible treatment assignments
all_combinations <- combn(1:12, 5, simplify = FALSE)

# Calculate DiD for EVERY possible assignment
all_dids <- sapply(all_combinations, calc_did)

# Actual DiD
actual_did <- calc_did(treatment_indices)

# Exact p-value (one-tailed)
exact_p <- mean(all_dids >= actual_did)

# Results
cat("Exact DiD:", actual_did, "\nExact p-value:", exact_p)

hist(all_dids, breaks=seq(-10,10,0.5), col="skyblue",
     main="",ylim = c(0,100),xlim = c(-10,10),
     xlab="Possible DiD Estimates",
     cex.main = 2, cex.lab=2,cex.axis = 2)
abline(v=actual_did, col="red", lwd=2.5)


# Pre-period changes (2015 to 2016)
pre_changes <- year_2016 - year_2015

# Function to calculate pre-trend DiD
calc_pretrend_did <- function(treated_indices) {
  mean(pre_changes[treated_indices]) - mean(pre_changes[-treated_indices])
}

# Actual pre-trend DiD (true treatment assignment)
actual_pretrend_did <- calc_pretrend_did(which(D == 1))

# Generate all possible treatment assignments (same as before)
all_combinations <- combn(1:12, 5, simplify = FALSE)

# Calculate pre-trend DiD for EVERY possible assignment
all_pretrend_dids <- sapply(all_combinations, calc_pretrend_did)

# Two-tailed p-value for pre-trend test
exact_pretrend_p <- mean(abs(all_pretrend_dids) >= abs(actual_pretrend_did))

# Lower-tailed p-value for pre-trend test
exact_pretrend_p <- mean((all_pretrend_dids) <= (actual_pretrend_did))
exact_pretrend_p

hist(all_pretrend_dids, breaks=seq(-10,10,0.5), col="lightgreen",
     main="",ylim = c(0,100),xlim = c(-10,10),
     xlab="Possible Pre-Trend Estimates",
     cex.main = 2, cex.lab=2,cex.axis = 2)
abline(v=actual_pretrend_did, col="red", lwd=2.5)
abline(v=-actual_pretrend_did, col="red", lty=2.5)  # Two-tailed critical region



########### DiD
library(tidyr)
library(dplyr)
library(ggplot2)

# Create dataframe
wic_data <- data.frame(
  site = sites,
  treated = D,
  y2015 = year_2015,
  y2016 = year_2016,
  y2017 = year_2017
)

# Convert to long format
long_data <- wic_data %>%
  pivot_longer(cols = starts_with("y"), 
               names_to = "year",
               values_to = "retention") %>%
  mutate(year = as.numeric(gsub("y", "", year)),
         post = as.numeric(year >= 2017))  # Treatment effect kicks in year 2017

library(fixest) # For efficient fixed effects models

fe_model <- feols(retention ~ treated*post | site + year, data = long_data)
summary(fe_model)


# Create a "time to treatment" variable
long_data <- long_data %>%
  mutate(post1 =  (year == 2017), pre1 = (year == 2015)) 

# Event-study regression
event_study <- feols(retention ~ pre1*treated + post1*treated | site + year, data = long_data) 

summary(event_study)


coefs <- summary(event_study)$coefficients
plot_data <- data.frame(
  term = c("Pre(2015)", "Post(2017)"),
  estimate = coefs,
  ci_lower = coefs - 1.96*summary(event_study)$se, 
  ci_upper = coefs + 1.96*summary(event_study)$se,
  period = c(-1,1)
)

# Add implementation year (2016) as reference point (effect = 0)
plot_data <- rbind(
  data.frame(term = "Implementation (2016)", estimate = 0, 
             ci_lower = NA, ci_upper = NA, period = 0),
  plot_data)

# Create plot with connected line
ggplot(plot_data, aes(x = period, y = estimate)) +
  # Reference elements
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  
  # Connecting line (only between points, excluding CI ranges)
  geom_line(aes(group = 1), color = "#0072B2", alpha = 0.5, linewidth = 0.8) +
  
  # Points and error bars (filter out NA for reference point)
  geom_point(data = filter(plot_data, !is.na(ci_lower)), 
             size = 3, color = "#0072B2") +
  geom_errorbar(data = filter(plot_data, !is.na(ci_lower)),
                aes(ymin = ci_lower, ymax = ci_upper), 
                width = 0.2, color = "#0072B2", linewidth = 1) +
  
  # Reference point (2016)
  geom_point(data = filter(plot_data, is.na(ci_lower)), 
             size = 3, shape = 21, fill = "white", color = "#0072B2") +
  
  # Axis and labels
  scale_x_continuous(
    breaks = c(-1, 0, 1), 
    labels = c("Pre (2015)", "Implementation (2016)", "Post (2017)"),
    limits = c(-1.5, 1.5)
  ) +
  labs(
    # title = "Event-Study DiD: WIC Retention Effects",
    # subtitle = "Connected line shows temporal pattern",
    x = "",
    y = "Treatment Effect (%)",
    # caption = "Error bars show 95% confidence intervals\nReference point at 2016 normalized to 0"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )




####################### ARM

# WIC Retention Data (%)
sites <- c("Burlington", "Bennington*", "White River", "Rutland*", "Springfield",
           "Newport*", "St. Johnsbury*", "Barre", "Brattleboro", "St. Albans*",
           "Middlebury*", "Morrisville*")
year_2015 <- c(68.1, 74.8, 70.8, 69.9, 66.9, 81.5, 79.9, 66.8, 72.6, 73.2, 84.0, 83.8)
year_2017 <- c(62.3, 63.3, 66.3, 59.2, 64.7, 59.4, 70.0, 61.9, 67.6, 59.0, 72.4, 77.9)

# Create treatment indicator (0 = control/*, 1 = treatment)
D <- c(0,1,0,0,0,1,0,1,0,0,1,0)

year_2016 <- c(62.6, 73.7, 63.3, 65.9, 62.5, 73.7, 77.3, 61.4, 70.3, 70.4, 79.7, 84.3)


library(combinat) # For combn() function

# Data from your table (2015 vs 2017 changes)
changes <- year_2017 - year_2015
treatment_indices <- which(D == 1)

# Function to calculate DiD for given treatment assignment
calc_did <- function(treated_indices) {
  mean(changes[treated_indices]) - mean(changes[-treated_indices])
}

# Generate all possible treatment assignments
all_combinations <- combn(1:12, 4, simplify = FALSE)

# Calculate DiD for EVERY possible assignment
all_dids <- sapply(all_combinations, calc_did)

# Actual DiD
actual_did <- calc_did(treatment_indices)

# Exact p-value (one-tailed)
exact_p <- mean(all_dids >= actual_did)

# Results
cat("Exact DiD:", actual_did, "\nExact p-value:", exact_p)

hist(all_dids, breaks=seq(-10,10,0.5), col="skyblue",
     main="",ylim = c(0,100),xlim = c(-10,10),
     xlab="Possible DiD Estimates",
     cex.main = 2, cex.lab=2,cex.axis = 2)
abline(v=actual_did, col="red", lwd=2.5)


# Pre-period changes (2015 to 2016)
pre_changes <- year_2016 - year_2015

# Function to calculate pre-trend DiD
calc_pretrend_did <- function(treated_indices) {
  mean(pre_changes[treated_indices]) - mean(pre_changes[-treated_indices])
}

# Actual pre-trend DiD (true treatment assignment)
actual_pretrend_did <- calc_pretrend_did(which(D == 1))

# Generate all possible treatment assignments (same as before)
all_combinations <- combn(1:12, 4, simplify = FALSE)

# Calculate pre-trend DiD for EVERY possible assignment
all_pretrend_dids <- sapply(all_combinations, calc_pretrend_did)

# Two-tailed p-value for pre-trend test
exact_pretrend_p <- mean(abs(all_pretrend_dids) >= abs(actual_pretrend_did))

hist(all_pretrend_dids, breaks=seq(-10,10,0.5), col="lightgreen",
     main="",ylim = c(0,100),xlim = c(-10,10),
     xlab="Possible Pre-Trend DiD Estimates",
     cex.main = 2, cex.lab=2,cex.axis = 2)
abline(v=actual_pretrend_did, col="red", lwd=2.5)
abline(v=-actual_pretrend_did, col="red", lty=2.5)  # Two-tailed critical region



########### DiD
library(tidyr)
library(dplyr)

# Create dataframe
wic_data <- data.frame(
  site = sites,
  treated = D,
  y2015 = year_2015,
  y2016 = year_2016,
  y2017 = year_2017
)

# Convert to long format
long_data <- wic_data %>%
  pivot_longer(cols = starts_with("y"), 
               names_to = "year",
               values_to = "retention") %>%
  mutate(year = as.numeric(gsub("y", "", year)),
         post = as.numeric(year >= 2017))  # Treatment occurred after 2016

library(fixest) # For efficient fixed effects models

fe_model <- feols(retention ~ treated*post | site + year, data = long_data)
summary(fe_model)


# Create a "time to treatment" variable
long_data <- long_data %>%
  mutate(post1 =  (year == 2017), pre1 = (year == 2015)) 

# Event-study regression
event_study <- feols(retention ~ pre1*treated + post1*treated | site + year, data = long_data) 

summary(event_study)


coefs <- summary(event_study)$coefficients
plot_data <- data.frame(
  term = c("Pre(2015)", "Post(2017)"),
  estimate = coefs,
  ci_lower = coefs - 1.96*summary(event_study)$se, 
  ci_upper = coefs + 1.96*summary(event_study)$se,
  period = c(-1,1)
)

# Add implementation year (2016) as reference point (effect = 0)
plot_data <- rbind(
  data.frame(term = "Implementation (2016)", estimate = 0, 
             ci_lower = NA, ci_upper = NA, period = 0),
  plot_data)

# Create plot with connected line
ggplot(plot_data, aes(x = period, y = estimate)) +
  # Reference elements
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  
  # Connecting line (only between points, excluding CI ranges)
  geom_line(aes(group = 1), color = "#0072B2", alpha = 0.5, linewidth = 0.8) +
  
  # Points and error bars (filter out NA for reference point)
  geom_point(data = filter(plot_data, !is.na(ci_lower)), 
             size = 3, color = "#0072B2") +
  geom_errorbar(data = filter(plot_data, !is.na(ci_lower)),
                aes(ymin = ci_lower, ymax = ci_upper), 
                width = 0.2, color = "#0072B2", linewidth = 1) +
  
  # Reference point (2016)
  geom_point(data = filter(plot_data, is.na(ci_lower)), 
             size = 3, shape = 21, fill = "white", color = "#0072B2") +
  
  # Axis and labels
  scale_x_continuous(
    breaks = c(-1, 0, 1), 
    labels = c("Pre (2015)", "Implementation (2016)", "Post (2017)"),
    limits = c(-1.5, 1.5)
  ) +
  labs(
    # title = "Event-Study DiD: WIC Retention Effects",
    # subtitle = "Connected line shows temporal pattern",
    x = "",
    y = "Treatment Effect (%)",
    # caption = "Error bars show 95% confidence intervals\nReference point at 2016 normalized to 0"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )
