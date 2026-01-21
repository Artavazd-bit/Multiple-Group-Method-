library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(stringr)
library(skimr)

res <- read.csv2("simresults/simresults_calc2.csv")

test <- res %>% 
          count(test, warning) %>% 
          mutate(rel_freq = n/sum(n))

res2 <- res

res2$issue <- if(!is.na(res2$error) | !is.na(res2$warning) | !is.na(res2$constrained_warning) | !is.na(res2$unconstrained_warning) | !is.na(res2$constrained_error) | !is.na(res2$unconstrained_error)) 
  
  
  
res2 <- res2 %>% 
        mutate(issue = case_when(
          !is.na(error) ~ 1,
          !is.na(warning) ~ 1,
          !is.na(constrained_warning) ~ 1,
          !is.na(unconstrained_warning) ~ 1,
          !is.na(constrained_error) ~ 1,
          !is.na(unconstrained_error) ~ 1,
          TRUE ~ 0
        ))


test <- res2  %>% 
          count(test, sample_size, correlation, type, issue) %>% 
          mutate(rel_freq = n/sum(n)*100) %>% 
          arrange(rel_freq)

issues <- test[test$issue == 1,]


write.csv2(issues, "issues.csv")


numeric_summary_long <- res %>%
  group_by(type, correlation, sample_size, test) %>%
  summarise(
    across(
      c(stat, constrained, unconstrained, FL_1, FL_2),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        var = ~var(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE),
        p05 = ~quantile(.x, 0.05, na.rm = TRUE),
        p25 = ~quantile(.x, 0.25, na.rm = TRUE),
        p50 = ~quantile(.x, 0.50, na.rm = TRUE),
        p75 = ~quantile(.x, 0.75, na.rm = TRUE),
        p95 = ~quantile(.x, 0.95, na.rm = TRUE)
      )
    ),
    n = n(),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -c(type, correlation, sample_size, test, n),
    names_to = c("variable", "statistic"),
    names_sep = "_(?=[^_]+$)",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = statistic,
    values_from = value
  )

numeric_summary_long_filt <- numeric_summary_long %>% 
  filter(!is.na(mean))

warning_table <- res %>%
  count(type, correlation, sample_size, test, warning) %>%
  arrange(type, correlation, sample_size, test, desc(n))


error_table <- res %>%
  count(type, correlation, sample_size, test, error) %>%
  arrange(type, correlation, sample_size, test, desc(n))

error_table_filt <- error_table %>% 
  filter(!is.na(error))

warning_table_filt <- warning_table %>% 
  filter(!is.na(warning))

p_error_facet <- error_table_filt %>%
  mutate(
    error_short = str_trunc(error, 50),
    error_short = fct_reorder(error_short, n, .fun = sum)
  ) %>%
  ggplot(aes(x = n, y = error_short, fill = test)) +
  geom_col(position = "dodge") +
  facet_grid(type + sample_size ~ correlation,
             labeller = labeller(
               sample_size = ~paste("n =", .),
               correlation = ~paste("r =", .)
             ),
             scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Error Frequency by Condition and Type",
    x = "Count",
    y = "Error Message",
    fill = "Test"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    strip.background = element_rect(fill = "grey90", color = NA),
    legend.position = "bottom"
  )

# --- Warning Plot with type as facet row ---
p_warning_facet <- warning_table_filt %>%
  mutate(
    warning_short = str_trunc(warning, 50),
    warning_short = fct_reorder(warning_short, n, .fun = sum)
  ) %>%
  ggplot(aes(x = n, y = warning_short, fill = test)) +
  geom_col(position = "dodge") +
  facet_grid(type + sample_size ~ correlation,
             labeller = labeller(
               sample_size = ~paste("n =", .),
               correlation = ~paste("r =", .)
             ),
             scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Warning Frequency by Condition and Type",
    x = "Count",
    y = "Warning Message",
    fill = "Test"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    strip.background = element_rect(fill = "grey90", color = NA),
    legend.position = "bottom"
  )



res %>% 
  group_by(test, type, correlation, sample_size) %>% 
  summarise(across(where(is.numeric)), list(mean = mean, sd = sd, na.rm = TRUE)
  )

skimr::skim(res)

skim <- res %>%
  group_by(test, type, correlation, sample_size) %>%
  skim(stat, FL_1, FL_2, constrained, unconstrained)


skim2 <- skim[!is.na(skim$numeric.p75), ]


res3 <- res[res$test == "FL",]

res4 <- res3[!is.na(res3$warning),]

nthroot = function(x,n) {(abs(x)^(1/n))*sign(x)}



res[!is.na(res$warning) & res$test== "MGA",]

data <-  lavaan::simulateData(model = simModels$model[1],
                              sample.nobs = 50, # Number of observations.
                              skewness = NULL,
                              kurtosis = NULL,
                              seed = 220649496, # Set random seed.
                              empirical = FALSE, # Logical. If TRUE, the implied moments (Mu and Sigma) specify the empirical not population mean and covariance matrix.
                              return.type = "data.frame"
)




