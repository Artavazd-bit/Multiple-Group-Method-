library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(stringr)

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
