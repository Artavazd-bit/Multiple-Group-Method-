library(dplyr)

res <- read.csv2("simresults_calc2.csv")

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
