library(lavaan)
library(semTools)
library(doParallel)
library(foreach)
library(dplyr)
library(ggplot2)
library(skimr)

source("2026_01_20_setup_v1.R")

res <- read.csv2("simresults/simres_v1.csv")
res2 <- res

all_warnings <- unique(res$warning)
all_constrained_warning <- unique(res$warning_constrained)

all_error <- unique(res$error)
all_constrained_error <- unique(res$error_constrained)
# there are no errors
# i will focus on the warnings

# check whether there are cases that are NA but without warning
na_without_warning <- res2[is.na(res2$stat) & re2s$warning == "list()" & !(res2$test %in% c("conphi", "fl")),]
na_without_warning2 <- res2[is.na(res2$stat)  & !(res2$test %in% c("conphi", "fl")),]
nrow(na_without_warning2) - nrow(na_without_warning2[is.na(na_without_warning2$warning)])

# there are no na when the waning msg 
# "list(list(message = monoblock1 is negative, call = fun(...)), list(message = monoblock2 is negative, call = fun(...)))"
# is saved, thus i will be more precise with the warning than msg just mono1 or mono2 are negative
test <- res2[is.na(res2$stat) & res2$warning %in% c("list(list(message = monoblock1 is negative, call = fun(...)), list(message = monoblock2 is negative, call = fun(...)))"   ),]

#the same but i check the conphi and fl cases because the stat is saved in another column
nrow(res2[is.na(res2$chisq_constrained) & res2$test == "conphi" & res2$warning_constrained == "list()",])
unique(res2$warning_constrained[is.na(res2$chisq_constrained) & res2$test == "conphi"])
# in some cases where chisq constrained is not na, there are warnings that i still consider to be issues
# "some estimated ov variances are negative"
# "Could not compute standard errors!"
unique(res2$warning_constrained[!is.na(res2$chisq_constrained) & res2$test == "conphi"])

nrow(res2[is.na(res2$chisq_unconstrained) & res2$test == "conphi" & res2$warning == "list()",])
unique(res2$warning[is.na(res2$chisq_unconstrained) & res2$test == "conphi"])
# In some cases where chisq unconstrained is not na, there are warnings that I still consider to be issues.
# warnings with  "some estimated ov variances are negative", "Could not compute standard errors!"
# will still be flagged as issues
unique(res2$warning[!is.na(res2$chisq_unconstrained) & res2$test == "conphi"])

# now for fl: 
# i want to check whether there are cases that have a value for ratio1 but not for ratio2 and vice versa 
nrow(res2[!is.na(res2$ave_cor_ratio_1) & is.na(res2$ave_cor_ratio_2),]) 
nrow(res2[is.na(res2$ave_cor_ratio_1) & !is.na(res2$ave_cor_ratio_2),]) 
# no issues regarding that

# now check whether when ratio1 cant be calculated what are the warnigns: 
unique(res2$warning[res2$test == "fl" & is.na(res2$ave_cor_ratio_1)])
unique(res2$warning[res2$test == "fl" & !is.na(res2$ave_cor_ratio_1)])

# fl can be calculated in every case
# still i consider 
# "some estimated ov variances are negative"
# "Could not compute standard errors!"
# Model estimation FAILED!

# further analyis
skim(res[res$test == "fl",])
skim(res[res$test == "conphi",])
skim(res[res$test == "mga",])
skim(res[res$test == "htmt_cov",])
skim(res[res$test == "htmt_cor",])
skim(res[res$test == "htmt_2",])


# i suspect that htmt2 and mga have exactly the same missings
htmt2mga <- res[res$test %in% c("mga", "htmt_2") & is.na(res$stat),]
length(unique(htmt2mga$seed))

htmt2mga2 <- res[res$test %in% c("mga", "htmt_2") & !is.na(res$stat),]

htmt2mga2$res <- (htmt2mga2$stat - htmt2mga2$correlation)^2

htmt2mgacomp <- htmt2mga2 %>% 
                group_by(test, sample_size, correlation, type) %>% 
                summarise(sum = sum(res))

ggplot(htmt2mgacomp, aes(x = factor(sample_size), y = sum, fill = type)) +
  geom_col(position = "dodge") +
  facet_grid(test ~ correlation) +
  scale_y_continuous(limits = c(0, 70000)) + 
  labs(x = "Sample Size", y = "Sum") +
  theme_minimal()


relevant_warnings <- c("list(list(message = monoblock2 is negative, call = fun(...)))",
                       "list(list(message = monoblock1 is negative, call = fun(...)))",
                       "nominator is NaN",
                       "NaNs produced",
                       "some estimated ov variances are negative",
                       "Could not compute standard errors!",
                       "Model estimation FAILED!")


res2$calciss <- grepl(paste(relevant_warnings, collapse = "|"), res2$warning) 

res2$calciss <- res2$calciss | (res2$test %in% c("htmt_cov", "htmt_cor", "htmt_2") & is.na(res2$stat))


issuetable <- res2 %>% 
              group_by(test, sample_size, correlation, type) %>% 
              summarise(issue_rate = mean(calciss))

write.csv2(issuetable, "issuetable.csv")

