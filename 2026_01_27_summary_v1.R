library(lavaan)
library(semTools)
library(doParallel)
library(foreach)
library(dplyr)
library(ggplot2)
library(skimr)
library(stringr)

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
na_without_warning <- res2[is.na(res2$stat) & res2$warning == "list()" & !(res2$test %in% c("conphi", "fl")),]
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
skim(res[res$test == "mga" & res$correlation == 1 & res$type == "harsh" & res$sample_size == 50,])
skim(res[res$test == "htmt_cov",])
skim(res[res$test == "htmt_cor",])
skim(res[res$test == "htmt_2" & res$correlation == 1 & res$type == "harsh" & res$sample_size == 50, ])


# i suspect that htmt2 and mga have exactly the same missings
htmt2mga <- res[res$test %in% c("mga", "htmt_2") & is.na(res$stat),]
length(unique(htmt2mga$seed))

htmt2mga2 <- res[res$test %in% c("mga", "htmt_2") & !is.na(res$stat),]

htmt2mga2$res <- (htmt2mga2$stat - htmt2mga2$correlation)^2

htmt2mgacomp <- htmt2mga2 %>% 
                group_by(test, sample_size, correlation, type) %>% 
                summarise(sum = mean(res)) 

ggplot(htmt2mgacomp, aes(x = factor(sample_size), y = sum, fill = type)) +
  geom_col(position = "dodge") +
  facet_grid(test ~ correlation) +
  labs(x = "Sample Size", y = "Sum") +
  theme_minimal()

# i want to plot the errors so im gonna select all cases where the warning is an issue 
relevant_warnings <- c("list(list(message = monoblock2 is negative, call = fun(...)))",
                       "list(list(message = monoblock1 is negative, call = fun(...)))",
                       "nominator is NaN",
                       "NaNs produced",
                       "some estimated ov variances are negative",
                       "Could not compute standard errors!",
                       "Model estimation FAILED!")

# grepl 
res2$calciss <- grepl(paste(relevant_warnings, collapse = "|"), res2$warning) 
res2$calciss_constrained <- grepl(paste(relevant_warnings, collapse = "|"), res2$warning_constrained)

# i see that i dont get all cases where htmts are na so im gonna add an OR statement 
res2$calciss <- res2$calciss | (res2$test %in% c("htmt_cov", "htmt_cor", "htmt_2") & is.na(res2$stat))
res2$calciss <- res2$calciss | res2$calciss_constrained

issuetable <- res2 %>% 
              group_by(test, sample_size, correlation, type) %>% 
              summarise(issue_rate = mean(calciss))


res2 %>%
  group_by(test, sample_size, correlation, type) %>%
  summarise(issue_rate = mean(calciss), .groups = "drop") %>%
  ggplot(aes(x = factor(sample_size), y = factor(correlation), fill = issue_rate)) +
  geom_tile() +
  geom_text(aes(label = round(issue_rate, 3)), color = "white") +
  facet_grid(test ~ type) +
  scale_fill_viridis_c() +
  theme_minimal()

# it seems that htmt_cor and htmt_cov have no issues whatsoever in mild conditions 
htmt_cor_mild <- res2[res2$test == "htmt_cor" & res2$type == "mild",]
htmt_cov_mild <- res2[res2$test == "htmt_cov" & res2$type == "mild",]

skim(htmt_cor_mild)
skim(htmt_cov_mild)
# indeed there are no missing values in mild conditions
# but are there any warning msg? 
htmt_cor_mild_warning <- unique(htmt_cor_mild$warning)
# no warnings for cor
htmt_cov_mild_warning <- unique(htmt_cov_mild$warning)
# no warnings for cov

unique(res2$warning[res2$calciss == TRUE & res2$stat == "conphi"])


res2$realwarning <- NA
res2$realwarning[res2$calciss == TRUE] <- res2$warning[res2$calciss == TRUE] 

res2_freq <- res2 %>%
  count(
    type,
    sample_size,
    correlation,
    test,
    realwarning,
    name = "n"
  ) %>%
  group_by(type, sample_size, correlation, test) %>%
  mutate(rel_freq = n / sum(n)) %>%
  ungroup()

res2_freq$warning_msg <- NA
res2_freq$warning_msg[res2_freq$realwarning %in% c("list()", NA)] <- "No warning"
res2_freq$warning_msg[res2_freq$realwarning == "list(list(message = lavaan->lav_lavaan_step11_estoptim():  \\n   Model estimation FAILED! Returning starting values., call = NULL))"] <- "Model_estimation_failed" 
res2_freq$warning_msg[res2_freq$realwarning == "list(list(message = lavaan->lav_model_vcov():  \\n   Could not compute standard errors! The information matrix could not be \\n   inverted. This may be a symptom that the model is not identified., call = NULL), list(message = lavaan->lav_object_post_check():  \\n   some estimated ov variances are negative, call = NULL))"] <- "could not compute standard errors, some estimated ov variances are negative"
res2_freq$warning_msg[res2_freq$realwarning == "list(list(message = lavaan->lav_object_post_check():  \\n   some estimated ov variances are negative, call = NULL))"] <- "some ov variances are negative"
res2_freq$warning_msg[res2_freq$realwarning == "list(list(message = monoblock1 is negative, call = fun(...)))"] <- "monoblock1 is negative"
res2_freq$warning_msg[res2_freq$realwarning == "list(list(message = monoblock2 is negative, call = fun(...)))"]<- "monoblock2 is negative"
res2_freq$warning_msg[res2_freq$realwarning == "list(list(message = NaNs produced, call = sqrt(diag(Lzero))))" ] <- "variances are negative"
res2_freq$warning_msg[res2_freq$realwarning == "list(list(message = lavaan->lav_object_post_check():  \\n   covariance matrix of latent variables is not positive definite ; use \\n   lavInspect(fit, \\cov.lv\\) to investigate., call = NULL))"] <- "No warning"
res2_freq$warning_msg[res2_freq$realwarning == "list(list(message = lavaan->lav_object_post_check():  \\n   some estimated ov variances are negative, call = NULL), list(message = lavaan->lav_object_post_check():  \\n   covariance matrix of latent variables is not positive definite ; use \\n   lavInspect(fit, \\cov.lv\\) to investigate., call = NULL))"] <- "some ov variances are negative"


res2_freq_only_warnings <- res2_freq
res2_freq_only_warnings <- res2_freq_only_warnings[res2_freq_only_warnings$warning_msg != "No warning" , ]

res2_freq_only_warnings_fl_conphi <-  res2_freq_only_warnings[res2_freq_only_warnings$test %in% c("conphi", "fl"),]
res2_freq_only_warnings_htmt_mga <-  res2_freq_only_warnings[res2_freq_only_warnings$test %in% c("mga", "htmt_cov", "htmt_cor", "htmt_2"),]


plot <- res2_freq_only_warnings_fl_conphi %>%
  ggplot(aes(y = warning_msg, x = rel_freq, fill = test)) +
  geom_col(position = "dodge") +
  facet_grid(
    type + sample_size ~ correlation,
    labeller = labeller(
      sample_size = ~paste("n =", .),
      correlation = ~paste("r =", .)
    )
  ) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Relative Frequency of Warning Messages",
    x = "Relative Frequency",
    y = "Warning Message",
    fill = "Method"
  ) +
  theme_minimal()


plot2 <- res2_freq_only_warnings_htmt_mga %>%
  ggplot(aes(y = warning_msg, x = rel_freq, fill = test)) +
  geom_col(position = "dodge") +
  facet_grid(
    type + sample_size ~ correlation,
    labeller = labeller(
      sample_size = ~paste("n =", .),
      correlation = ~paste("r =", .)
    )
  ) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Relative Frequency of Warning Messages",
    x = "Relative Frequency",
    y = "Warning Message",
    fill = "Method"
  ) +
  theme_minimal()

