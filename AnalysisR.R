##########################################################################
## Analysis script: Permissive construction (working title)             
## April 17, 2021
##
##########################################################################








# Preamble ----------------------------------------------------------------

library(tidyverse)
library(brms)
library(bayestestR)



set.seed(4791)








# Custom functions --------------------------------------------------------

# intraclass correlation

icc <- function(x) {
  v <- (x^2 / ((x^2) + ((pi^2)/3) ) )
  v
  }





# Corpus Analysis ---------------------------------------------------------

# 1. Pre-processing

# 1.1 Read file into R

permis <- read_csv(file = "./Data/dat.csv")

# If the file is in a different folder (e.g., a primary folder), use: 
# Permis <- read_csv(file = "dat.csv")



# 1.2 Drop columns not needed for current analysis

dat <- permis %>% 
  select(!c(starts_with("Match"), starts_with("Tagged"))
         )



# 1.3 Check NAs

dat %>% 
  filter(across(.cols = everything(), 
                .fns = ~ !is.na(.)
                )
         )


# The two columns: BNC_domain and Com_chanel contain all of the 62 NAs 
# All other columns do not contain any NAs

dat %>% 
  filter(is.na(Com_channel) | is.na(BNC_domain)) %>% 
  summarize(n = n())



# 1.4 Extract verbs from the column V2 and store them in V2_re with word()
# By default, word() splits 1st element with a white space word(X, sep = fixed(" "))
# Then, count a string length with str_length()

dat <- dat %>% 
  mutate(V2_re        = word(V2),
         V2_re_length = str_length(V2_re)
         ) %>% 
  relocate(c(V2_re, V2_length, V2_re_length), .after = V2)
  




# 2 Calculate association statistics

# 2.1 Export a data frame with V1 and V2_re columns for collostructional analysis

write_delim(dat %>% 
              select(V1, V2_re) %>% 
              as.data.frame(), 
            file = "./Data/Verb_list.txt",
            delim = "\t"
            )



# 2.2 Read back in results of collostructional analysis (conducted with Gries' R script)

dca <- read_csv("./Data/DCA_results.csv")



# 2.3 Calculate Attraction and Reliance

# 2.3.1 We will calculate these two scores separately for let & permit

dca_let <- dca %>% 
  filter(Prefer == "let" | Prefer == "no_preference")


dca_let <- dca_let %>% 
  mutate(Attract = RawF_let/818,
         Relianc = RawF_let/(RawF_let + RawF_permit)
         ) %>% 
  rowwise() %>% 
  mutate(Minimum = min(Attract, Relianc)
         ) %>% 
  ungroup() %>% 
  mutate(across(.cols = c(Attract,Relianc, Minimum), .fns = ~round(., digits = 3))
         )


dca_permit <- dca %>% 
  filter(Prefer == "permit") 


dca_permit <- dca_permit %>% 
  mutate(Attract = RawF_permit/818,
         Relianc = RawF_permit/(RawF_permit + RawF_let)
         ) %>% 
  rowwise() %>% 
  mutate(Minimum = min(Attract, Relianc)
         ) %>% 
  ungroup() %>% 
  mutate(across(.cols = c(Attract,Relianc, Minimum), .fns = ~round(., digits = 3))
         )


dca <- bind_rows(dca_permit, dca_let)

rm(dca_let, dca_permit)



# 2.4 Take absolute values of Delta P
# We now have scores indicating how words are attracted to their preferred construction

dca <- dca %>% 
  mutate(across(.cols = starts_with("Dp_"), abs)
         )



# 2.5 Join the two tibbles together

dat <- left_join(dat, 
                 dca %>% select(-starts_with(c("Raw", "ExpF", "Prefer"))),
                 by = c("V2_re" = "Words")
                 )

rm(dca)





# 3 Prepare the data frame for regression

# 3.1 Create a new column "Text_rnd" (for varying intercepts)

# 3.1.1 We began by re-coding files keeping those with frequencies of 5 or greater
# while converting files with frequencies < 5 as "OTH" (similar to what Levshina did)

dat <- dat %>% 
  group_by(Text_id) %>% 
  mutate(n_file = n()
         ) %>% 
  mutate(Text_rnd = if_else(n_file < 5, "OTH", Text_id)
         ) %>% 
  ungroup() %>% 
  select(!n_file)


# 3.1.2 We re-coded "OTH" because there were more than 1000 rows with this factor level
# Info for each text_id is available from "textidentifier.csv"
# We created categories based on: http://www.natcorp.ox.ac.uk/docs/URG/bibliog.html

iden <- read_csv(file = "./Data/textidentifier.csv") 


# we split data frame into two smaller ones and recode "dat_oth"

dat_com <- dat %>% 
  filter(Text_rnd != "OTH")

dat_oth <- dat %>% 
  filter(Text_rnd == "OTH")


dat_oth <- left_join(dat_oth, iden, by = "Text_id") %>% 
  mutate(Text_rnd = case_when(Text_types == "academic journal" ~ "ACJ",
                              Text_types == "book"             ~ "BOK", 
                              Text_types == "essays"           ~ "ESY",
                              Text_types == "interview"        ~ "INT",
                              Text_types == "leaflet"          ~ "LFT",
                              Text_types == "lecture/ training"~ "LCT",
                              Text_types == "magazine"         ~ "MAG",
                              Text_types == "medical consultation"        ~ "MED",
                              Text_types == "meeting"         ~ "MET",
                              Text_types == "newspaper"       ~ "NEW",
                              Text_types == "presentation/ public speech" ~ "PRE",
                              Text_types == "seminar/ conference"         ~ "SEM",
                              Text_types %in% c("court hearing", "debate", "other - spoken", 
                                                "tv or radio broadcast")  ~ "OTS",
                              Text_types %in% c("document", "grant abstract", "letters", "newsletter", 
                                                "other - written", "report", "thesis")         ~ "OTW"
                              )
         ) %>% 
  select(!Text_types)


# Then bind the two data frames back together

dat <- bind_rows(dat_com, dat_oth) %>% 
  relocate(Text_rnd, .after = Text_id)


rm(iden, dat_com, dat_oth)



# 3.2 Create a new column "V2_rnd" (for varying intercepts)
# Following Levshina, we coded verbs with freq < 5 as "other"

dat <- dat %>% 
  group_by(V2_re) %>% 
  mutate(n_verb = n()
         ) %>% 
  mutate(V2_rnd = if_else(n_verb < 5, "other", V2_re)
         ) %>% 
  ungroup() %>% 
  relocate(V2_rnd, .after = V2_re) %>% 
  select(!n_verb)





# 4 Bayesian regression modeling

# 4.1 Model 0

# 4.1.1 We began by setting priors with weakly informative prior Normal(0,1) and exponential(1)
# We carried out prior predictive checks at this step as well but code is provided below

mod_prior <- c(set_prior("normal(0, 1)", class = "Intercept"),
               set_prior("exponential(1)", class = "sd")
               )


# 4.1.2 we ran the intercept-only model to assess variability that exists in the data
# We do this by calculating intra-class correlation (ICC)

mod_empt <- brm(V1 ~ 1 + (1|Text_rnd) + (1|V2_rnd),
                data = dat, 
                family = bernoulli(link = "logit"),
                prior  = mod_prior,
                chains = 4,
                iter   = 4000,
                warmup = 1000,
                seed   = 1112,
                cores  = parallel::detectCores() 
                )  


icc(summary(mod_empt)$random[[1]][1] )
icc(summary(mod_empt)$random[[2]][1] )

# ANSWER Intra-class correlation: BNC files = 0.87; Verbs2 = 0.45
# 87% of the chances of "permit" being used is explained by between-file differences
# 45% of the chances of "permit" being used is explained by between-verb differences

rm(mod_empt, mod_prior)



# 4.2 Model 1: Continuous and categorical predictors 

# 4.2.1 Re-code predictors
# V2_permittee_control: convert "unclear" (n = 194) to NA

dat %>% count(V2_permittee_control)

dat <- dat %>% 
  mutate(V2_permittee_control_re = case_when(V2_permittee_control == "no control" ~ "no_con",
                                             V2_permittee_control == "unclear"    ~ NA_character_,
                                             TRUE ~ V2_permittee_control)
         ) %>% 
  relocate(V2_permittee_control_re, .after = V2_permittee_control)


# 4.3.2 Center continuous predictors
# V2_re_length (level-1/text-level predictor): grand-mean center (NOTE: mean = 4.75; range = 2, 13)
# V1_distance (level-1/text-level predictor): grand-mean center (NOTE: mean = 1.92; range = 1, 20)

dat <- dat %>% 
  mutate(V1_distance_c = V1_distance - mean(V1_distance),
         V2_length_c = V2_re_length - mean(V2_re_length)
         ) %>% 
  relocate(V2_length_c, .after = V2_re_length) %>% 
  relocate(V1_distance_c, .after = V1_distance) %>% 
  mutate(across(.cols  = Dp_word_given_con:Minimum,
                .fns   = ~. - mean(.),
                .names = "{.col}_c")
         )


# 4.3.3 Check frequency counts

dat %>% count(V2_permittee_control_re)
dat %>% count(V2_valency) 
dat %>% count(Permittee_semantics)
dat %>% count(Permitter_semantics)
dat %>% count(Com_channel)
dat %>% count(BNC_domain)
dat %>% count(Domain_use) 
dat %>% count(V1_imperative)



#################################### NOTE #################################### 
# Level-one (text/sentence) predictors: 
# V2 = V2_permittee_control_re, V2_valency, V2_length_c, association scores 
# V1 = V1_imperative, V1_distance (from V1 to V2)                          
# N = Permittee_semantics, Permitter_semantics

# Level-two (file) predictors
# Domain_use, Com_channel
############################################################################## 



# 4.3.4 Set-up contrasts

dat <- dat %>% 
  mutate(across(.cols = c(V2_permittee_control_re, 
                          V2_valency, 
                          V1_imperative,
                          Permittee_semantics, 
                          Permitter_semantics,
                          Domain_use, 
                          Com_channel),
                .fns = as.factor)
         )


contrasts(dat$V2_permittee_control_re) <- contr.sum(2)
contrasts(dat$V2_valency) <- contr.sum(2)
contrasts(dat$V1_imperative) <- contr.sum(2)
contrasts(dat$Permittee_semantics) <- contr.sum(2)
contrasts(dat$Permitter_semantics) <- contr.sum(2)
contrasts(dat$Domain_use) <- contr.sum(2)
contrasts(dat$Com_channel) <- contr.sum(2)


# 4.3.5 Specify priors for the model
  
mod_prior <- c(set_prior("normal(0,1)", class = "Intercept"),
               set_prior("normal(0,1)", class = "b"),
               set_prior("exponential(1)", class = "sd")
               )


# 4.3.6 Conduct prior predictive checks with sample_prior = "only" in code
# Then, extract posterior samples from the model and plot the samples to see 
# what the model expects from the priors

# We used some of the predictors to run this model (continuous and categorical)

priorcheck <- brm(V1 ~ 1 + V2_permittee_control_re + V2_length_c + V1_imperative + 
                    V1_distance_c + Domain_use + (1|Text_rnd) + (1|V2_rnd),
                  data   = dat, 
                  family = bernoulli(link = "logit"),
                  prior  = mod_prior,
                  sample_prior = "only",
                  chains = 4,
                  iter   = 4000,
                  warmup = 1000,
                  seed   = 191,
                  cores  = parallel::detectCores()
                  )


sam_priorchk <- posterior_samples(priorcheck, 
                                  pars = c("^b_", "^sd_"),
                                  add_chain = TRUE)


sam_priorchk %>%
  select(!c(chain, iter)) %>% 
  mutate(across(.cols = everything(), 
                .fns  = brms::inv_logit_scaled
                ) 
         ) %>% 
  ggplot() + 
  # geom_density(aes(x = sd_Text_rnd__Intercept) )
  # geom_density(aes(x = sd_V2_rnd__Intercept) )
  # geom_density(aes(x = b_Intercept) )
  # geom_density(aes(x = b_V1_imperative1) )
  # geom_density(aes(x = b_V1_distance_c) )
  geom_density(aes(x = b_V2_length_c) )


# We can see, for example, the "fixed effect" term for the intercept ranges from probability of 0 to 1
# with the peak lying around 0.5, indicating that before seeing data, the model thinks the probability of 
# 0.5 is the most plausible

rm(priorcheck, sam_priorchk)



# 4.3.7 Run models with all predictors + one association measure

mod1 <- brm(V1 ~ 1 + V2_permittee_control_re + V2_valency + V2_length_c + V1_imperative +
              V1_distance_c + Permittee_semantics + Permitter_semantics + Dp_word_given_con_c +
              Domain_use + Com_channel + (1|Text_rnd) + (1|V2_rnd),
            data   = dat,
            family = bernoulli(link = "logit"),
            prior  = mod_prior,
            chains = 4,
            iter   = 4000,
            warmup = 1000,
            seed   = 1150,
            cores  = parallel::detectCores(),
            control = list(adapt_delta = 0.99),
            save_pars = save_pars(all = TRUE)
            ) 

mod2 <- update(mod1,
               formula. = ~. - Dp_word_given_con_c + Dp_con_given_word_c,
               newdata  = dat 
               )

mod3 <- update(mod1,
               formula. = ~. - Dp_word_given_con_c + Coll_strength_c,
               newdata  = dat
               )

mod4 <- update(mod1, 
               formula. = ~. - Dp_word_given_con_c + Attract_c,
               newdata  = dat
               )

mod5 <- update(mod1,
               formula. = ~. - Dp_word_given_con_c + Relianc_c,
               newdata  = dat
               )

mod6 <- update(mod1, 
               formula. = ~. - Dp_word_given_con_c + Minimum_c,
               newdata  = dat
               )


# 4.3.8 Add LOO-IC estimates to each model and compare all models

l_mod1 <- add_criterion(mod1, criterion = c("loo", "waic"), moment_match = TRUE)
l_mod2 <- add_criterion(mod2, criterion = c("loo", "waic"), moment_match = TRUE)
l_mod3 <- add_criterion(mod3, criterion = c("loo", "waic"), moment_match = TRUE)
l_mod4 <- add_criterion(mod4, criterion = c("loo", "waic"), moment_match = TRUE)
l_mod5 <- add_criterion(mod5, criterion = c("loo", "waic"), moment_match = TRUE)
l_mod6 <- add_criterion(mod6, criterion = c("loo", "waic"), moment_match = TRUE)

loo_compare(l_mod1, l_mod2, l_mod3, l_mod4, l_mod5, l_mod6, criterion = "loo")
loo_compare(l_mod1, l_mod2, l_mod3, l_mod4, l_mod5, l_mod6, criterion = "waic")


# ANSWER: It was found that model with collostructional strength had the smallest
# LOOIC score (though the estimates of it overlapped with zero)

rm(mod1, mod2, mod3, mod4, mod5, mod6)


# 4.3.9 Rerun the model with collo. strength as the association measure

modf <- brm(V1 ~ 1 + V2_permittee_control_re + V2_valency + V2_length_c + V1_imperative +
              V1_distance_c + Permittee_semantics + Permitter_semantics + Coll_strength_c +
              Domain_use + Com_channel + (1|Text_rnd) + (1|V2_rnd),
            data   = dat,
            family = bernoulli(link = "logit"),
            prior  = mod_prior,
            chains = 4,
            iter   = 4000,
            warmup = 1000,
            seed   = 1619,
            cores  = parallel::detectCores(),
            control = list(adapt_delta = 0.99),
            save_pars = save_pars(all = TRUE)
            ) 


summary(modf)$fixed
summary(modf)$random[[1]]
summary(modf)$random[[2]]

icc(summary(modf)$random[[1]][1])
icc(summary(modf)$random[[2]][1])


# 4.3.10 Rerun the model without Permittee_semantics
# Estimates overlapped with zero; thus, the predictor was dropped
# Estimates of collostructional strength also overlapped with zero. Despite this,
# we kept the predictor in since we're interested in how this info predicted human response

modr <- update(modf, 
               formula. = ~. - Permittee_semantics,
               newdata  = dat
               )


summary(modr)$fixed


# 4.3.11 Compare the two models

l_modf <- add_criterion(modf, "loo", moment_match = TRUE)
l_modr <- add_criterion(modr, "loo", moment_match = TRUE)

loo_compare(l_modf, l_modr)


# 4.3.12 Make predictions of the response values < logit(theta) > 

# We used predict() to obtain predicted values of the response variable  
# taking into account variance components and residuals included in the model
# e.g., sd_Intercept_V2, sd_Intercept_file, etc.

# We filtered out NA from data and bind columns from these two tibbles

dat_est <- bind_cols(dat %>% 
                       filter(across(.cols = everything(),
                                     .fns  = ~!is.na(.)
                                     )
                              ),
                     predict(modr, summary = T) %>% 
                       as_tibble() %>% 
                       select(Estimate)
                     )


# We binned probabilities into 5 distinct groups

dat_est <- dat_est %>% 
  mutate(Estimate = round(Estimate, digits = 3)
         ) %>% 
  mutate(across(.cols = Estimate,
                .fns  = ~case_when(. <= 0.199             ~ "Bin_1",
                                   . >=  0.2 & . <= 0.399 ~ "Bin_2",
                                   . >=  0.4 & . <= 0.599 ~ "Bin_3",
                                   . >=  0.6 & . <= 0.799 ~ "Bin_4",
                                   . >=  0.8              ~ "Bin_5"
                                   ),
                .names = "{.col}_b"
                  )
         ) 








# 4.3.13 Randomly select sentences

set.seed(122021)

stimuli <- dat_est %>% 
  filter(V2_type == "verb") %>% 
  filter(Estimate_b %in% c("Bin_2", "Bin_3")) %>% 
  group_by(Estimate_b, V1_imperative, V2_permittee_control_re, V2_valency, 
           Permitter_semantics, Com_channel, Domain_use) %>% 
  slice_sample(n = 4, replace = FALSE) %>% 
  ungroup() %>% 
  mutate(No = row_number()) %>% 
  relocate(No)

stimuli <- stimuli %>% 
  arrange(Estimate_b, V2_permittee_control_re, V2_valency, Permitter_semantics, Com_channel, Domain_use, V1_imperative)

stimuli %>%
  ggplot() + 
  geom_point(aes(x = No, y = Estimate)
             )


write_csv(stimuli, file = "bin2_3_moreexp.csv")






