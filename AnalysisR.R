





# Preamble ----------------------------------------------------------------

library(tidyverse)
library(brms)


set.seed(4791)





# Preamble ----------------------------------------------------------------

# intraclass correlation

icc <- function(x) {
  v <- (x^2 / ((x^2) + ((pi^2)/3) ) )
  v
  }





# Corpus Analysis ---------------------------------------------------------

# 1. Pre-processing

# 1.1 Read file into R

permis <- read_csv(file = "./Data/dat.csv")

# If you don't store this file in "Data" folder but instead in main folder, use: 
# Permis <- read_csv(file = "dat.csv")


# 1.2 Drop columns not needed for current analysis

dat <- permis %>% 
  select(!starts_with("Match") & !starts_with("Tagged"))


# 1.3 Check various counts

dat %>% 
  filter(across(.cols = everything(), 
                .fns = ~ !is.na(.)
                )
         )


dat %>% 
  count(V1)

dat %>% 
  count(V2_type)

dat %>% 
  count(V2_permittee_control)



# 2 Calculate association statistics

# 2.1 Export a data frame with V1 and V2 columns for collostructional analysis

write_delim(dat %>% 
              select(V1, V2) %>% 
              as.data.frame(), 
            file = "./Data/Verb_list.txt",
            delim = "\t"
            )


# 2.2 Read into R the results of collostructional analysis

dca <- read_csv("./Data/DCA_results.csv")



# 2.3 Join the two tibbles together

Permis <- left_join(Permis, 
                    dca %>% 
                      select(-starts_with("Raw"), -starts_with("Exp")) %>% 
                      rename(Delta_p_word_given_con = Delta_p_word_given_constr1,
                             Delta_p_con_given_word = Delta_p_constr1_given_word,
                             V2_preference          = Preference) %>% 
                      mutate(across(where(is.numeric), abs)),
                    by = c("V2" = "Words")
                    )




#Obtain frequency counts

# NOTE
# (1) Pivot a count table into a wide format with columns of counts for "permit" and "let"
# (2) Convert NA into 0 with replace_na() inside across(); we use . inside replace_na() to
# reference current columns; thus, ~replace_na(., 0)

FreqTab <- Permis %>% 
  group_by(V2) %>% 
  count(V1) %>% 
  pivot_wider(names_from = V1, values_from = n) %>% 
  mutate(across(where(is.numeric), ~replace_na(., 0)
                )
         )

a <- FreqTab %>% 
  mutate(permit_y = permit > let,
         prefer   = if_else(permit_y == TRUE, "permit", "let")
         )



TotPermit <- sum(FreqTab$permit)

FreqTab %>% 
  mutate(a = permit, 
         b = let,
         c = TotPermit - a) %>% 
  mutate(Attraction = a/(a+c),
         Reliance = a/(a+b)) %>% 
  select(V2, Attraction, Reliance)








# 3 Prepare the data frame for regression

# 3.1 Crate a new column for file names "Text_rnd" and for verbs "V2_rnd"
# Following Levshina, we converted files with fewer than 5 occurrences to "OTH" (other)
# We also recoded verbs with fewer than 5 occurrences to "other"

dat <- dat %>% 
  group_by(Text_id) %>% 
  mutate(n_file = n()
         ) %>% 
  mutate(Text_rnd = if_else(n_file < 5, "OTH", Text_id)
         ) %>% 
  ungroup() %>% 
  relocate(Text_rnd, .after = Text_id) %>% 
  group_by(V2) %>% 
  mutate(n_verb = n()
         ) %>% 
  mutate(V2_rnd = if_else(n_verb < 5, "other", V2)
         ) %>% 
  ungroup() %>% 
  relocate(V2_rnd, .after = V2) %>% 
  select(!c(n_file, n_verb))



# 3.2 Re-code predictors
# V2_type: combine "idiom" and "phrasal v" into one category called "verb_phrase"
# V2_permittee_control: convert "unclear" to NA

dat <- dat %>% 
  mutate(V2_type_re = case_when(V2_type == "idiom"        ~ "verb_phrase",
                                V2_type == "phrasal verb" ~ "verb_phrase",
                                TRUE ~ V2_type),
         V2_permittee_control_re = case_when(V2_permittee_control == "no control" ~ "no_con",
                                             V2_permittee_control == "unclear"    ~ NA_character_,
                                             TRUE ~ V2_permittee_control)
         ) %>% 
  relocate(V2_type_re, .after = V2_type) %>% 
  relocate(V2_permittee_control_re, .after = V2_permittee_control)



# 3.3 Center continuous variable
# V2_length (level-1/text-level predictor): grand-mean center (NOTE: mean V2 length = 5.3)

dat <- dat %>% 
  mutate(V2_length_c = V2_length - mean(V2_length))




dat %>% count(V2_type_re)
dat %>% count(V2_permittee_control_re)
dat %>% count(V2_valency) 
dat %>% count(Permittee_semantics)
dat %>% count(Permitter_semantics)
dat %>% count(Com_channel)
dat %>% count(BNC_domain)
dat %>% count(Domain_use) #HERE
dat %>% count(V1_imperative)









# Model 0: Empty model
# We begin by setting priors with weakly informative prior Normal(0,1) and exponential(1)

mod_prior <- c(set_prior("normal(0, 1)", class = "Intercept"),
               set_prior("exponential(1)", class = "sd")
               )


# We then conduct prior predictive check with  sample_prior = "only" in code
# Then, we extract posterior samples (intercept and varying intercept components) 
# We plot the samples to see samples generated solely from the prior

PriorCheck_ModEmpt <- brm(V1 ~ 1 + (1|Text_rnd) + (1|V2_rnd),
                          data = dat, 
                          family = bernoulli(link = "logit"),
                          prior  = mod_prior,
                          sample_prior = "only",
                          chains = 4,
                          iter   = 4000,
                          warmup = 1000,
                          seed   = 1112,
                          cores  = parallel::detectCores() 
                          )


md_priorch <- posterior_samples(PriorCheck_ModEmpt, 
                                pars = c("^b_", "^sd_"),
                                add_chain = TRUE)


md_priorch %>%
  select(!c(chain, iter)) %>% 
  mutate(across(.cols = everything(), 
                .fns  = brms::inv_logit_scaled
                ) ) %>% 
  ggplot() + 
  geom_density(aes(x = sd_Text_rnd__Intercept)) #+
  geom_density(aes(x = b_Intercept))
    
# We can see, for example, the "fixed" effect term for the intercept ranges from probability of 0 to 1
# with the peak lying around 0.5, indicating that before seeing data, our priors think the probability of 
# 0.5 is the most plausible

  
# Finally, we run the intercept-only model to assess variability that exists in the data
# We do this by calculating intra-class correlation (ICC)
  

Mod_empt <- brm(V1 ~ 1 + (1|Text_rnd) + (1|V2_rnd),
                data = dat, 
                family = bernoulli(link = "logit"),
                prior  = mod_prior,
                chains = 4,
                iter   = 4000,
                warmup = 1000,
                seed   = 1112,
                cores  = parallel::detectCores() 
                )  


icc(summary(Mod_empt)$random[[1]][1] )
icc(summary(Mod_empt)$random[[2]][1] )

# ANSWER Intra-class correlation: BNC files = 0.91; Verbs2 = 0.49
# 91% of the chances of "permit" being used is explained by between-file differences
# 49% of the chances of "permit" being used is explained by between-verb differences

# We then fit a model with predictors to see if they can explain this variablity




# Level-one (text/sentence) predictors: 
# V2 = V2_type_re, V2_permittee_control_re, V2_valency, V2_length_c
# V1 = V1_imperative
# Permittee_semantics, Permitter_semantics

# Level-two (file) predictors
# Domain_use, Com_channel



# Set-up contrasts

dat <- dat %>% 
  mutate(across(.cols = c(V2_type_re, V2_permittee_control_re, V2_valency,
                          V1_imperative, Permittee_semantics, Permitter_semantics,
                          Domain_use, Com_channel),
                .fns = as.factor)
         )


contrasts(dat$V2_type_re) <- contr.sum(2)
contrasts(dat$V2_permittee_control_re) <- contr.sum(2)
contrasts(dat$V2_valency) <- contr.sum(2)
contrasts(dat$V1_imperative) <- contr.sum(2)
contrasts(dat$Permittee_semantics) <- contr.sum(2)
contrasts(dat$Permitter_semantics) <- contr.sum(2)
contrasts(dat$Domain_use) <- contr.sum(2)
contrasts(dat$Com_channel) <- contr.sum(2)




# Specify priors for the model
  
mod1_prior <- c(set_prior("normal(0,1)", class = "Intercept"),
                set_prior("normal(0,1)", class = "b"),
                set_prior("exponential(1)", class = "sd")
                )


# Run the model

mod <- brm(V1 ~ V2_type_re + V2_permittee_control_re + V2_valency + V2_length_c  +
             V1_imperative + Permittee_semantics + Permitter_semantics +
             Domain_use + Com_channel + (1|Text_rnd) + (1|V2_rnd),
             data   = dat,
             family = bernoulli(link = "logit"),
             prior  = mod1_prior,
             chains = 4,
             iter   = 4000,
             warmup = 1000,
             seed   = 1331,
             cores  = parallel::detectCores(),
             control = list(adapt_delta = 0.99)
             ) 


summary(mod)$fixed



# Mutate new columns with estimates from the model 
# Then, sum all these estimates up and convert the sum to probability

dat_est <- dat %>% 
  mutate(Est_Int             = -2.3629171,
         Est_V2_type         = if_else(V2_type_re == "verb", 0.2624965, -0.2624965
                                       ),
         Est_V2_permitee_con = case_when(V2_permittee_control_re == "control" ~ 0.8788518,
                                         V2_permittee_control_re == "no_con"  ~ -0.8788518,
                                         is.na(V2_permittee_control_re) ~ NA_real_
                                         ),
         Est_V2_valency      = if_else(V2_valency == "active", -1.2648889, 1.2648889
                                       ),
         Est_V1_imperative   = if_else(V1_imperative == "imperative", -2.1841174, 2.1841174
                                       ),
         Est_Permittee_seman = if_else(Permittee_semantics == "animate", -0.2186778, 0.2186778
                                       ),
         Est_Permitter_seman = if_else(Permitter_semantics == "animate", -1.6444657, 1.6444657
                                       ),
         Est_Domain_use      = case_when(Domain_use == "personal" ~ -0.4653962, 
                                         Domain_use == "public" ~ 0.4653962,
                                         is.na(Domain_use) ~ NA_real_
                                         ),
         Est_Com_channel     = if_else(Com_channel == "spoken", -1.8841794, 1.8841794),
         Est_V2_length       = case_when(V2_length == 2 ~ 0.1212567 * -3,
                                         V2_length == 3 ~ 0.1212567 * -2,
                                         V2_length == 4 ~ 0.1212567 * -1,
                                         V2_length == 5 ~ 0.1212567 * 0,
                                         V2_length == 6 ~ 0.1212567 * 1,
                                         V2_length == 7 ~ 0.1212567 * 2,
                                         V2_length == 8 ~ 0.1212567 * 3,
                                         V2_length == 9 ~ 0.1212567 * 4,
                                         V2_length == 10 ~ 0.1212567 * 5,
                                         V2_length == 11 ~ 0.1212567 * 6,
                                         V2_length == 12 ~ 0.1212567 * 7,
                                         V2_length == 13 ~ 0.1212567 * 8,
                                         V2_length == 14 ~ 0.1212567 * 9,
                                         V2_length == 15 ~ 0.1212567 * 10,
                                         V2_length == 21 ~ 0.1212567 * 16
                                         )
         )

dat_est <- dat_est %>%
  rowwise() %>% 
  mutate(Est_sum = sum(c_across(cols = starts_with("Est")), na.rm = TRUE),
         Est_prob = brms::inv_logit_scaled(Est_sum)
         ) %>% 
  mutate(Est_prob = round(Est_prob, digits = 3))



