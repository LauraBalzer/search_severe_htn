##############
# Code to reproduce effectiveness analyses in 
# "SEARCH integrated HIV/hypertension community health worker-led intervention in rural East Africa"
#
# Laura B. Balzer, PhD MPhil
# laura.balzer@berkeley.edu
# Lead Statistician for SEARCH
##############

rm(list=ls())

library('SuperLearner')
source('estimation.R')
source('stage2.R')
source('tmle.R')
source('aps.R')

# read in data and process for estimation code
df <- read.csv('search_severe_htn_shared.csv')
df$A <- df$studyarm
df$id <- 1:nrow(df) # unique id
df$U <- df$alpha <- 1 # dummy variables used in TMLE

# specify outcome variables
outcome_var <- c('htncontrol_24', 'sevhtn_24','retained_24', 'sbp_24',
                 'htncontrol_48', 'sevhtn_48',  'retained_48',  'sbp_48' ,'tic')

OUT <- NULL

# for each outcome variable, run the primary effectiveness analysis (TMLE) + unadjusted effect estimator
# see the Statistical Analysis Plan (SAP) for further details

for(j in outcome_var){
  set.seed(1)
  this_outcome <- j
  print(this_outcome)
  temp <- run_analyses(O=df, outcome_var=this_outcome)
  write.csv(temp, row.names=F, file=paste0(this_outcome, '.csv'))
  temp <- data.frame(cbind(outcome=this_outcome, temp$est))
  OUT <- rbind(OUT, temp)
}

write.csv(OUT, file='combined_results.csv')

# Table 2 in paper
OUT[OUT$Approach=='Primary', c("outcome","INT.est", "INT.CI.lo","INT.CI.hi", 
                               "CON.est", "CON.CI.lo", "CON.CI.hi", 
                               "EFFECT.est", "EFFECT.CI.lo", "EFFECT.CI.hi","EFFECT.pval")]
                               