run_analyses <- function(O, outcome_var='htncontrol_24', goals='RD',
                        verbose=F){
  # get settings 
  SETTINGS <- get_settings(outcome_var=outcome_var)
  
  Ycs <- rownames(SETTINGS$settings)
  
  # now actually do estimation 
  OUT <- NULL 
  
  for(k in 1:length(Ycs)){
    
    this_Yc <- Ycs[k]
    # if small subgroup than always do unadjusted
    this_Yc_temp <- ifelse(nrow(O)<51,  'Unadj', this_Yc)
    
    print(paste0('--------------------', this_Yc_temp, ' Analysis--------------------'))
    
    W <- c('age', 'age_ge60', 'male', 'sbp_enroll', 'grade2_enroll','grade3_enroll', 'kenya')
    
    dt <- make_O(O=O, this_setting=SETTINGS$settings[this_Yc_temp,],
                 verbose=verbose)[,c('id', 'alpha', 'U', W, 'A',  'Delta','Y')]
    
    x <- do_estimation_htn(dt=dt, SETTINGS=SETTINGS,
                           this_Yc=this_Yc_temp, goals=goals, 
                           verbose=verbose)
    
    # overwrite with original name
    x$Approach <- this_Yc
    OUT <- rbind(OUT,x ) 
    
  }
  
  list(est=OUT, SETTINGS=SETTINGS)
  
}




get_settings <- function(outcome_var){
  
  W <- c('age', 'age_ge60', 'male', 'grade2_enroll','grade3_enroll', 'kenya')
  
  # adaptive prespecification candidate algorithms for the fancy version 
  cand.Qform.fancy <- c('glm', 'step') 
  cand.gform.fancy <- c('glm', 'step')
  
  # SuperLearner library for missing data 
  SL.library <- list( 'SL.mean','SL.glm','SL.step')
  
  # CROSS-VALIDATION
  V=10 # unless fewer than 40 indpt units, default to LOOV
  
  # OUTCOME = Hypertension control 
  if(outcome_var %in% c('htncontrol_24', 'htncontrol_48') ){
    settings=data.frame(rbind(
      Pri = cbind( missing='failure', effect='tmle'),
      Unadj = cbind(missing='failure', effect='unadj')
    ))
    rownames(settings) <- c('Primary', 'Unadj')
    
    # OUTCOME = moderate-to-severity 
  } else if(outcome_var %in% c('sevhtn_24', 'sevhtn_48')){
    settings=data.frame(rbind(
      TMLE = cbind( missing='adjust', effect='tmle'),
      Unadj = cbind(missing='exclude', effect='unadj')
    ))
    rownames(settings) <- c('Primary', 'Unadj')
    
    # OUTCOME = retention OR OR time-in-care
  }else if(outcome_var %in% c('retained_24', 'retained_48','tic')){
    # no missingness; no TMLE  for adjustment 
    settings=data.frame(rbind(
      Pri = cbind( missing='NA', effect='tmle'),
      Unadj = cbind(missing='NA', effect='unadj')
    ))
    rownames(settings) <- c('Primary', 'Unadj')
    
    # OUTCOME=average systolic 
  } else if(outcome_var %in% c('sbp_24', 'sbp_48')){
    settings=data.frame(rbind(
      Pri = cbind( missing='impute', effect='tmle', impute.var='sbp_enroll'),
      Unadj = cbind(missing='impute', effect='unadj', impute.var='sbp_enroll')
    ))
    rownames(settings) <- c('Primary', 'Unadj')
  }
  
  settings <- cbind(outcome_var=outcome_var, settings)
  
  list(settings=settings, W=W, V=V, 
       cand.Qform.fancy=cand.Qform.fancy, cand.gform.fancy=cand.gform.fancy,
       SL.library=SL.library, simple.AP=F,
       one.sided=F, alt.smaller=NA)
}


make_O <- function(O, this_setting,verbose=F){
  
  if(verbose) print( paste0('# in initial dataset: ', nrow(O)) )
 
  
  # outcome variable 
  Y <- O[,this_setting$outcome_var]
  
  if(this_setting$missing=='failure'){
    # count any missing as failures (e.g., unsuppressed)
    if(verbose) print(paste0('Counting ', sum(is.na(Y)), ' missing is failure'))
    Y[is.na(Y)] <- 0
  } else if (this_setting$missing=='impute'){
    # impute with selected variable
    if(verbose) print(paste0('Imputing ', sum(is.na(Y)) ))
    Y[is.na(Y)] <- O[is.na(Y), this_setting$impute.var]
  }
  
  # assumed to be 0 (not measured) unless evidence otherwise
  Delta <- rep(0, nrow(O) )	
  # measured endpoint
  Delta[!is.na(Y)] <- 1
  
  # code Delta as Censoring
  Delta <- BinaryToCensoring(is.uncensored=Delta)
  
  dt <- data.frame(O, Delta, Y)
  
  if(this_setting$missing=='exclude'){
    # subset dt to be on measured only
    if(verbose) print(paste0('Excluding ', sum(dt$Delta=='censored'), ' missing endpt'))
    dt <- dt[dt$Delta=='uncensored', ]
  }                        
  
  dt
}



do_estimation_htn <- function(dt, SETTINGS, this_Yc, goals='aRR', verbose=F){
  
  these.cols <- c(paste0('Txt', c('.est', '.CI.lo', '.CI.hi')),
                  paste0('Con', c('.est', '.CI.lo', '.CI.hi')), 
                  c('est', 'CI.lo', 'CI.hi', 'se', 'pval', 'QAdj', 'gAdj'))
  
  
  if(sum(dt$Delta=='uncensored')==nrow(dt)){
    
    # if no missingness, then run stage2 
    if(verbose) print('No missingness')
    
    adj <-  QAdj <- gform <- NULL
    do.data.adapt <- F
    do.unadjusted <- F
    
    if(SETTINGS$settings[this_Yc, 'effect']=='unadj') {
      # unadjusted effect estimator 
      if(verbose) print('Unadjusted effect estimator')
      do.unadjusted <- T
      
    } else if(SETTINGS$settings[this_Yc, 'effect']=='grade'){
      if(verbose) print('Forced adjustment for grade')
      QAdj <- c('grade2_enroll','grade3_enroll')
      
    } else if(SETTINGS$settings[this_Yc, 'effect']=='sbp'){
      if(verbose) print('Forced adjustment for sbp')
      QAdj <- c('sbp_enroll')
      
    } else{
      # adjusted effect estimator
      if(verbose) print('Adjusted effect estimator')
      do.data.adapt <- T
      do.unadjusted <- F
      if(nrow(dt)<41 | SETTINGS$simple.AP){
        # if N<= 40, then do simple adaptive prespec to search among GLMS
        if(verbose) print('Doing simple AP')
        adj <- get.cand.adj(all.cand=SETTINGS$W, cand.Qform.fancy=NULL, cand.gform.fancy=NULL)
      }else{
        # running fancy version with other ML methods
        adj <- get.cand.adj(all.cand=SETTINGS$W, 
                            cand.Qform.fancy=SETTINGS$cand.Qform.fancy, 
                            cand.gform.fancy=SETTINGS$cand.gform.fancy)
      }
    }
    
    est <- Stage2(goal=goals, data.input=dt, target='indv', 
                  do.unadjusted=do.unadjusted,
                  QAdj=QAdj, gAdj=QAdj, 
                  do.data.adapt=do.data.adapt, 
                  remove.pscore=F,
                  cand.QAdj=adj$cand.QAdj, cand.Qform=adj$cand.Qform,
                  cand.gAdj=adj$cand.gAdj, cand.gform=adj$cand.gform,
                  V=SETTINGS$V,
                  break.match=T, one.sided=SETTINGS$one.sided, 
                  alt.smaller=SETTINGS$alt.smaller, 
                  verbose=F)[,these.cols]
    
  } else{
    
    # if missingness, then adjust using LTMLE 
    if(verbose) print('Adjusting for missingness')
    
    id <- dt$id
    dt <- dt[,c(SETTINGS$W,'A', 'Delta', 'Y')]
    
    out <- suppressMessages ( ltmle(data=dt, Anodes='A',  Cnodes='Delta',
                                    Ynodes='Y', abar=list(1,0), 
                                    SL.library=SETTINGS$SL.library, 
                                    SL.cvControl =list(V=SETTINGS$V), 
                                    estimate.time = F, variance.method='ic', id=id) )
    
    degree.free <- length(unique(id))-2
    est <- format.ltmle.output(out, goals=goals, degree.free=degree.free, 
                               SETTINGS=SETTINGS)  
  }
  
  
  colnames(est) <- c(paste0('INT', c('.est', '.CI.lo', '.CI.hi')),
                     paste0('CON', c('.est', '.CI.lo', '.CI.hi')),
                     paste0('EFFECT',c('.est', '.CI.lo', '.CI.hi', '.se', '.pval', '.QAdj', '.gAdj'))
  )
  
  cbind(Approach=this_Yc, N=nrow(dt), est)
}





format.ltmle.output <- function(out, goals, degree.free, SETTINGS){
  
  temp <- summary(out)$effect.measures
  INT <- temp$treatment
  CON <- temp$control
  ATE <- temp$ATE
  RR <- temp$RR
  
  if(goals=='aRR'){
    if(is.null(temp$RR)){
      ee <- data.frame( cbind(est=NA, CI.lo=NA, CI.hi=NA, se=NA, pval=NA))
    } else{
      ee <- get.inference(goal='aRR', psi.hat=log(RR$estimate), 
                          se=RR$std.dev, df=degree.free, 
                          one.sided=SETTINGS$one.sided, 
                          alt.smaller=SETTINGS$alt.smaller)
    }
    
  }else if (goals=='RD'){
    ee <- get.inference(goal='RD', psi.hat=ATE$estimate, 
                        se=ATE$std.dev, df=degree.free, 
                        one.sided=SETTINGS$one.sided, 
                        alt.smaller=SETTINGS$alt.smaller)
  } else{
    print('Error in transforming LTMLE results')
  }
  p <- data.frame(cbind( INT$estimate, INT$CI, 
                         CON$estimate, CON$CI,
                         ee[,c('est','CI.lo', 'CI.hi', 'se','pval')],
                         QAdj=NA, gAdj=NA
  ))
  colnames(p) <- c(paste0('INT', c('.est', '.CI.lo', '.CI.hi')),
                   paste0('CON', c('.est', '.CI.lo', '.CI.hi')),
                   paste0('EFFECT',c('.est', '.CI.lo', '.CI.hi', '.se', '.pval', '.QAdj', '.gAdj'))
  )
  p
}
