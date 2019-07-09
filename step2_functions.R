
#-----------------------------------------------------------------------------------------------------------------------------
# baron & kenny estimator for natural direct and indirect effects 
#-----------------------------------------------------------------------------------------------------------------------------

baron_kenny <- function(obsdat, ymodel, mmodel) {
  
  # fit outcome model 
  fit1 <- glm(formula=ymodel, data=obsdat)
  # fit mediator model
  fit2 <- glm(formula=mmodel, data=obsdat)
  
  # bootstrap for variance of NIE  
  
  bk_boot_nie <- function(iteration) {
    boot_df <- obsdat[sample(nrow(obsdat), replace=T),]
    
    b_fit1 <- glm(formula=ymodel, data=boot_df)
    b_fit2 <- glm(formula=mmodel, data=boot_df) 
    
    b_nie <- coef(summary(b_fit1))[3]*coef(summary(b_fit2))[2]
    return(b_nie)
  }
  
  nie_b_list <- lapply(1:1000, function(x) bk_boot_nie(x))
  nie_b <- do.call(rbind, nie_b_list)
  
  # direct effect results
  bk_nde <- data.frame(nde = coef(summary(fit1))[2], 
                       nde_lb = coef(summary(fit1))[2] - 1.96*coef(summary(fit1))[2,2], 
                       nde_ub = coef(summary(fit1))[2] + 1.96*coef(summary(fit1))[2,2])
  
  # indirect effect results
  bk_nie <- data.frame(nie = coef(summary(fit1))[3]*coef(summary(fit2))[2], 
                       nie_lb = quantile(nie_b,0.025), nie_ub = quantile(nie_b,0.975))
  
  
  bk <- cbind(bk_nde, bk_nie) 
  return(bk) 
}

#-----------------------------------------------------------------------------------------------------------------------------
# inverse odds ratio weighting estimator for natural direct and indirect effects on the risk difference scale 
#-----------------------------------------------------------------------------------------------------------------------------

iorw <- function(obsdat, ammodel, ymmodel) {
  
  # calculate inverse odds weights
  am_fit <- glm(formula = ammodel, family="quasibinomial", data=obsdat)
  obsdat$p <- predict(am_fit, type="response")
  obsdat$iow <- (1 - obsdat$p)/obsdat$p 
  # assign unexposed weight of 1
  obsdat$iow[obsdat$a==0] <- 1
  
  # direct effect estimation 
  nde_fit <- glm(formula = ymmodel, family="gaussian", weights=iow, data=obsdat)
  te_fit <- glm(formula = ymmodel, family="gaussian", data=obsdat)
  
  # bootstrap for inference 
  iorw_boot <- function(iteration) {
    boot_df <-  obsdat[sample(nrow(obsdat), replace=T),]
    
    am_fit <- glm(formula=ammodel, family="quasibinomial", data=boot_df)
    p <- predict(am_fit, type="response")
    
    iow <- (1 - p)/p 
    iow[boot_df$a==0] <- 1
    
    nde_fit <- glm(formula=ymmodel, family="gaussian", weights=iow, data=boot_df)
    te_fit <- glm(formula=ymmodel, family="gaussian", data=boot_df)
    
    iorw_b <- c(NDE = coef(nde_fit)[2], NIE = coef(te_fit)[2] - coef(nde_fit)[2])
    return(iorw_b)
  }
  
  iorw_b_list <- lapply(1:250, function(x) iorw_boot(x))
  iorw_b_df <- data.frame(do.call(rbind, iorw_b_list))
  
  # direct effect results
  iorw_nde <- data.frame(nde = coef(nde_fit)[2], 
                         nde_lb = coef(nde_fit)[2] - 1.96*sqrt(var(iorw_b_df$NDE)), 
                         nde_ub = coef(nde_fit)[2] + 1.96*sqrt(var(iorw_b_df$NDE)))
  
  # indirect effect results
  iorw_nie <- data.frame(nie = coef(te_fit)[2] - coef(nde_fit)[2], 
                         nie_lb = (coef(te_fit)[2] - coef(nde_fit)[2]) - 1.96*sqrt(var(iorw_b_df$NIE)), 
                         nie_ub = (coef(te_fit)[2] - coef(nde_fit)[2]) + 1.96*sqrt(var(iorw_b_df$NIE)))
  
  
  iorw <- cbind(iorw_nde, iorw_nie)
  return(iorw)
}
#-----------------------------------------------------------------------------------------------------------------------------
# tmle
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
#     
# This is a function to estimate stochastic direct and indirect effects when there are direct effects of the 
#   exposure on the mediator and the outcome, and there is a mediator-outcome confounder affected by prior exposure.

# The stochastic direct effect is defined as E(Y_{1,g_0}) - E(Y_{0,g_0})
#   and the stochastic indirect effect is defined as E(Y_{1,g_1}) - E(Y_{1,g_0})
#   where Y_{a,g_a*} is the counterfactual value of Y when a=a and the mediator m is drawn from 
#   the distribution g under a*. 
#
# Inputs to the function: 
#     1. obsdat 
#       - should be a data frame that includes an exposure variable (called a with values 0/1), 
#              a mediator-outcome confounder affected by prior exposure (called z with values 0/1), 
#              a mediator variable (called m with values 0/1), and an outcome variable (called y with values 0/1). 
#       - there can be other covariates that can be named anything and have any values.
#     2. amodel 
#       - the parametric model for a that includes the parents of the exposure. For example: "a ~ w1 + w2 + w3"
#     3. zmodel 
#       - the parametric model for z that includes the parents of the mediator-outcome confounder 
#              affected by prior exposure. For example: "z ~ a + w1 + w2 + w3" 
#     4. mmodel 
#       - the parametric model for m that includes the parents of the mediator. 
#     5. ymodel 
#       - the parametric model for y that includes the parents of the outcome. 
#     6. qmodel 
#       - the covariates to marginalize over. For example, "w1 + w2 + w3"
#  
#      
#       there needs to be attached to obsdat the gma0 and gma1 values 
#-----------------------------------------------------------------------------------------------------------------------------
medtmle_intermedvar <- function(obsdat, amodel, zmodel, mmodel, ymodel, qmodel) {
  
  dfm1 <- dfm0 <- tmpdat <- obsdat
  dfm1$m <- 1 
  dfm0$m <- 0
  
  # calculate a weights
  psa1 <- I(tmpdat$a==1)/predict(glm(formula=amodel, family="binomial", data=tmpdat), type="response")
  psa0 <- I(tmpdat$a==0)/(1-predict(glm(formula=amodel, family="binomial", data=tmpdat), type="response"))
  
  # calculate m weights
  mazw <-predict(glm(formula=mmodel, family="binomial", data=tmpdat), type="response")
  psm<- I(tmpdat$m==1)*mazw + I(tmpdat$m==0)*(1-mazw)
  
  
  # first generate initial estimate of y in full pop, and then under m=0 and m=1 
  tmpdat$qyinit<-cbind(predict(glm(formula=ymodel, family="binomial", data=tmpdat), newdata=tmpdat, type="response"), 
                       predict(glm(formula=ymodel, family="binomial", data=tmpdat), newdata=dfm0, type="response"),
                       predict(glm(formula=ymodel, family="binomial", data=tmpdat), newdata=dfm1, type="response"))
  
  #-------------------------------------------------------------------------------------------------------------
  #   get E(Y_{1 , g_0})
  #-------------------------------------------------------------------------------------------------------------
  # calculate weights when a=1 and gm when a=0 
  tmpdat$ha1gma0 <- ((I(tmpdat$m==1)*tmpdat$gma0 + I(tmpdat$m==0)*(1-tmpdat$gma0))/psm) * psa1 
  
  # targeting step using weights and logit(q1) as an offset
  # epsilon is calculated using full data, then applied to estimate of y under m=0 and m=1
  epsilona1gma0 <- coef(glm(y ~  1 , weights=ha1gma0, offset=(qlogis(qyinit[,1])), family="quasibinomial", data=tmpdat))
  tmpdat$qyupm0a1gma0 <- plogis(qlogis(tmpdat$qyinit[,2]) + epsilona1gma0)
  tmpdat$qyupm1a1gma0 <- plogis(qlogis(tmpdat$qyinit[,3]) + epsilona1gma0)
  
  # marginalize over m 
  tmpdat$Qa1gma0 <- tmpdat$qyupm0a1gma0*(1-tmpdat$gma0) + tmpdat$qyupm1a1gma0*tmpdat$gma0
  
  # regress Q on covariates among those with a=1
  Qa1g0fit <- glm(formula = paste("Qa1gma0",qmodel, sep="~"), data=tmpdat[tmpdat$a==1,], family="quasibinomial")
  Qa1g0 <- predict(Qa1g0fit, type="response", newdata=tmpdat)
  
  # update again (necessary when a is nonrandom)
  epsilonza1gma0 <-coef(glm(Qa1gma0~ 1 , weights=psa1, offset=qlogis(Qa1g0), family="quasibinomial", data=tmpdat))
  Qzupa1gma0 <-plogis(qlogis(Qa1g0) + epsilonza1gma0)
  
  # calculate tmle for y under a=1 and gm where a=0
  tmlea1m0 <- mean(Qzupa1gma0)
  
  #-------------------------------------------------------------------------------------------------------------
  #   get E(Y_{0 , g_0})
  #-------------------------------------------------------------------------------------------------------------
  # calculate weights when a=0 and gm when a=0 
  tmpdat$ha0gma0 <- ((tmpdat$m*tmpdat$gma0 + (1-tmpdat$m)*(1-tmpdat$gma0))/psm) * psa0 
  
  # targeting step 
  epsilona0gma0 <- coef(glm(y ~  1 , weights=ha0gma0, offset=(qlogis(qyinit[,1])), family="quasibinomial", data=tmpdat))
  tmpdat$qyupm0a0gma0 <- plogis(qlogis(tmpdat$qyinit[,2]) + epsilona0gma0)
  tmpdat$qyupm1a0gma0 <- plogis(qlogis(tmpdat$qyinit[,3]) + epsilona0gma0)
  
  # marginalize over m 
  tmpdat$Qa0gma0 <- tmpdat$qyupm0a0gma0*(1-tmpdat$gma0) + tmpdat$qyupm1a0gma0*tmpdat$gma0
  
  # regress Q on covariates among those with a=0
  Qa0g0fit <- glm(formula = paste("Qa0gma0",qmodel, sep="~"), data=tmpdat[tmpdat$a==0,], family="quasibinomial")
  Qa0g0 <- predict(Qa0g0fit, type="response", newdata=tmpdat)
  
  # update again (necessary when a is nonrandom)
  epsilonza0gma0 <- coef(glm(Qa0gma0~ 1 , weights=psa0, offset=qlogis(Qa0g0), family="quasibinomial", data=tmpdat))
  Qzupa0gma0 <-plogis(qlogis(Qa0g0) + epsilonza0gma0)
  
  # calculate tmle for y under a=0 and gm where a=0
  tmlea0m0 <- mean(Qzupa0gma0)
  
  #-------------------------------------------------------------------------------------------------------------
  #   get E(Y_{1 , g_1})
  #-------------------------------------------------------------------------------------------------------------
  # calculate weights when a=1 and gm when a=1 
  tmpdat$ha1gma1<- ((tmpdat$m*tmpdat$gma1 + (1-tmpdat$m)*(1-tmpdat$gma1))/psm) * psa1 
  
  # targeting step
  epsilona1gma1<-coef(glm(y ~  1 , weights=ha1gma1, offset=(qlogis(qyinit[,1])), family="quasibinomial", data=tmpdat))
  tmpdat$qyupm0a1gma1<-plogis(qlogis(tmpdat$qyinit[,2]) + epsilona1gma1)
  tmpdat$qyupm1a1gma1<-plogis(qlogis(tmpdat$qyinit[,3]) + epsilona1gma1)
  
  # marginalize over m
  tmpdat$Qa1gma1<-tmpdat$qyupm0a1gma1*(1-tmpdat$gma1) + tmpdat$qyupm1a1gma1*tmpdat$gma1
  
  # regress Q on covariates among those with a=1
  Qa1g1fit <- glm(formula = paste("Qa1gma1",qmodel, sep="~"), data=tmpdat[tmpdat$a==1,], family="quasibinomial")
  Qa1g1 <- predict(Qa1g1fit, type="response", newdata=tmpdat)
  
  # update again (necessary when a is nonrandom)
  epsilonza1gma1 <-coef(glm(Qa1gma1~ 1 , weights=psa1, offset=qlogis(Qa1g1), family="quasibinomial", data=tmpdat))
  Qzupa1gma1 <- plogis(qlogis(Qa1g1) + epsilonza1gma1)
  
  # calculate tmle for y under a=1 and gm where a=1
  tmlea1m1 <- mean(Qzupa1gma1)
  
  #-------------------------------------------------------------------------------------------------------------
  #   estimate variances using EIC 
  #-------------------------------------------------------------------------------------------------------------
  
  tmpdat$qyupa1g0<-plogis(qlogis(tmpdat$qyinit[,1]) + epsilona1gma0)
  tmpdat$qyupa1g1<-plogis(qlogis(tmpdat$qyinit[,1]) + epsilona1gma1)
  tmpdat$qyupa0g0<-plogis(qlogis(tmpdat$qyinit[,1]) + epsilona0gma0)
  
  # EIC for E(Y_{1, g_0})
  eic1a1g0<-tmpdat$ha1gma0 * (tmpdat$y - tmpdat$qyupa1g0)
  eic2a1g0<-psa1*(tmpdat$Qa1gma0 - Qzupa1gma0)
  eic3a1g0<-Qzupa1gma0 -tmlea1m0
  
  eica1g0<-eic1a1g0 + eic2a1g0 + eic3a1g0
  
  # EIC for E(Y_{1, g_1})
  eic1a1g1<-tmpdat$ha1gma1 * (tmpdat$y - tmpdat$qyupa1g1)
  eic2a1g1<-psa1*(tmpdat$Qa1gma1- Qzupa1gma1)
  eic3a1g1<-Qzupa1gma1-tmlea1m1
  
  eica1g1<-eic1a1g1 + eic2a1g1 + eic3a1g1
  
  # EIC for E(Y_{0, g_0})
  eic1a0g0<-tmpdat$ha0gma0 * (tmpdat$y - tmpdat$qyupa0g0)
  eic2a0g0<-psa0*(tmpdat$Qa0gma0 - Qzupa0gma0)
  eic3a0g0<-Qzupa0gma0 -tmlea0m0
  
  eica0g0<-eic1a0g0 + eic2a0g0 + eic3a0g0
  
  
  #-------------------------------------------------------------------------------------------------------------
  #   combine to get mediation parameters, variances, and 95% confidence intervals  
  #-------------------------------------------------------------------------------------------------------------
  
  sde_tmle <- tmlea1m0-tmlea0m0
  sde_eic <- eica1g0 - eica0g0
  var_sde_eic <- var(sde_eic)/nrow(tmpdat)
  
  sie_tmle <- tmlea1m1-tmlea1m0
  sie_eic <- eica1g1 - eica1g0
  var_sie_eic <- var(sie_eic)/nrow(tmpdat)
  
  results <- data.frame(cbind(sde=sde_tmle, sde_var=var_sde_eic, sie=sie_tmle, sie_var=var_sie_eic, 
                              sde_lb = sde_tmle - 1.96*sqrt(var_sde_eic), sde_ub = sde_tmle + 1.96*sqrt(var_sde_eic), 
                              sie_lb = sie_tmle - 1.96*sqrt(var_sie_eic), sie_ub = sie_tmle + 1.96*sqrt(var_sie_eic)))
  
  return(results)
}

#-----------------------------------------------------------------------------------------------------------------------------
#     
# This is a function to estimate stochastic direct and indirect effects when there are direct effects of the 
#   exposure on the mediator and the outcome, and there is no mediator-outcome confounder affected by prior exposure.

# The stochastic direct effect is defined as E(Y_{1,g_0}) - E(Y_{0,g_0})
#   and the stochastic indirect effect is defined as E(Y_{1,g_1}) - E(Y_{1,g_0})
#   where Y_{a,g_a*} is the counterfactual value of Y when a=a and the mediator m is drawn from 
#   the distribution g under a*. 
#
# Inputs to the function: 
#     1. obsdat 
#       - should be a data frame that includes an exposure variable (called a with values 0/1), 
#              a mediator variable (called m with values 0/1), and an outcome variable (called y with values 0/1). 
#       - there can be other covariates that can be named anything and have any values.
#     2. amodel 
#       - the parametric model for a that includes the parents of the exposure. For example: "a ~ w1 + w2 + w3"
#     3. mmodel 
#       - the parametric model for m that includes the parents of the mediator. 
#     4. ymodel 
#       - the parametric model for y that includes the parents of the outcome. 
#     5. qmodel 
#       - the covariates to marginalize over. For example, "w1 + w2 + w3"
#  
#      
#       
#-----------------------------------------------------------------------------------------------------------------------------
medtmle_nointermedvar <- function(obsdat, amodel, mmodel, ymodel, qmodel) {
  
  dfa1 <- dfa0 <- dfm1 <- dfm0 <- tmpdat <- obsdat
  
  dfa1$a <- dfm1$m <-  1 
  dfa0$a <- dfm0$m <-  0 
  
  # calculate a weights 
  psa1 <- I(tmpdat$a==1)/predict(glm(formula=amodel, family="binomial", data=tmpdat), type="response")
  psa0 <- I(tmpdat$a==0)/(1-predict(glm(formula=amodel, family="binomial", data=tmpdat), type="response"))
  
  # calculate m weights 
  ma <- predict(glm(formula=mmodel, family="binomial", data=tmpdat), type="response")
  psm <- (ma*tmpdat$m) + ((1-ma)*(1-tmpdat$m))
  
  
  # first generate initial estimate of y in full pop, and then under m=0 and m=1 
  tmpdat$qyinit <- cbind(predict(glm(formula=ymodel, family="binomial", data=tmpdat), newdata=tmpdat, type="response"), 
                         predict(glm(formula=ymodel, family="binomial", data=tmpdat), newdata=dfm0, type="response"),
                         predict(glm(formula=ymodel, family="binomial", data=tmpdat), newdata=dfm1, type="response"))
  
  #-------------------------------------------------------------------------------------------------------------
  #   get E(Y_{1 , g_0})
  #-------------------------------------------------------------------------------------------------------------
  # calculate weights when a=1 and gm when a=0 
  tmpdat$ha1gma0 <- ((tmpdat$m*tmpdat$gma0 + (1-tmpdat$m)*(1-tmpdat$gma0))/psm) * psa1 
  
  # targeting step using weights and logit(q1) as an offset
  # epsilon is calculated using full data, then applied to estimate of y under m=0 and m=1
  epsilona1gma0 <-coef(glm(y ~  1 , weights=ha1gma0 , offset=(qlogis(qyinit[,1])), family="quasibinomial", data=tmpdat))
  tmpdat$qyupm0a1gma0 <- plogis(qlogis(tmpdat$qyinit[,2]) + epsilona1gma0)
  tmpdat$qyupm1a1gma0 <- plogis(qlogis(tmpdat$qyinit[,3]) + epsilona1gma0)
  
  # marginalize over m 
  tmpdat$Qa1gma0 <-tmpdat$qyupm0a1gma0*(1-tmpdat$gma0) + tmpdat$qyupm1a1gma0*tmpdat$gma0
  
  # regress Q on W among those with A=1 and predict
  Qzfita1gma0 <- glm(formula=paste("Qa1gma0", qmodel, sep="~"), data=tmpdat[tmpdat$a==1,], family="quasibinomial")
  Qza1gma0 <-predict(Qzfita1gma0, type="response", newdata=tmpdat)
  
  # update again (necessary when a is nonrandom)
  epsilonza1gma0 <- coef(glm(Qa1gma0 ~ 1 , weights=psa1, offset=qlogis(Qza1gma0), family="quasibinomial", data=tmpdat))
  Qzupa1gma0 <- plogis(qlogis(Qza1gma0) + epsilonza1gma0)
  
  # calculate tmle for y under a=1 and gm where a=0
  tmlea1m0 <- mean(Qzupa1gma0)
  
  #-------------------------------------------------------------------------------------------------------------
  #   get E(Y_{0 , g_0})
  #-------------------------------------------------------------------------------------------------------------
  # calculate weights when a=0 and gm when a=0 
  tmpdat$ha0gma0 <- ((tmpdat$m*tmpdat$gma0 + (1-tmpdat$m)*(1-tmpdat$gma0))/psm) * psa0 
  
  # targeting step 
  epsilona0gma0 <-coef(glm(y ~  1 , weights=ha0gma0, offset=(qlogis(qyinit[,1])), family="quasibinomial", data=tmpdat))
  tmpdat$qyupm0a0gma0 <-plogis(qlogis(tmpdat$qyinit[,2]) + epsilona0gma0)
  tmpdat$qyupm1a0gma0 <-plogis(qlogis(tmpdat$qyinit[,3]) + epsilona0gma0)
  
  # marginalize over m 
  tmpdat$Qa0gma0 <-tmpdat$qyupm0a0gma0*(1-tmpdat$gma0) + tmpdat$qyupm1a0gma0*tmpdat$gma0
  
  # regress Q on W among those with A=0 and predict
  Qzfita0gma0 <- glm(formula=paste("Qa0gma0", qmodel, sep="~"), data=tmpdat[tmpdat$a==0,], family="quasibinomial")
  Qza0gma0 <- predict(Qzfita0gma0, type="response", newdata=tmpdat)
  
  # update again (necessary when a is nonrandom)
  epsilonza0gma0 <- coef(glm(Qa0gma0 ~ 1 , weights=psa0, offset=qlogis(Qza0gma0), family="quasibinomial", data=tmpdat))
  Qzupa0gma0 <- plogis(qlogis(Qza0gma0) + epsilonza0gma0)
  
  # calculate tmle for y under a=0 and gm where a=0
  tmlea0m0 <- mean(Qzupa0gma0)
  
  #-------------------------------------------------------------------------------------------------------------
  #   get E(Y_{1 , g_1})
  #-------------------------------------------------------------------------------------------------------------
  # calculate weights when a=1 and gm when a=1 
  tmpdat$ha1gma1 <- ((tmpdat$m*tmpdat$gma1 + (1-tmpdat$m)*(1-tmpdat$gma1))/psm) * psa1 
  
  # targeting step
  epsilonma1gma1<-coef(glm(y ~  1 , weights=ha1gma1, offset=(qlogis(qyinit[,1])), family="quasibinomial", data=tmpdat))
  tmpdat$qyupm0a1gma1<-plogis(qlogis(tmpdat$qyinit[,2]) + epsilonma1gma1)
  tmpdat$qyupm1a1gma1<-plogis(qlogis(tmpdat$qyinit[,3]) + epsilonma1gma1)
  
  # marginalize over m
  tmpdat$Qa1gma1<-tmpdat$qyupm0a1gma1*(1-tmpdat$gma1) + tmpdat$qyupm1a1gma1*tmpdat$gma1
  
  # regress Q on W among those with A=1 and predict
  Qzfita1gma1<-glm(formula=paste("Qa1gma1", qmodel, sep="~"), data=tmpdat[tmpdat$a==1,], family="quasibinomial")
  Qza1gma1<-predict(Qzfita1gma1, type="response", newdata=tmpdat)
  
  # update again (necessary when a is nonrandom)
  epsilonza1gma1<-coef(glm(Qa1gma1~ 1 , weights=psa1, offset=qlogis(Qza1gma1), family="quasibinomial", data=tmpdat))
  Qzupa1gma1<-plogis(qlogis(Qza1gma1) + epsilonza1gma1)
  
  # calculate tmle for y under a=1 and gm where a=1
  tmlea1m1 <- mean(Qzupa1gma1)
  
  #-------------------------------------------------------------------------------------------------------------
  #   estimate variances using EIC  
  #-------------------------------------------------------------------------------------------------------------
  
  tmpdat$qyupa1g0<-plogis(qlogis(tmpdat$qyinit[,1]) + epsilona1gma0)
  tmpdat$qyupa1g1<-plogis(qlogis(tmpdat$qyinit[,1]) + epsilonma1gma1)
  tmpdat$qyupa0g0<-plogis(qlogis(tmpdat$qyinit[,1]) + epsilona0gma0)
  
  # EIC for E(Y_{1, g_0})
  eic1a1g0 <- tmpdat$ha1gma0 * (tmpdat$y - tmpdat$qyupa1g0)
  eic2a1g0 <- psa1*(tmpdat$Qa1gma0 - Qzupa1gma0)
  eic3a1g0 <- Qzupa1gma0 -tmlea1m0
  
  eica1g0 <- eic1a1g0 + eic2a1g0 + eic3a1g0
  
  # EIC for E(Y_{1, g_1})
  eic1a1g1 <- tmpdat$ha1gma1 * (tmpdat$y - tmpdat$qyupa1g1)
  eic2a1g1 <- psa1*(tmpdat$Qa1gma1- Qzupa1gma1)
  eic3a1g1 <- Qzupa1gma1-tmlea1m1
  
  eica1g1 <- eic1a1g1 + eic2a1g1 + eic3a1g1
  
  # EIC for E(Y_{0, g_0})
  eic1a0g0 <- tmpdat$ha0gma0 * (tmpdat$y - tmpdat$qyupa0g0)
  eic2a0g0 <- psa0*(tmpdat$Qa0gma0 - Qzupa0gma0)
  eic3a0g0 <- Qzupa0gma0 - tmlea0m0
  
  eica0g0 <- eic1a0g0 + eic2a0g0 + eic3a0g0
  
  #-------------------------------------------------------------------------------------------------------------
  #   combine to get mediation parameters, variances, and 95% confidence intervals
  #-------------------------------------------------------------------------------------------------------------
  
  sde_tmle <- tmlea1m0-tmlea0m0
  sde_eic <- eica1g0 - eica0g0
  var_sde_eic <- var(sde_eic)/nrow(tmpdat)
  
  sie_tmle <- tmlea1m1-tmlea1m0
  sie_eic <- eica1g1 - eica1g0
  var_sie_eic <- var(sie_eic)/nrow(tmpdat)
  
  results <- data.frame(cbind(sde=sde_tmle, sde_var=var_sde_eic, sie=sie_tmle, sie_var=var_sie_eic, 
                              sde_lb = sde_tmle - 1.96*sqrt(var_sde_eic), sde_ub = sde_tmle + 1.96*sqrt(var_sde_eic), 
                              sie_lb = sie_tmle - 1.96*sqrt(var_sie_eic), sie_ub = sie_tmle + 1.96*sqrt(var_sie_eic)))
  
  return(results)
}

#-----------------------------------------------------------------------------------------------------------------------------
#     
#  power function using equation from Vittinghof et al 
#
#-----------------------------------------------------------------------------------------------------------------------------

mpower_eq <- function(dat, mmodel, ymodel, MonY) {
  
  # define vars for power calculation
  # calculate multiple correlation 
  corfit <- glm(formula=mmodel, data=dat)
  mhat <- predict(corfit)
  rho <- cov(mhat, dat$m)/(sqrt(var(mhat))*sqrt(var(dat$m)))
  
  # set type 1 error = 0.975
  z_a <- qnorm(0.975) 
  # set type 2 error = 0.2, which implies power=0.8
  z_g <- qnorm(0.8)
  
  lmfit <- glm(formula=ymodel, data=dat)
  b2 <- MonY 
  
  sigma_1 <- var(dat$a)
  sigma_2 <- var(dat$m)
  
  sigma_e <- var(residuals(lmfit))
  
  req_ss <- ((z_a + z_g)^2*sigma_e)/((b2^2*sigma_2)*(1-rho^2))
  
  power <- pnorm(sqrt((n*b2^2*sigma_2*(1-rho^2))/sigma_e) - z_a)
  
  results <-  cbind(lb = NA, ub=NA, ed = power, coverage = 1.05, parameter="NIE")
  return(results)
}