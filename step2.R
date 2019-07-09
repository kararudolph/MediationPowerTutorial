source("step2_functions.R")

# baron & kenny estimator for natural direct and indirect effects 
bk_results <- baron_kenny(obsdat, ymodel, mmodel)

# inverse odds ratio weighting estimator for natural direct and indirect effects on the risk difference scale 
iorw_results <- iorw(obsdat, ammodel, ymmodel) 

#  power function using equation from Vittinghof et al 
analysticEqn_results <- mpower_eq(fulldat, mmodel, ymodel, MonY) #tktk fulldat or obsdat; NA results?

#-----------------------------------------------------------------------------------------------------------------------------
#
# tmle for  DGMs 1-4:
#
#-----------------------------------------------------------------------------------------------------------------------------
tmle_results <- medtmle_nointermedvar(obsdat, amodel, mmodel, ymodel, qmodel) 

#-----------------------------------------------------------------------------------------------------------------------------
#
# tmle for  DGMs 5-12:
#
#-----------------------------------------------------------------------------------------------------------------------------
tmle_results <- medtmle_intermedvar(obsdat, amodel, zmodel, mmodel, ymodel, qmodel) 


  