# run power simulations

source("step3_functions.R")

#-----------------------------------------------------------------------------------------------------------------------------
#
# for  DGMs 1-4:
#
#-----------------------------------------------------------------------------------------------------------------------------

all_results_list_a <- lapply(1:simN, function(x) mpower_nointermedvar(x,superN, n, fulldat=fulldat, truth=truth, 
                                                                      amodel=amodel, mmodel=mmodel, 
                                                                      ymodel=ymodel, qmodel=qmodel, 
                                                                      ammodel=ammodel, ymmodel=ymmodel))
# combine results
all_results_t <- do.call(rbind, all_results_list_a)
# make sure these variables are numeric
all_results_t$ed <- as.numeric(levels(all_results_t$ed))[all_results_t$ed]
all_results_t$coverage <- as.numeric(levels(all_results_t$coverage))[all_results_t$coverage]

# summarize power and coverage 
power_results <- all_results_t %>% group_by(parameter, method) %>% summarise(power = mean(ed), coverage=mean(coverage))
power_results$effect_size <- ifelse(power_results$parameter=="NDE",truth$nde,truth$nie)


#-----------------------------------------------------------------------------------------------------------------------------
#
# for  DGMs 5-12:
#
#-----------------------------------------------------------------------------------------------------------------------------
all_results_list_a <- lapply(1:simN, function(x) mpower_intermedvar(x,superN, n, fulldat=fulldat, truth=truth,
                                                                    amodel=amodel, mmodel=mmodel,
                                                                    ymodel=ymodel, qmodel=qmodel, ammodel=ammodel,
                                                                    ymmodel=ymmodel,zmodel=zmodel,
                                                                    ammodel_noz=ammodel_noz, mmodel_noz=mmodel_noz,
                                                                    ymodel_noz=ymodel_noz, ymmodel_noz=ymmodel_noz))
# combine results
all_results_t <- do.call(rbind, all_results_list_a)
# make sure these variables are numeric
all_results_t$ed <- as.numeric(levels(all_results_t$ed))[all_results_t$ed]
all_results_t$coverage <- as.numeric(levels(all_results_t$coverage))[all_results_t$coverage]
# summarize power and coverage
power_results <- all_results_t %>% group_by(parameter, method) %>% summarise(power = mean(ed), coverage=mean(coverage))
power_results$effect_size <- ifelse(power_results$parameter=="NDE",truth$sde,truth$sie)

