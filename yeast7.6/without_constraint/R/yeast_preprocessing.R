library(sybil)
library(sybilSBML)
library(lattice)
library(cplexAPI)
library(parallel)
library(methods)
library(ggplot2)


SYBIL_SETTINGS("SOLVER", "cplexAPI")
SYBIL_SETTINGS("TOLERANCE",1E-06)
tol <- SYBIL_SETTINGS("TOLERANCE")

#------------------------------------------------------------------------------#
#                              sources                                         #
#------------------------------------------------------------------------------#
path<-"Data_scripts_models/yeast7.6/common/"
fluxPercent<-seq(0,1,0.005)# flux reduction
source(paste0(path,'apply_constraint.R'))# WT flux(pFBA)
source(paste0(path,'add_Ex_reaction.R'))# add Exchange reactions to the model
source(paste0(path,'onegenedel_030816_FBA.R'))# new onegenedel function without constraint[with new flux Percent parameter].
source(paste0(path,'sysBiolAlg_fba1.R'))
source(paste0(path,'optimizer1.R'))

#------------------------------------------------------------------------------#
#                              loading Model                                  #
#------------------------------------------------------------------------------#
path1<-"Data_scripts_models/yeast7.6/data/"
print(load(paste0(path1,"yeast_7_6.Rdata")))


#------------------------------------------------------------------------------#
#                   WT flux with PFBA                                          #
# to constraint the flux[specific flux Percent] of its corropsponding reactions
# (for specific gene)
#------------------------------------------------------------------------------#


mod2=mod2irrev(yeast)

Emodel<-constraint_minmial(mod2)

#------------------------------------------------------------------------------#
#           essential components for yeast model                                 #
#------------------------------------------------------------------------------#

Ebio<-biomass_reactsADDALL1(mod2,match(react_id(yeast)[2984],react_id(mod2)))
#Ebio[[1]]# model with  a newly added exchange reactions
ReacEx<-Ebio[[2]]#essential components [exchange reactions]

RNX_id<-match(ReacEx,react_id(Ebio[[1]]))


#------------------------------------------------------------------------------#
#           send orginal bounds for the model to the optimizer                 #
#------------------------------------------------------------------------------#

fdWT <- .generateFluxdels(mod2,seq(1:length(allGenes(mod2))),1)

fdWT[["mdlb"]][sapply(fdWT[["mdlb"]], is.null)] <- 0
fdWT[["mdub"]][sapply(fdWT[["mdub"]], is.null)] <- 0












