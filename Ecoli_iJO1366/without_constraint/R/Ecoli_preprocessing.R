library(sybil)
library(sybilSBML)
library(lattice)
library(cplexAPI)
library(parallel)
library(methods)



SYBIL_SETTINGS("SOLVER", "cplexAPI")
SYBIL_SETTINGS("TOLERANCE",1E-06)
tol <- SYBIL_SETTINGS("TOLERANCE")

#------------------------------------------------------------------------------#
#                              sources                                         #
#------------------------------------------------------------------------------#
path<-"Data_scripts_models/Ecoli_iJO1366/common/"
fluxPercent<-seq(0,1,0.005)#flux reduction
source(paste0(path,'apply_constraint.R'))# WT flux(pFBA)
source(paste0(path,'add_Ex_reaction.R'))# add Exchange reactions to the model
source(paste0(path,'onegenedel_030816_FBA.R'))# new function
source(paste0(path,'sysBiolAlg_fba1.R'))
source(paste0(path,'optimizer1.R'))

#------------------------------------------------------------------------------#
#                              loading Model                                  #
#------------------------------------------------------------------------------#
path1<-"Data_scripts_models/Ecoli_iJO1366/data/"
print(load(paste0(path1,"iJO1366.Rdata")))

mod2=mod2irrev(Ecoli)
#------------------------------------------------------------------------------#
#                   WT flux with PFBA                                          #
# to constraint the flux[specific flux Percent] of its corropsponding reactions
# (for specific gene)
#------------------------------------------------------------------------------#


Emodel<-constraint_minmial(mod2)



#------------------------------------------------------------------------------#
#           essential components for yeast model                                 #
#------------------------------------------------------------------------------#

print(load(paste0(path1,"essential_comps_without_ATP.Rdata")))


Ebio<-biomass_reactsADDALL(mod2)
ReacEx<-sub( '(?<=.{0})', 'EX_' ,met_id(mod2)[essential_comps], perl=TRUE )
ReacEx<-c(ReacEx,react_id(mod2)[match("ATPM",react_id(mod2))])
RNX_id<-match(ReacEx,react_id(Ebio[[1]]))


#------------------------------------------------------------------------------#
#           send orginal bounds for the model for optimizer                    #
#------------------------------------------------------------------------------#

fdWT <- .generateFluxdels(mod2,seq(1:length(allGenes(mod2))),1)

fdWT[["mdlb"]][sapply(fdWT[["mdlb"]], is.null)] <- 0
fdWT[["mdub"]][sapply(fdWT[["mdub"]], is.null)] <- 0










