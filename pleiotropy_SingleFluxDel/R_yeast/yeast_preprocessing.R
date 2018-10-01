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
path<-"/Users/deya/Desktop/perpaer/Data_scripts_models/yeast7.6/commen/"
fluxPercent<-seq(0,1,0.005)
source(paste0(path,'apply_constraint.R'))# WT flux(pFBA)
source(paste0(path,'add_Ex_reaction.R'))# add Exchange reactions to the model
source(paste0(path,'onefluxdel.R'))

#------------------------------------------------------------------------------#
#                              loading Model                                  #
#------------------------------------------------------------------------------#
path1<-"/Users/deya/Desktop/perpaer/Data_scripts_models/yeast7.6/data/"
print(load(paste0(path1,"yeast_7_6.Rdata")))

######################
#biomass 5
lowbnd(yeast)[2984]<-0
uppbnd(yeast)[2984]<-1000
######################
#turn off biomass 8=
lowbnd(yeast)[3493]<-0
uppbnd(yeast)[3493]<-0



#------------------------------------------------------------------------------#
#                   WT flux with PFBA                                          #
# to constraint the flux[specific flux Percent] of its corropsponding reactions
# (for specific gene)
#------------------------------------------------------------------------------#


Emodel<-constraint_minmial(yeast)# I am only intersted on  WT fluxes

#------------------------------------------------------------------------------#
#           essential components for yeast model                                 #
#------------------------------------------------------------------------------#


Ebio<-biomass_reactsADDALL1(Emodel,2984)#

#Ebio[[1]]# model with  a newly added exchange reactions
ReacEx<-Ebio[[2]]#essential components

RNX_id<-match(ReacEx,react_id(Ebio[[1]]))



#------------------------------------------------------------------------------#
#                              important reactions                       #
#------------------------------------------------------------------------------#

exposWT<-react_pos(findExchReact(Ebio[[1]]))
RXNWT<-seq(1:react_num(Ebio[[1]]))
RXNWT1<-RXNWT[-c(exposWT)]













