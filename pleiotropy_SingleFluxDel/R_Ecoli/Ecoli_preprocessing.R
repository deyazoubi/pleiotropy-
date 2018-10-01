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
path<-"/Users/deya/Desktop/perpaer/Data_scripts_models/Ecoli_iJo/commen/"
fluxPercent<-seq(0,1,0.005)
source(paste0(path,'apply_constraint.R'))# WT flux(pFBA)
source(paste0(path,'add_Ex_reaction.R'))# add Exchange reactions to the model
source(paste0(path,'onefluxdel.R'))

#------------------------------------------------------------------------------#
#                              loading Model                                  #
#------------------------------------------------------------------------------#
path1<-"/Users/deya/Desktop/perpaer/Data_scripts_models/Ecoli_iJo/data/"
print(load(paste0(path1,"iJO1366.Rdata")))


#------------------------------------------------------------------------------#
#                   WT flux with PFBA                                          #
# to constraint the flux[specific flux Percent] of its corropsponding reactions
# (for specific gene)
#------------------------------------------------------------------------------#



A<-748
Emodel<-constraint_minmial(Ecoli)
# new update #1 on 23.06.2015 case 1
lowbnd(Emodel)[A]<-3.15
uppbnd(Emodel)[A]<-1000


#------------------------------------------------------------------------------#
#           essential components for yeast model                                 #
#------------------------------------------------------------------------------#

print(load(paste0(path1,"essential_comps_without_ATP.Rdata")))

Ebio<-biomass_reactsADDALL(Emodel)
ReacEx<-sub( '(?<=.{0})', 'EX_' ,met_id(Emodel)[essential_comps], perl=TRUE )

ReacEx<-c(ReacEx,react_id(Emodel)[A])#Add ATPM instead of max EX_ATP  max ATPM

RNX_id<-match(ReacEx,react_id(Ebio[[1]]))


#------------------------------------------------------------------------------#
#                              important reactions                       #
#------------------------------------------------------------------------------#

exposWT<-react_pos(findExchReact(Ebio[[1]]))
RXNWT<-seq(1:react_num(Ebio[[1]]))
RXNWT1<-RXNWT[-c(exposWT)]












