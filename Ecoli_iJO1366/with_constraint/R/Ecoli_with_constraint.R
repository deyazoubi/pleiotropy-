
#------------------------------------------------------------------------------#
#                              calculating pleiotropy                          #
#------------------------------------------------------------------------------#

# The following function simulates the effect of all possible gene deletions
# in the production of particular biomass components

affectedBiomassComponent <- function (comp,model,j,imp) {
    
    mod<-changeObjFunc(model, react=comp,obj_coef=1)
    wt1<-optimizeProb(mod)
    wt<-lp_obj(invisable_values(wt1 ))
    print(wt)
    sgds <- onegenedel(mod,alg="fba1",geneList=allGenes(mod)[imp],flux_Percent=j)
    ogd<-lp_obj(invisable_values(sgds ))
    red <- (abs(wt)-abs(ogd))/abs(wt) >= 1E-4
    gc()
    return ( red )
}

run_pl<-function(comp,model,fluxpercent,imp){
    
    pl <- mcmapply ( affectedBiomassComponent,comp, MoreArgs =list(model,fluxpercent,imp ), SIMPLIFY = TRUE ,mc.cores=2)
    gc()
    return(pl)
}






#------------------------------------------------------------------------------#
#simulated mutations of different severity by restricting the flux through
#all reactions catalyzed by the gene to a fixed fraction of the wildtype flux,
#which I reduced in 0.5%-steps
#------------------------------------------------------------------------------#

fluxPercent<-seq(0,1,0.005)

for(j in 1:length(fluxPercent)){
    
    #   If Epl [i, j] == true, then the deletion of gene i has reduced production of the metabolite j
    Epls<-run_pl(RNX_id,Ebio[[1]],fluxPercent[j],seq(1:length(allGenes(Emodel))))
    
    
    
    # Count for each gene, for how many biomass components,
    # the maximum production rate was reduced by the deletion of the gene.
    EplSum<-rowSums(Epls)
    
    
    fname = format(Sys.time(), paste0(paste0(j,"_at_"),paste0("WT_matrix_Ecoli","%Y%m%d_%H%M.Rdata")))
    save(Epl,file=paste0("results/",fname))
    save.image()
    
    
    
    fname = format(Sys.time(), paste0(paste0(j,"_at_"),paste0("WT_SUM_Perc_Ecoli","%Y%m%d_%H%M.Rdata")))
    save(EplSum,file=paste0("results/",fname))
    save.image()
    
    
}













