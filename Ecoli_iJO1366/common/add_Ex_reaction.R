#Ecoli
biomass_reactsADDALL<-function(model)
{
    
    nam<-which(obj_coef(model)!=0)# biomass reaction number
    
    Ex_biomass <- shrinkMatrix(model, j = nam, tol = 10^-10)
    
    Ex_biomass<-subset.matrix(Ex_biomass, Ex_biomass[,1] <tol)#influx biomass components (reactants)
    
    rns<-rownames(Ex_biomass)
    ReacEx<-sub( '(?<=.{0})', 'EX_' , rns, perl=TRUE )# new react id
    
    
    for(i in 1:length(ReacEx)){
        
        
        model <-addReact(model,id=ReacEx[i],met=rns[i],Scoef=-1,reversible =T,lb=0,ub=1000)
        print(optimizeProb(model))
        print(ReacEx[i])
    }
    
    uppbnd(model)[match("EX_atp[c]",react_id(model))]<-0 #
    return (list( model,ReacEx))
}

