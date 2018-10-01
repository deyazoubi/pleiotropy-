# added a new exchange reaction to the model 
biomass_reactsADDALL1<-function(model,biomassR)
{
    if(biomassR==T){
        nam<-which(obj_coef(model)!=0)# biomass reaction number
    }else{
        nam<-  biomassR
    }
    Ex_biomass <- shrinkMatrix(model, j = nam, tol = 10^-10)
    
    Ex_biomass<-subset.matrix(Ex_biomass, Ex_biomass[,1] <tol)#influx biomass components (reactants)
    
    rns<-rownames(Ex_biomass)
    
    ReacEx<-sub( '(?<=.{0})', 'EX_' , rns, perl=TRUE )# new react id
    
    
    for(i in 1:length(ReacEx)){
        
        
        model <-addReact(model,id=ReacEx[i],met=rns[i],Scoef=-1,reversible =T,lb=0,ub=1000)
        print(optimizeProb(model))
        print(ReacEx[i])
    }
    
    model<-addReact(model,id="Max_ATP_ADP",met=met_id(model)[c(324,614,605,288,1045)],Scoef=c(-1,-1,1,1,1),reversible =F,lb=0,ub=1000)
    
    ReacEx<- ReacEx[-match("EX_0434",ReacEx)]#ATP
    ReacEx<- ReacEx[ -match("EX_0803",ReacEx)]#H2O
    ReacEx<- ReacEx[ -match("EX_1467",ReacEx)]#sulphate
   
    ReacEx<-c(ReacEx,"Max_ATP_ADP")
    
    uppbnd(model)[match("EX_0434",react_id(model))]<-0 # ex atp 3496
    #uppbnd(model)[3504]<-0 # ex h2o
    
    return (list( model,ReacEx))
}







