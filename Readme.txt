This file gives a brief description of the contents of the Data_scripts_models folder.

Ecoli_iJO1366 - functions for modeling E. coli
	common - R functions called in the individual simulations
	data - Rdata files with the iJO1366 model and essential metabolites
	with_constraint - R files used for the simulations that constrain fluxes by WT fluxes 
		-results
			PLWT_Ecoli:Rdata file with a wildtype result at different wildtype flux ratio 
			PL_NADPH_Ecoli:Rdata file with the result when NADPH is made freely available at different wildtype flux ratio
			Free_met_at_full-knockout: Rdata files with the results when we allowed several cofactors freely available at full knockout 
			
	without_constraint - R files used for the simulations that do NOT constrain fluxes by WT fluxes
	 	-results
			PLWT_Ecoli_without_constraint:Rdata file with a wildtype result at different wildtype flux ratio[ without applying the constraint] 
			PL_NADPH_Ecoli_without_constraint:Rdata file with the result when NADPH is made freely available at different wildtype 
							  flux ratio[ without applying the constraint] 
			Free_met_at_full-knockout: Rdata files with the results when we allowed several cofactors freely available at full knockout 
						  [without applying the constraint] 
	

yeast7.6 -	functions for modeling yeast
	common - R functions called in the individual simulations
	data - Rdata files with the yeast7.6 model and essential metabolites
	with_constraint - R files used for the simulations that constrain fluxes by WT fluxes 
		-results
			PLWT_yeast:Rdata file with a wildtype result at different wildtype flux ratio 
			PL_NADPH_yeast:Rdata file with the result when NADPH is made freely available at different wildtype flux ratio
			Free_met_at_full-knockout: Rdata files with the results when we allowed several cofactors freely available at full knockout  results

	without_constraint - R files used for the simulations that do NOT constrain fluxes by WT fluxes 
		-results
			PLWT_yeast_without_constraint:Rdata file with a wildtype result at different wildtype flux ratio[ without applying the constraint] 
			PL_NADPH_yeast_without_constraint:Rdata file with the result when NADPH is made freely available at different wildtype 
							  flux ratio[ without applying the constraint] 
			Free_met_at_full-knockout: Rdata files with the results when we allowed several cofactors freely available at full knockout 
						  [without applying the constraint] 

	



pleiotropy_SingleFluxDel -functions for calculating the maximal pleiotropy observed when blocking the individual reactions for which this gene is essential
	R_yeast - R files used for the simulations for yeast  
		-results
	R_Ecoli - R files used for the simulations for Ecoli 
		-results
modularity - functions for calculating modularity
	matlab - Matlab file used for calculating modularity 
	data_plot - data files used for calculating modularity 

plot- R files used for the plot the figures shown in the paper

Supplemental table S1
	-essential biomass components for both Ecoli and yeast

Supplemental table S2
	-currency metabolites taken from (Fritzemeier  et al, 2017)
		-we mapped all of them to be matched with both E. coli and S. cerevisiae
		-in this table you will find the Name of the currency metabolites ,Chemical equation,E. coli model equation,number of Reactions in E. coli, S. 			 cerevisiae model equation and number of Reactions in S. cerevisiae
		
