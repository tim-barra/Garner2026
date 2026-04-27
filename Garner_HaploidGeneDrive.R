##A model of a fungal plant pathogen on a plant crop in a single field
##We will model uninfected (S), infected (I) plants using S and I to match SI model for epidemiology

##Code to run simulations in Garner, Alkemade and Barraclough 2026 'Evaluating the potential to control fungal plant pathogens using gene drive technology'

library(stringr)
library(fBasics)
library(data.table)

##------------------------------------------------------------------------------------------------------------
##Scenario
	scenario<-1

##matrix to store results in
	grand.res<-NULL

	params<-1:99
	yield_params<-matrix(NA,nrow=length(params),ncol=50)
	
	param_effect<-rep(seq(0,0.2,length.out=11),9)
	        
	param_sex<-c(rep(0.4,33),rep(0.9,66))
	param_sexresist<-c(rep(0.05,66),rep(0.85,33))
	param_selfish<-rep(rep(c(10,20,100),each=11),3)
	
for (param in params) {	

##1) Define parameter values 

##Total time to run the simulation
	nseasons<-50
##Maximum number of asexual generations completed on the crop plant each season
	ngens<-5
##Mutation Rate
	mut=c(mut_A=0.001,mut_B=0.001,mut_C=0.001,mut_D=0.001,mut_E=0.001,mut_F=0.00)
##Proportionate yield of infected plants
	Inf_Yield<-0.1
##Toggle whether to additionally plot the within-season dynamics for the final season - by default plots seasons
	plot.within.season<-FALSE

##Wild-type (ABCDEF) parameter values for asexual stages
	##A) Fecundity (Number of Offspring) in soil/debris in spring
	##B) Infection of plants
	##C) Growth on plants
	##D) Proportion of Sex Offspring at end of season
	##E) Winter Survival
	##F) Gene Drive locus

	wt.params<-c(Fecundity=1.3,Infection=0.9,Growth=1.5, Pro_sex=param_sex[param], Survival=0.9,Blank=1)
	wt.alleles<-c("A","B","C","D","E","F")

##Effects of antifungal control on parameter values = selection coefficient of WT in fungicide environment
	control_effect<-c(control_effect_A=0.2,control_effect_B=0.2,control_effect_C=0.2,control_effect_D=0.2,control_effect_E=0.2,control_effect_F=0.2)
##Efficacy of resistance, i.e. how close does resistance get proportionately to restoring original parameter value
	efficacy_resistance<-c(efficacy_resistance_A=0.01,efficacy_resistance_B=0.01,efficacy_resistance_C=0.01,efficacy_resistance_D=0.01,efficacy_resistance_E=0.01,efficacy_resistance_F=0.01) 
##Cost of resistance = selection coefficient of resistant allele in wildtype environment without fungicide
	cost_resistance<-c(cost_resistance_A=0.05,cost_resistance_B=0.05,cost_resistance_C=0.05,cost_resistance_D=param_sexresist[param],cost_resistance_E=0.05,cost_resistance_F=0.05) 
	
	my.model<-5
	
##Gene drive parameters - TIM ADDED INOCULUM SIZE OF GENE DRIVE SPORES
		gene_drive.params<- c(Lethality = 0.0 ,Selfishness = 20, Fecundity_Inhibition = 0.2, Infection_Inhibition = param_effect[param], Growth_Inhibition = 0.2, Sex_Inhibition = 0.2, Survival_Inhibition=0.2, Inoculum = 50)
   	
    ##NB Compare with Burt 2003, epsilon = Selfishness/(Selfishness+1), and Selfishness = epsilon/(1-epsilon)
	##Where epsilon is probability of F converted to f in Homing Endonuclease Gene
	##Selfishness = 9 equates to epsilon = 0.9, default used in Burt paper (10 = 0.91)
	##Selfishness = 20 equates to epsilon = 0.952381. (19 = epsilon 0.95)
	##Selfishness = 99 (or 100 rounding) equates to epsilon =0.99, cited as realistic high value in yeast in Burt


##Assigns when the prescribed control measure is applied
	gene_drive<-array(FALSE,dim=c(nseasons,1))
	fungicide<-array(1,dim=c(nseasons,ngens,5))	

##The set of possible models to consider
	model.choices<-c("No_control","One_Fungi_20","Two_Fungi_Simult","Two_Fungi_Alt","HaploidGeneDrive","DiploidGeneDrive","HaploidFungicide","DiploidFungicide")
##change my.model to select a model

	
	#my.model<-8
	which.model<-model.choices[my.model]

	## a) No control
	 if (which.model=="No_control") {
			}
			
	##Fungicide models from here

	## b) One fungicide, 20 seasons
	 if (which.model=="One_Fungi_20") {
	  which.seasons<-20:40
	  which.fungicide<-1
	  fungicide[which.seasons,, which.fungicide]<-2
    	 }
    	
   	##c) Two fungicides simultaneously
	 if (which.model=="Two_Fungi_Simult") {
	  which.seasons<-20:50
	  which.fungicide<-c(2,3)
	  fungicide[which.seasons,, which.fungicide]<-2
    	}   	

	## d) Alternating fungicides every 5 years
	 if (which.model=="Two_Fungi_Alt") {
	  which.fungicide<-c(2,3)
	  fungicide[c(20:24,30:34,40:44),, which.fungicide[1]]<-2
	  fungicide[c(20:24,30:34,40:44)+5,, which.fungicide[2]]<-2
    	}

	##Gene drive models from here

    ## e) Haploid gene drive - need to set lethality to 0, which parameter is affected and inhibition amount above.
     if (which.model=="HaploidGeneDrive") {
     	##Gene drive parameters - TIM ADDED INOCULUM SIZE OF GENE DRIVE SPORES
		gene_drive.params<- c(Lethality = 0.0 ,Selfishness = param_selfish[param], Fecundity_Inhibition = 0.2, Infection_Inhibition = param_effect[param], Growth_Inhibition = 0.2, Sex_Inhibition = 0.2, Survival_Inhibition=0.2, Inoculum = 50)
		##Selector for which stage the gene_drive is affecting in haploid model
		gene_drive_control<-c(2)
		##when is the gene drive started
		  gene_drive[20:nseasons]<-TRUE
    	}
 
    	
    ##f) Diploid gene drive -  need to set lethality to >0, e.g. 0.99 and inhibition parameters to zero.
  if (which.model=="DiploidGeneDrive") {
     	##Gene drive parameters - TIM ADDED INOCULUM SIZE OF GENE DRIVE SPORES
		gene_drive.params<- c(Lethality = 0.99 ,Selfishness = 20, Fecundity_Inhibition = 0.0, Infection_Inhibition = 0.0, Growth_Inhibition = 0.0, Sex_Inhibition = 0.0, Survival_Inhibition=0.0, Inoculum = 50)
		##Selector for which stage the gene_drive is affecting in haploid model - REQUIRED BUT NOT USED IN THIS VERSION
		gene_drive_control<-c(2)		  
		gene_drive[20:nseasons]<-TRUE
    	}
    
	## g) Haploid gene drive plus one fungicide
   	if (which.model=="HaploidFungicide") {
    	##Gene drive parameters - TIM ADDED INOCULUM SIZE OF GENE DRIVE SPORES
		gene_drive.params<- c(Lethality = 0.0 ,Selfishness = 20, Fecundity_Inhibition = 0.2, Infection_Inhibition = 0.2, Growth_Inhibition = 0.2, Sex_Inhibition = 0.2, Survival_Inhibition=0.2, Inoculum = 50)
		##Selector for which stage the gene_drive is affecting in haploid model
		gene_drive_control<-c(2)
		##when is the gene drive started
		  gene_drive[20:nseasons]<-TRUE
	  		which.seasons<-20:40
	  		which.fungicide<-c(3)
	  		fungicide[which.seasons,, which.fungicide]<-2
    	}
    	
     ##h) Diploid gene drive + fungicide
  if (which.model=="DiploidFungicide") {
     	##Gene drive parameters - TIM ADDED INOCULUM SIZE OF GENE DRIVE SPORES
		gene_drive.params<- c(Lethality = 0.99 ,Selfishness = 20, Fecundity_Inhibition = 0.0, Infection_Inhibition = 0.0, Growth_Inhibition = 0.0, Sex_Inhibition = 0.0, Survival_Inhibition=0.0, Inoculum = 50)
		##Selector for which stage the gene_drive is affecting in haploid model - REQUIRED BUT NOT USED IN THIS VERSION
		gene_drive_control<-c(2)		  
		gene_drive[20:nseasons]<-TRUE
			which.fungicide<-c(2)
	  		fungicide[which.seasons,, which.fungicide]<-2
    	}
   	
    
#----------------------------------------------------------------------------------------------------------------

##2) Construct the mutation table

##TIM JAN 2026- rearranged order of some lines as wasn't working from cold
	Genotype.combos<-expand.grid(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1))
	Genotype.ids<-apply(expand.grid(c("A","a"),c("B","b"),c("C","c"),c("D","d"),c("E","e"),c("F","f")),1,function(x) paste(x,collapse=""))
	##Set up mutation matrix
	mut.list<-list()
	for (i in (1:6)) {
	  mut.list[[i]]<-outer(substring(Genotype.ids,i,i),substring(Genotype.ids,i,i),"==")
	  mut.list[[i]][mut.list[[i]]]<-(1-mut[i])
	  mut.list[[i]][!mut.list[[i]]]<-mut[i]
	}
	mut.matrix<-mut.list[[1]]*mut.list[[2]]*mut.list[[3]]*mut.list[[4]]*mut.list[[5]]*mut.list[[6]]
	

#-----------------------------------------------------------------------------------------
		
##3) Build the parameter array for asexual steps, initially with wildtype values
	##Each matrix in the array contains the associated coefficients for each different locus
	Asex.params<-array(dim=c(2,length(Genotype.ids),length(wt.params)))
	dimnames(Asex.params)[[1]]<-c("None","Control")
	##Row 1 = values when fungicide for that life stage is absent
	##Row 2 = values when fungicide for that life stage is present
	dimnames(Asex.params)[[2]]<-Genotype.ids
	##Columns = each genotype in Genotype.ids
	dimnames(Asex.params)[[3]]<-names(wt.params)
	for (i in (1:length(wt.params))) {Asex.params[,,i]<-wt.params[i]
	Asex.params[1,!grepl(wt.alleles[i],Genotype.ids,fixed=TRUE),i]<-wt.params[i]*(1-cost_resistance[i])}
	for (i in (1:5)) {
	  Asex.params[2, grepl(wt.alleles[i],Genotype.ids,fixed=TRUE),i]<-wt.params[i]*(1-control_effect[i])  
	  ##Values for wild-type allele: first number = value without control, second number = value with control
	  Asex.params[2, !grepl(wt.alleles[i],Genotype.ids,fixed=TRUE),i]<-wt.params[i]*(1-efficacy_resistance[i])    
	  ##Values for alternative/resistant allele: first number = value without control, second number = value with control
	}

##Alter the parameters based on any gene drive effects on haploid life-cycle
	if (any(gene_drive)) {
		for (i in gene_drive_control) {
	  Asex.params[, grepl("f",Genotype.ids,fixed=TRUE),i]<-Asex.params[, grepl("f",Genotype.ids,fixed=TRUE),i]*c(1-gene_drive.params[i+2],1-gene_drive.params[i+2])
	} #Reduces the fitness for the chosen gene_driven effects
	}
	
##Build parameter array for sexual spores
	Sex.params<-Asex.params
	##Currently we assume survival is proportionately higher and fecundity is proportionately lower, irrespective of control
	##TIM JAN 2026- need to think about justification for this a bit, but in any case, 
	##could change so there is a trade-off between survival and fecundity but no net benefit, 
	##plus avoid >100% survival sexual spores, hence i.e. changed to 0.909 and 1.1 to work with defaults. 
	sex_effects<-c(1/1.1,1,1,1,1.1,1)
	Sex.params<-sweep(Sex.params,MARGIN=3,STATS=sex_effects,FUN="*")

#----------------------------------------------------------------

##4) Calculate the random mating and gamete table for sexual reproduction

	alleles<-list()
	for (i in (1:6)) {
	  alleles[[i]]<-strsplit(outer(substring(Genotype.ids,i,i),substring(Genotype.ids,i,i),paste)," ")}
	
	count.genotypes<-function(x) {
	  genos<-outer(alleles[[1]][[x]],alleles[[2]][[x]],paste)
	  genos<-gsub(" ","",as.vector(outer(genos,alleles[[3]][[x]],paste)))
	  genos<-gsub(" ","",as.vector(outer(genos,alleles[[4]][[x]],paste)))
	  genos<-gsub(" ","",as.vector(outer(genos,alleles[[5]][[x]],paste)))
	  genos<-gsub(" ","",as.vector(outer(genos,alleles[[6]][[x]],paste)))
	  return(Genotype.ids%in%genos/sum(Genotype.ids%in%genos))}
	
	RM.table<-t(simplify2array(lapply(1:length(alleles[[1]]),count.genotypes)))
	colnames(RM.table)<-Genotype.ids
	
	allele.array<-array(unlist(alleles),dim=c(2,length(alleles[[1]]),length(alleles)))
	rownames(RM.table)<-paste(apply(allele.array[1,,],1,function(x) paste(x,collapse="")),".",
	                          apply(allele.array[2,,],1,function(x) paste(x,collapse="")),sep="")
	
	##4b) Gene drive alterations to the mating table
	#Rescales f columns to have a greater proportion from sexual reproduction
	RM.table[,grepl("f",Genotype.ids,fixed=TRUE)]<- gene_drive.params["Selfishness"]*RM.table[,grepl("f",Genotype.ids,fixed=TRUE)]
	#Renormalises the table to ensure all rowsums
	RM.table<-RM.table/apply(RM.table,1,sum)
	#Applies lethality to the 'ff' haploids if applicable
	RM.table[str_count(rownames(RM.table),"f")==2,]<-(1-gene_drive.params["Lethality"]^2)*RM.table[str_count(rownames(RM.table),"f")==2,]
	##N.B. For consistency with literature the scalar is 1-e^2 as opposed to simply just redefining the scalar

#---------------------------------------------------------------------------------------------------------

##5) Assign initial values

	##Initial Crop and spore Densities across genotypes
	Initial_Crop<-10000
	#Initial_Spore_Asex<-rep(1.00,length(Genotype.ids))
	#Initial_Spore_Asex[33:64]<-0.0001
	##TIM CHANGED 2026 - start with 100% wild-type alleles
	Initial_Spore_Asex<-rep(0,length(Genotype.ids))
	Initial_Spore_Asex[1]<-100
	Initial_Spore_Sex<-Initial_Spore_Asex
	names(Initial_Spore_Asex)<-Genotype.ids
	names(Initial_Spore_Sex)<-Genotype.ids
##Records of Yield, Load and Spores
	Yield<-array(NA,nseasons)
	Load<-array(NA,nseasons)
	Asex_spores<-matrix(NA,nrow=nseasons,ncol=2^6)
	Sex_spores<-matrix(NA,nrow=nseasons,ncol=2^6)

##Set-up for season 1 from some hypothetical YEAR 0
	 N_Asex <- Initial_Spore_Asex
	 N_Sex<-Initial_Spore_Sex
	  
#--------------------------------------------------------------------------------------------

##6) Run model - step through a set of discrete life stages of fungus within each loop of a yearly cycle##
 
 ##array to plot example of within season dynamics in final season
 within_season<-matrix(NA,nrow=31,ncol=6)
 colnames(within_season)<-c("S","I","Nasex","Nsex","Stage","Gen")

##Start of season loop
 
for (season in (1: nseasons)) {

##A) Start of the year, new field of crop planted, initially all uninfected
	#print(season)
	S<-Initial_Crop
	I<-0
	
	within_season[1,]<-c(S,I,sum(N_Asex),sum(N_Sex),1,1)
	
	#Builds the fungicide set-up for the year
	
	#Sets up the gene_drive if it is enabled for this season TIM ALTERRED 2026 - NOT SURE WHY DOING THIS EVERY LOOP
		##INSTEAD, JUST NEED TO INTRODUCE f ALLELES AT START OF GENE DRIVE
	if (any(gene_drive)){
		if (season == which(gene_drive)[1]) {
		
		N_Asex["ABCDEf"]<-gene_drive.params["Inoculum"]
		
	  #Resets the parameter array
	  #for (i in (1:length(wt.params))) {Asex.params[,,i]<-wt.params[i]
	  	##CHECK HERE - NOT SURE IT IS RIGHT
	  #Asex.params[1,!grepl(wt.alleles[i],Genotype.ids,fixed=TRUE),i]<-wt.params[i]*(1-cost_resistance[i])}
	  #Now applies the fungicide effects
	  #for (i in 1:6) {
	  #  Asex.params[2, grepl(wt.alleles[i],Genotype.ids,fixed=TRUE),i]<-wt.params[i]*(1-control_effect[i])  
	  #  Asex.params[2, !grepl(wt.alleles[i],Genotype.ids,fixed=TRUE),i]<-wt.params[i]*(1-efficacy_resistance[i])}
  	#for (i in gene_drive_control) {
  	 # Asex.params[, grepl("f",Genotype.ids,fixed=TRUE),i]<-Asex.params[, grepl("f",Genotype.ids,fixed=TRUE),i]*c(1-#gene_drive.params[i+2],1-gene_drive.params[i+2])
  	#}
  	
	}
	}
	#Finally sets up the sexual spore parameters - TIM 2026, DON'T NEED THIS?
	#Sex.params<-Asex.params
	#Sex.params<-sweep(Sex.params,MARGIN=3,STATS=sex_effects,FUN="*")
	
	##Record numbers of spores after first planting 
	Asex_spores[season,]<-N_Asex
	Sex_spores[season,]<-N_Sex
	#print(N_Asex)
	#print(N_Sex) 
	
##B) Spores produced by surviving fungal material from the previous year
	N_Asex <- N_Asex * Asex.params[fungicide[season,1,1],,"Fecundity"]
	N_Sex  <- N_Sex  *  Sex.params[fungicide[season,1,1],,"Fecundity"] #Modeled as directly proportional

#Mutation of spores as they reproduce
 	N_Asex <- as.vector(N_Asex%*%mut.matrix)
 	N_Sex <- as.vector(N_Sex%*%mut.matrix)

##TIM JAN 2026- I think we should combine all the spores here for infection
##				otherwise density-dependence is not working properly (each infects more than should
##				if combined together), and no legacy of origin of spores by this stage.

	N_Asex <- N_Asex + N_Sex

	within_season[2,]<-c(S,I,sum(N_Asex),0,1,2)


##C) Initial Infections from the soil 
	#Modeled as proportional to spore density and number of uninfected crops
	
	I_Asex <- N_Asex * Asex.params[fungicide[season,1,2],,"Infection"] * S / (S + I)   ## or S/(S + I) to match SIR models
	## COMMENT OUT JAN 2026 I_Sex <- N_Sex * Sex.params[fungicide[season,1,2],,"Infection"] * S/(S + I)   ## or S/(S + I) to match SIR models
	# Calculates the new infections
	
	New_I<- I_Asex #Genotype array of the new infections
	I <- sum(New_I) #Total number of new infections
	S <- S - I #Remaining healthy crop
	
	within_season[3,]<-c(S,I, sum(New_I),0,2,3)

	
##D) Asexual sporulation and infection generations occur as fungus spreads on the crop

##WITHIN THIS LOOP ASSUME THAT ALL REPRODUCTION IS ASEXUAL - SO NO SEXUAL SPORES FORMED DURING THE GEN LOOP OR PRESENT AT THIS TIME
##NEED PARAMETER/FITNESS THAT DEPENDS ON PRESENCE/ABSENCE OF CONTROL MEASURE, SO CAN TOGGLE ON/OFF OR ALTERNATE 

for (gen in (1:ngens)) {
#Production of new spores per infected plant
	N_Asex <- New_I * Asex.params[fungicide[season,gen,3],,"Growth"]
#Mutation of spores as they reproduce
 	N_Asex <- as.vector(N_Asex%*%mut.matrix)
	
#New infections arise from asexual spores	
  New_I <- N_Asex * Asex.params[fungicide[season,gen,2],,"Infection"] * S/(S + I)   ## or S/(S + I) to match SIR models
	I <- I + sum(New_I)
	S <- S - sum(New_I)
	
	within_season[3+gen,]<-c(S,I, sum(New_I),0,3,3+gen)

	
	} ##end of generations loop
	
	##No need to update these during the Asexual generations, but still need the correct values for the fungicide state at the end 
	Sex.params<-Asex.params
	Sex.params<-sweep(Sex.params,MARGIN=3,STATS=sex_effects,FUN="*")
##E) Crop Harvesting and Load Recording

	Yield[season]<- S + Inf_Yield * I
	Load[season]<-I/Initial_Crop
	
	within_season[3+ngens+1,]<-c(0,0, sum(New_I),0,4,3+ngens+1)

	
##F) Some fungal material produced by last infection cycle overwinters as sexual or as asexual, depending on genotype

	N_Asex <- New_I * (1-Asex.params[fungicide[season,ngens,4],,"Pro_sex"])
	N_Sex <- New_I * (Asex.params[fungicide[season,ngens,4],,"Pro_sex"])
	
	within_season[3+ngens+2,]<-c(0,0, sum(N_Asex),sum(N_Sex),5,3+ngens+2)

	
##G) Sexual Reproduction at the end of the season
#N.B Final sum of density will be same, but ratios of genotypes could change within N_Sex
	
	tmp.Sex_Spore.prop<-N_Sex%o%N_Sex
	tmp.Sex_Spore.prop<-tmp.Sex_Spore.prop/sum(tmp.Sex_Spore.prop)
	Sex_Spore.prop<-colSums(RM.table*as.vector(tmp.Sex_Spore.prop))
	N_Sex<-sum(N_Sex)* Sex_Spore.prop
##Mutations occur during spore formation 
 	N_Sex <- as.vector(N_Sex%*%mut.matrix)	

 	N_Asex <- as.vector(N_Asex%*%mut.matrix)	

##H) Over-winter survival  
	
 	N_Asex <- N_Asex * Asex.params[fungicide[season,ngens,5],,"Survival"]  
 	N_Sex <- N_Sex * Sex.params[fungicide[season,ngens,5],,"Survival"] #Currently modelled as directly proportional
 	
 	names(N_Asex)<-Genotype.ids 
 	names(N_Sex)<-Genotype.ids 
 	
 	within_season[3+ngens+3,]<-c(0,0, sum(N_Asex),sum(N_Sex),6,3+ngens+3)
 		  	
	} ##end of season loop

##add column names
colnames(Asex_spores)<-Genotype.ids
colnames(Sex_spores)<-Genotype.ids

#------------------------------------------------------------------------------------------------------------

##7) PLOT DENSITIES AND ALLELE FREQUENCIES	

##set color-blind palette for loci
col.pal<-c("#D55E00","#E69F00","#F0E442","#009E73","#0072B2","#CC79A7")
cx<-1.2

##Make filename so that keep figures for different versions
tmp<-""
if (exists("which.fungicide")&(max(fungicide)==2)) tmp<-paste0(".",paste0(which.fungicide,collapse=""))
if (exists("gene_drive_control")&any(gene_drive)) tmp<-paste0(tmp,".", paste0( gene_drive_control,collapse=""))

file.name<-paste0("Figures/Model",my.model,"_param",param,tmp,".pdf")

pdf(file=file.name,height=5.5,width=2.5)

par(mfcol=c(3,1),mar=c(2,4,2,1),oma=c(2,0,0,0))

##PLOT OVER SEASONS
comb.spores<-Asex_spores+ Sex_spores
plot(Yield/Initial_Crop,type="l",main="",xlab='Season',ylab="% maximum",las=1,ylim=c(0,1),cex.axis= cx,cex.lab= cx)
lines(Load,type="l",main="Load",xlab='Season',las=1,lty=3)
legend("bottomleft", bty="n",inset=.05, legend=c("Yield","Load"),lty=c(1,3),cex= cx)
plot(rowSums(Asex_spores),ylim=c(0,max(c(rowSums(Asex_spores),rowSums(Sex_spores)))),col="darkblue",type="l",ylab="Fungal density",main="",xlab='Season',las=1,cex.axis= cx,cex.lab= cx)
lines(rowSums(Sex_spores),col="red")
legend("bottomleft", bty="n",inset=.05,legend=c("Asex","Sex"),col=c('darkblue','red'),lty=1,cex= cx)
##combine genotype counts of asexual and sexual spores
##CALCULATE AND PLOT ALLELE FREQUENCES
allele.freqs<-matrix(NA,nrow=nseasons,ncol=length(wt.alleles))
for (i in (1:length(wt.alleles))) {
  allele.freqs[,i]<-rowSums(comb.spores[,grepl(wt.alleles[i],colnames(comb.spores))])/rowSums(comb.spores)
}
dimnames(allele.freqs)[[2]]<-wt.alleles
matplot(allele.freqs,type="l",lty=1,las=1,main="",col=col.pal,xlab='Season',ylab='Allele Proportion',cex.axis=cx,cex.lab=cx)
legend("bottomleft", bty="n",inset=.05,legend=wt.alleles,fill=col.pal,horiz=FALSE)
mtext(text = "Season",side=1,line=2.5,cex=0.7*cx)

dev.off()


# 7a) Optional, needs switching on - PLOT WITHIN SEASON


file.name2<-gsub("Model","WithinSeason",file.name)

if (plot.within.season) {

  pdf(file= file.name2, height=5.5,width=2.5)
  
  par(mfcol=c(3,1),mar=c(2,4,2,1),oma=c(2,0,0,0))
  matplot(within_season[1:11,1:2],type="s",ylab="Plant Density",xlab="step",main="",col=1,lty=c(1,3),xaxt="n",cex.axis=cx,cex.lab=cx)
  legend("topright", bty="n",inset=.05,legend=c("Susceptible","Infected"),lty=c(1,3))
  polygon(c(1,2,2,1),c(0,0,max(within_season[1:11,1:2])*1.1, max(within_season[1:11,1:2]))*1.1,col=gray(0.5,0.1),border=NA)
  polygon(c(3,8,8,3),c(0,0,max(within_season[1:11,1:2])*1.1, max(within_season[1:11,1:2]))*1.1,col=gray(0.5,0.1),border=NA)
  polygon(c(9,10,10,9),c(0,0,max(within_season[1:11,1:2])*1.1, max(within_season[1:11,1:2]))*1.1,col=gray(0.5,0.1),border=NA)
  mtext(text=1:6,at=c(1.5,2.5,5.5,8.5,9.5,10.5),side=1,line=0.2,cex=0.7*cx)
  mtext(text="Stage",side=1,line=1.2,cex=0.7*cx)
  
  matplot(within_season[1:11,3:4],type="s",ylab="Fungal Density",xlab="step",main="",col=c("darkblue","red"),lty=1,xaxt="n",cex.axis=cx,cex.lab=cx)
  legend("center", bty="n",inset=.05, legend=c("Asex","Sex"),col=c('darkblue','red'),lty=1,cex=cx)
  polygon(c(1,2,2,1),c(0,0,max(within_season[1:11,3:4])*1.1, max(within_season[1:11,3:4]))*1.1,col=gray(0.5,0.1),border=NA)
  polygon(c(3,8,8,3),c(0,0,max(within_season[1:11,3:4])*1.1, max(within_season[1:11,3:4]))*1.1,col=gray(0.5,0.1),border=NA)
  polygon(c(9,10,10,9),c(0,0,max(within_season[1:11,3:4])*1.1, max(within_season[1:11,3:4]))*1.1,col=gray(0.5,0.1),border=NA)
  mtext(text=1:6,at=c(1.5,2.5,5.5,8.5,9.5,10.5),side=1,line=0.2,cex=0.7*cx)
  mtext(text="Stage",side=1,line=1.2,cex=0.7*cx)
  
  ##check fungal densities include value immediately prior to asex/sex reproduction
  ##add generations/stages to the plot instead of steps
  dev.off()

}

yield_params[param,20:50]<-Yield[20:50]/Initial_Crop

#----
##CALCULATE EQUILIBRIUM SOLUTIONS FOR ECOLOGICAL CASE

##Tried doing this in SageMath - but it failed to find a solution for Nasex at equilibrium. 

#-------------------------------------------------------------------------------------------------------------------
##SUMMARY STATISTICS

Jumps<-c(0)
Behaviour<-array(dim=c(nseasons,3))
for (i in 2:(nseasons)){Jumps[i]<-(Yield[i]-Yield[i-1])}

##Behaviour Categorising - Cycles, Spikes, Converged, Decreasing
for (i in 2:nseasons) {
  if (abs(Jumps[i])<2) {Behaviour[i,1]<-'Converged'}
  if (Jumps[i]>2) {Behaviour[i,1]<-'Increasing'}
  if (Jumps[i]< -2) {Behaviour[i,1]<-'Decreasing'}
  if ((Jumps[i]*Jumps[i-1])<0) {Behaviour[i,2]<-'Cyclic'}
  else {Behaviour[i,2]<-'Monotone'}
  if (abs(Jumps[i])>2000) {Behaviour[i,3]<-'Spike'}
}

Data <- data.frame(yield = Yield,
                   differences = Jumps,
                   load = Load,
                   behaviour = Behaviour,
                   asex_total = rowSums(Asex_spores),
                   sex_total = rowSums(Sex_spores),
                   asex_Densities = Asex_spores,
                   sex_Densities = Sex_spores,
                   allele_Freqs = allele.freqs,
                   
                   row.names = 1:nseasons)
Summary_Data<-basicStats(Data[c('yield','load','asex_total','sex_total')])[c("Mean","Minimum","Maximum","Stdev"),]
#matplot(frollmean(Data[,"yield"],2),type="l",lty=1,main="Yield Rolling Avg",ylabel='Two step average')

Data[c('asex_Densities.ABCDEf','sex_Densities.ABCDEf','asex_Densities.ABCDEF','sex_Densities.ABCDEF')]
##Parameter Summary

Parameters <- data.frame(Mutation_Rates = mut,
                         W.T.Parameters = wt.params,
                         Sexual_Effects = sex_effects,
                         Gene_Drive_Effects = gene_drive.params[3:8],
                         Costs = cost_resistance,
                         Efficacies = efficacy_resistance,
                         Control_Effects = control_effect,
                         row.names= wt.alleles)
    Parameters<-Parameters[1:5,]
    Parameters<-Parameters[,-4]
                         

yield_res<-c(pre_yield= mean(Yield[5:20]), peak_yield= max(Yield[5:20]), yield_40=Yield[40], yield_50=Yield[50],mean_21_40=mean(Yield[21:40]), mean_41_50=mean(Yield[41:50]), mean_21_50=mean(Yield[21:50]))/Initial_Crop

Gene_drive<-array(0,5)
if (any(gene_drive)) {Gene_drive[gene_drive_control]<-1}
names(Gene_drive)<-paste0("GeneDrive_",1:5)

##set to NA for any that weren't used
gene_drive.params[(3:7)[which(Gene_drive==0)]]<-NA

Fungi_cont<-array(0,5)
if (any(fungicide==2)) {Fungi_cont[which.fungicide]<-1}
names(Fungi_cont)<-paste0("Fungicide_",1:5)

##set to NA for any that weren't used
Parameters[which(Fungi_cont==0),5:6]<-NA

##pull out allele frequencies at these seasons
allele.freqs20<-allele.freqs[20,]
names(allele.freqs20)<-paste("Freq20_",colnames(allele.freqs),sep="")
allele.freqs40<-allele.freqs[40,]
names(allele.freqs40)<-paste("Freq40_",colnames(allele.freqs),sep="")
allele.freqs50<-allele.freqs[50,]
names(allele.freqs50)<-paste("Freq50_",colnames(allele.freqs),sep="")

##single line summary of results
res<-c(scenario=scenario,model=my.model,unlist(Parameters),Fungi_cont,Num_fung=sum(Fungi_cont),gene_drive.params, Gene_drive,Num_drive=sum(Gene_drive),Infection_Yield=Inf_Yield,Number_of_Generations=ngens,
	yield_res, allele.freqs20 , allele.freqs40 ,allele.freqs50)
res<-t(res)

grand.res<-rbind(grand.res,res)

}  ##end my.model loop

	file.name4<-paste0("Figures/Model",my.model,"_params",tmp,".pdf")
	
	pdf(file=file.name4,height=6,width=6)
	
	par(mfcol=c(3,3),mar=c(2,2,2,1),oma=c(3,3,0,0),xpd=NA)
	
	##Mainly asex
	
	subs<-(param_sex==0.4)&(param_selfish==10)
	plot(param_effect[subs],1-grand.res[subs,"Freq50_F"],ylab="Frequency drive allele",xlab="",type="l",las=1,cex.axis=cx,cex.lab=cx,lty=3)
	subs<-(param_sex==0.4)&(param_selfish==20)
	lines(param_effect[subs],1-grand.res[subs,"Freq50_F"],type="l",lty=2)
	subs<-(param_sex==0.4)&(param_selfish==100)
	lines(param_effect[subs],1-grand.res[subs,"Freq50_F"],type="l",lty=1)

	subs<-(param_sex==0.4)&(param_selfish==10)
	plot(param_effect[subs],grand.res[subs,"mean_21_50"],ylim=c(0.5,1.0),ylab="% Yield",xlab="",type="l",las=1,cex.axis=cx,cex.lab=cx,lty=3)
	subs<-(param_sex==0.4)&(param_selfish==20)
	lines(param_effect[subs],grand.res[subs,"mean_21_50"],type="l",lty=2)
	subs<-(param_sex==0.4)&(param_selfish==100)
	lines(param_effect[subs],grand.res[subs,"mean_21_50"],type="l",lty=1)
	legend("left", bty="n",inset=.05,legend=c("100x","20x","10x"),lty=c(1:3),cex=cx)

	subs<-(param_sex==0.4)&(param_selfish==10)
	plot(param_effect[subs],grand.res[subs,"Freq50_D"],ylim=c(0,1.0),ylab="Frequency sex allele",xlab="q (effect on infection rate)",type="l",las=1,cex.axis=cx,cex.lab=cx,lty=3)
	subs<-(param_sex==0.4)&(param_selfish==20)
	lines(param_effect[subs],grand.res[subs,"Freq50_D"],type="l",lty=2)
	subs<-(param_sex==0.4)&(param_selfish==100)
	lines(param_effect[subs],grand.res[subs,"Freq50_D"],type="l",lty=1)

	##Mainly sex
	
	subs<-(param_sex==0.9)&(param_selfish==10)&(param_sexresist==0.05)
	plot(param_effect[subs],1-grand.res[subs,"Freq50_F"],ylab="",xlab="",type="l",las=1,cex.axis=cx,cex.lab=cx,lty=3)
	subs<-(param_sex==0.9)&(param_selfish==20)&(param_sexresist==0.05)
	lines(param_effect[subs],1-grand.res[subs,"Freq50_F"],type="l",lty=2)
	subs<-(param_sex==0.9)&(param_selfish==100)&(param_sexresist==0.05)
	lines(param_effect[subs],1-grand.res[subs,"Freq50_F"],type="l",lty=1)
	
	subs<-(param_sex==0.9)&(param_selfish==10)&(param_sexresist==0.05)
	plot(param_effect[subs],grand.res[subs,"mean_21_50"],ylim=c(0.5,1.0),ylab="",xlab="",type="l",las=1,cex.axis=cx,cex.lab=cx,lty=3)
	subs<-(param_sex==0.9)&(param_selfish==20)&(param_sexresist==0.05)
	lines(param_effect[subs],grand.res[subs,"mean_21_50"],type="l",lty=2)
	subs<-(param_sex==0.9)&(param_selfish==100)&(param_sexresist==0.05)
	lines(param_effect[subs],grand.res[subs,"mean_21_50"],type="l",lty=1)

	subs<-(param_sex==0.9)&(param_selfish==10)&(param_sexresist==0.05)
	plot(param_effect[subs],grand.res[subs,"Freq50_D"],ylim=c(0,1.0),ylab="",xlab="q (effect on infection rate)",type="l",las=1,cex.axis=cx,cex.lab=cx,lty=3)
	subs<-(param_sex==0.9)&(param_selfish==20)&(param_sexresist==0.05)
	lines(param_effect[subs],grand.res[subs,"Freq50_D"],type="l",lty=2)
	subs<-(param_sex==0.9)&(param_selfish==100)&(param_sexresist==0.05)
	lines(param_effect[subs],grand.res[subs,"Freq50_D"],type="l",lty=1)

	##Mainly sex
	
	subs<-(param_sex==0.9)&(param_selfish==10)&(param_sexresist==0.85)
	plot(param_effect[subs],1-grand.res[subs,"Freq50_F"],ylab="",xlab="",type="l",las=1,cex.axis=cx,cex.lab=cx,lty=3)
	subs<-(param_sex==0.9)&(param_selfish==20)&(param_sexresist==0.85)
	lines(param_effect[subs],1-grand.res[subs,"Freq50_F"],type="l",lty=2)
	subs<-(param_sex==0.9)&(param_selfish==100)&(param_sexresist==0.85)
	lines(param_effect[subs],1-grand.res[subs,"Freq50_F"],type="l",lty=1)
	
	subs<-(param_sex==0.9)&(param_selfish==10)&(param_sexresist==0.85)
	plot(param_effect[subs],grand.res[subs,"mean_21_50"],ylim=c(0.5,1.0),ylab="",xlab="",type="l",las=1,cex.axis=cx,cex.lab=cx,lty=3)
	subs<-(param_sex==0.9)&(param_selfish==20)&(param_sexresist==0.85)
	lines(param_effect[subs],grand.res[subs,"mean_21_50"],type="l",lty=2)
	subs<-(param_sex==0.9)&(param_selfish==100)&(param_sexresist==0.85)
	lines(param_effect[subs],grand.res[subs,"mean_21_50"],type="l",lty=1)
	
	subs<-(param_sex==0.9)&(param_selfish==10)&(param_sexresist==0.85)
	plot(param_effect[subs],grand.res[subs,"Freq50_D"],ylim=c(0,1.0),ylab="",xlab="q (effect on infection rate)",type="l",las=1,cex.axis=cx,cex.lab=cx,lty=3)
	subs<-(param_sex==0.9)&(param_selfish==20)&(param_sexresist==0.85)
	lines(param_effect[subs],grand.res[subs,"Freq50_D"],type="l",lty=2)
	subs<-(param_sex==0.9)&(param_selfish==100)&(param_sexresist==0.85)
	lines(param_effect[subs],grand.res[subs,"Freq50_D"],type="l",lty=1)
	
		
	dev.off()
	
	
	#file.name3<-paste0("Scenario",scenario,".csv")

write.csv(grand.res,file="HaploidGeneDrive.csv",row.names=F)