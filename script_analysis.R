##############################################
# Open the data
##############################################
#setwd("/Users/DGravel/Documents/Projects_On_Going/Niche_allometrique/analyse_I")
setwd("/Users/DGravel/Documents/Manuscripts/2012/Niche_MEE_Sept2012/analysis_II")
source("get_pars_Niche.R")

data = read.table("pp.txt")[,c(9,20,18,28,2,11,22)]
data[,5] = as.numeric(data[,5])
data[,1] = as.character(data[,1])
data[,2] = as.character(data[,2])
data[,6] = as.character(data[,6])
data[,7] = as.character(data[,7])

un_named_1 = data[,1]=="-999"
data[un_named_1,1] = data[un_named_1,6]
un_named_2 = data[,2]=="-999"
data[un_named_2,2] = data[un_named_2,7]

##############################################
# Figure 1: example
##############################################
i = 2
subdata = subset(data,data[,5]==i)
nprey = as.character(subdata[,1])
npred = as.character(subdata[,2])	
Bpred = log(subdata[,3],10)
Bprey = log(subdata[,4],10)

B_all = c(Bprey,Bpred)
dup = duplicated(B_all)
B = subset(B_all,dup==FALSE)
sp_all = c(nprey,npred)
sp = subset(sp_all,dup==FALSE)

library(quantreg)
mean_reg = lm(Bprey~Bpred)
qr95 = rq(Bprey~Bpred,tau = .95)
qr05 = rq(Bprey~Bpred,tau = .05)

quartz(height = 4.5, width = 3.5*1.5)
par(mar=c(5,4.5,2,1))
plot(Bpred,Bprey,xlab = "Log10 Body size of predators", ylab = "Log10 Body size of preys",xlim=range(Bpred),cex.lab =1.25, cex.axis=1.)
abline(mean_reg)
abline(qr95, col = "grey")
abline(qr05, col = "grey")

delta = mean_reg$coef[2]
b1 = mean_reg$coef[1]
b2 = delta	
n = B
c = b1 + b2*B
high = qr95$coef[1] + qr95$coef[2]*B
low = qr05$coef[1] + qr05$coef[2]*B

S = length(n)   	
pairs_niche = NULL
for(i in 1:S)
   	for(j in 1:S)
   	  if(n[i] > min(Bpred) & n[j]>low[i] && n[j]<high[i]) pairs_niche = rbind(pairs_niche,c(B[j],B[i]))
    	  
points(pairs_niche[,2],pairs_niche[,1],pch = 19, cex = 0.3)
setwd("/Users/DGravel/Desktop/ms_NicheAllometric/allometric_v2")
dev.copy2pdf(file = "Example.pdf")


##############################################
# Figure 2. Calculate TSS for all webs
##############################################

results = matrix(nr = 15, nc = 5)
stats = matrix(nr=15,nc=2)

# Loop over all food webs
for(web in 1:15) {
		
	# Prepare data
	subdata = subset(data,data[,5]==web)
	Bpred = log(subdata[,3],10)
	Bprey = log(subdata[,4],10)
	nprey = as.character(subdata[,1])
	npred = as.character(subdata[,2])	

	B_all = c(Bprey,Bpred)
	dup = duplicated(B_all)
	B = subset(B_all,dup==FALSE)
	sp_all = c(nprey,npred)
	sp = subset(sp_all,dup==FALSE)
		
	# Calculate the parameters
	pars = reg_Niche(Bprey,Bpred)
	pars_all = get_pars_Niche(pars,B)
		
	# Calculate the food web structure
	Lobs = Lobs_fn(Bprey,Bpred)
	Lpred = Lpred_fn(pars_all[,1],pars_all[,2],pars_all[,3],pars_all[,4],min(Bpred))
		
	# Calculate the TSS
	TSS = TSS_fn(Lpred,Lobs)
		
	# Record results
	results[web,1] = pars[[4]]
	results[web,2:5] = TSS
	stats[web,1]=length(B)
	stats[web,2]=length(Bpred)
	cat(web," S=",length(B)," L=",length(Bpred),'\n')	
}

r2 = results[-c(6,7,10,13),1]
a  = results[-c(6,7,10,13),2]
b  = results[-c(6,7,10,13),3]
c  = results[-c(6,7,10,13),4]
d  = results[-c(6,7,10,13),5]

# TSS
TSS = (a*d-b*c)/((a+c)*(b+d))

# Taux de succès des prédictions
FracPred = (a+d)/(a+b+c+d)

# Partition des fractions

Frac_a = (a)/(a+b+c+d)
Frac_d = (d)/(a+b+c+d)
Frac_bc = (b+c)/(a+b+c+d)

# Liens observés non prédits
Frac_c = (c)/(a+b+c+d)

quartz(height = 4.5, width = 3.5*1.5*2)
par(mar=c(6,6,2,1), mfrow = c(1,2),mgp=c(3,1,0))
plot(r2,TSS,xlab = "R2 of the Pred:Prey relationship",ylab = "TSS",cex.axis = 1., cex.lab = 1.25,ylim=c(0,1))
abline(lm(TSS~r2))
mtext(text="A)",side=3,line=0.5,adj=-0.1,cex=1.25)

plot(r2,Frac_a,xlab = "R2 of the Pred:Prey relationship",ylab = "Fractions",cex.axis = 1, cex.lab = 1.25,ylim=c(0,1),pch = 19)
mtext(text="B)",side=3,line=0.5,adj=-0.1,cex=1.25)
points(r2,Frac_bc,pch=1,cex=1)
points(r2,Frac_d,pch=9,cex=1)
abline(lm((Frac_a)~r2),lty=1)
abline(lm(Frac_bc~r2),lty=2)
abline(lm(Frac_d~r2),lty=3)
legend("topleft",pch=c(19,1,9),legend = c("a","b+c","d"),bty="n",cex=1,lty = c(1,2,3))

setwd("/Users/DGravel/Desktop/ms_NicheAllometric/allometric_v2")
dev.copy2pdf(file = "TSS_perweb.pdf")


##############################################
# Figure 3. Run the analysis with random subsampling
##############################################

random_subsample = function(Bprey,Bpred,nprey,npred,n_min=5,n_repl=100){
	## Output
	n_sp_rem = NULL
	p_sp_rem = NULL
	TSS_rem = NULL
	## Number of species to remove
	for(n_rem in n_min:length(Bprey)){
		for(repl in 1:n_repl){
			Idx_to_keep = sample(c(1:length(Bprey)),n_rem,replace=FALSE)
			t_Bprey = Bprey[Idx_to_keep]
			t_Bpred = Bpred[Idx_to_keep]
			t_nprey = nprey[Idx_to_keep]
			t_npred = npred[Idx_to_keep]

			## Get the parameters from the subset of species pairs
			pars = reg_Niche(t_Bprey,t_Bpred)

			## Evaluate the TSS for the original network
			B_all = c(Bprey,Bpred)
			B = subset(B_all,duplicated(B_all)==FALSE)
			sp_all = c(nprey,npred)
			sp = subset(sp_all,duplicated(B_all)==FALSE)		
			pars_all = get_pars_Niche(pars,B)
			Lpred = Lpred_fn(pars_all[,1],pars_all[,2],pars_all[,3],pars_all[,4],min(t_Bpred))
			Lobs = Lobs_fn(Bprey,Bpred)
			res = TSS_fn(Lpred,Lobs)
			a  = res[1]
			b  = res[2]
			c  = res[3]
			d  = res[4]
			current_TSS = (a*d-b*c)/((a+c)*(b+d))

			## Data output
			TSS_rem = c(TSS_rem,current_TSS)
			n_sp_rem = c(n_sp_rem,length(Bprey)-n_rem)
			p_sp_rem = c(p_sp_rem,(length(Bprey)-n_rem)/length(Bprey))
		}
	}
	return(list(S=n_sp_rem,P=p_sp_rem,TSS=TSS_rem))
}

web = 2
subdata = subset(data,data[,5]==web)
nprey = as.character(subdata[,1])
npred = as.character(subdata[,2])	
Bprey = log(subdata[,4],10)
Bpred = log(subdata[,3],10)
	
r_subsample = random_subsample(Bprey,Bpred,nprey,npred,n_min=5,n_repl=100)
Means = aggregate(r_subsample,by=list(P=r_subsample$P),mean)
Sds = aggregate(r_subsample,by=list(P=r_subsample$P),sd)

quartz(height = 4.5, width = 3.5*1.5)
par(mar=c(4.5,4.5,1,1))
plot(Means$P,Means$TSS,ylim=c(0,1),lty=1,lwd=2,type='l',xlab = "Fraction of links removed", ylab = "TSS",cex.lab =1.25)
lines(Means$P,Means$TSS-Sds$TSS,col = "grey")
lines(Means$P,Means$TSS+Sds$TSS, col = "grey")
setwd("/Users/DGravel/Desktop/ms_NicheAllometric/allometric_v2")
dev.copy2pdf(file = "TSS_sampling.pdf")



##############################################
# Figure 4. Application to the Med
##############################################
setwd("/Users/DGravel/Documents/Projects_On_Going/Niche_allometrique/analyse_I")
source("get_pars_Niche.R")

get_pairs = function(L,BodySize) {
	S = nrow(BodySize)
	pairs = NULL
	for(i in 1:S)
		for(j in 1:S)
			if(L[i,j]==1) pairs = rbind(pairs,c(BodySize[i,],BodySize[j,]))

	return(pairs) 
}

L_bonif = read.table("L_bonif.txt")
BodySize_bonif = read.table("BodySize_bonif.txt")
pairs_bonif = get_pairs(L_bonif,BodySize_bonif)

L_martha = read.table("L_martha.txt")
BodySize_martha = read.table("BodySize_martha.txt")
pairs_martha = get_pairs(L_martha,BodySize_martha)

L_italy = read.table("L_italy.txt")
BodySize_italy = read.table("BodySize_italy.txt")
pairs_italy = get_pairs(L_italy,BodySize_italy)

pairs_all = rbind(pairs_bonif,pairs_martha,pairs_italy)

pars = reg_Niche(log(pairs_all[,1],10),log(pairs_all[,2],10))

BodySize_all0 = log(read.table("BodySize_med.txt",header=T)[,9],10)
BodySize_rank = read.table("BodySize_med.txt",header=T)[,12]
BodySize_all = sort(BodySize_all0)
pars_all = get_pars_Niche(pars,BodySize_all)
		
# Remove links from the matrix for rare co-occurrences
load("Cooc_mat")
ID0 = row.names(Cooc_mat)
ID0 = as.numeric(substring(ID0,5,10))
S = max(ID0)
mat0 = Cooc_mat
mat1 = matrix(999,nr = S, nc = S)
for(i in 1:S) 
	for(j in 1:S) if(is.element(i,ID0) == TRUE & is.element(j,ID0) == TRUE) 
			mat1[i,j] = mat0[which(ID0==i),which(ID0==j)]		
IDsort = sort(ID0)
mat1 = mat1[c(1:S)%in%ID0==TRUE,c(1:S)%in%ID0==TRUE]	
mat1 = mat1[-c(64,102),-c(64,102)]

mat2 = mat1*0
for(i in 1:nrow(mat1)) 
	for(j in 1:nrow(mat1)) mat2[i,j] = mat1[which(BodySize_rank==i),which(BodySize_rank==j)]

# Calculate the food web structure
Lpred = Lpred_fn(pars_all[,1],pars_all[,2],pars_all[,3],pars_all[,4],min(log(pairs_all[,1],10)))
		
thresh = 0.001*max(mat2)
#thresh = 1
Lpred_real = Lpred
Lpred_real[mat2<thresh] = 0

# Plot the web
quartz(height = 4.5, width = 3.5*1.5*2)
par(mar=c(5,4.5,2,1), mfrow=c(1,2))

image(x=c(1:length(BodySize_all)),y=c(1:length(BodySize_all)),z=Lpred,xlab = "Predator body size rank", ylab = "Prey body size rank",cex.axis = 1,cex.lab = 1.25,col=c("black","white"))
mtext(text="A)",side=3,line=0.5,adj=-0.1,cex=1.25)

image(x=c(1:length(BodySize_all)),y=c(1:length(BodySize_all)),z=Lpred_real,xlab = "Predator body size rank", ylab = "Prey body size rank",cex.axis = 1,cex.lab = 1.25,col=c("black","white"))
mtext(text="B)",side=3,line=0.5,adj=-0.1,cex=1.25)

setwd("/Users/DGravel/Desktop/ms_NicheAllometric/allometric_v2")
dev.copy2pdf(file = "Exemple_MED.pdf")

##############################################
# Figure 5. Impact of bodysize distribution on degree distribution
##############################################
setwd("/Users/DGravel/Desktop/ms_NicheAllometric/analysis_II/Med")
MW_Med = read.table("MW_Med.txt")
BS_Med = as.matrix(read.table("BS_Med.txt"))

setwd("/Users/DGravel/Desktop/ms_NicheAllometric/analysis_II")
source("get_pars_Niche.R")

pairs_all = matrix(nr = sum(MW_Med), nc = 2)

S = nrow(MW_Med)
n  = 1
for(i in 1:S) for(j in 1:S) if(MW_Med[i,j]==1) {
	pairs_all[n,1] = BS_Med[i]
	pairs_all[n,2] = BS_Med[j]
	n = n+1
	}
 
pars = reg_Niche(log(pairs_all[,1],10),log(pairs_all[,2],10))

BS075 = BS_Med*0.6
Bsub = BS_Med
Bsub[BS_Med>quantile(BS_Med,0.6)] = 0

pars075 = get_pars_Niche(pars,BS075)
parssub = get_pars_Niche(pars,Bsub)

Lpred075 = Lpred_fn(pars075[,1],pars075[,2],pars075[,3],pars075[,4],min(BS_Med))
Lpred075[mat2<thresh] = 0

Lpredsub = Lpred_fn(parssub[,1],parssub[,2],parssub[,3],parssub[,4],min(BS_Med))
Lpredsub[mat2<thresh] = 0

Dall = (apply(MW_Med,2,sum)+apply(MW_Med,1,sum))
Dall = Dall[Dall>0]
D075 = (apply(Lpred075,2,sum)+apply(Lpred075,1,sum))
D075 = D075[D075>0]
Dsub = (apply(Lpredsub,2,sum)+apply(Lpredsub,1,sum))
Dsub = Dsub[Dsub>0]

cumprobAll = 1-cumsum(as.matrix(table(Dall)))/length(Dall)+(cumsum(as.matrix(table(Dall)))/length(Dall))[1]
DegAll = as.numeric(names(table(Dall)))

cumprob075 = 1-cumsum(as.matrix(table(D075)))/length(D075)+(cumsum(as.matrix(table(D075)))/length(D075))[1]
Deg075 = as.numeric(names(table(D075)))

cumprobsub = 1-cumsum(as.matrix(table(Dsub)))/length(Dsub)+(cumsum(as.matrix(table(Dsub)))/length(Dsub))[1]
Degsub = as.numeric(names(table(Dsub)))

quartz(height = 4.5, width = 3.5*1.5)
par(mar=c(4.5,4.5,1,1))
plot(DegAll,cumprobAll,type="l",log="y",lwd=1.5,ylim=range(cumprobAll,cumprob075),xlim=range(DegAll,Deg075),xlab = "Degree", ylab = "Cumulative distribution",cex.lab = 1.25,cex.axis = 1.)
lines(Deg075,cumprob075,lty=3,lwd=1.5)
lines(Degsub,cumprobsub,lty=2,lwd=1.5)

legend("bottomleft",bty="n",lty=c(1,3,2),lwd = 1.5, legend = c("All species", "40% reduction in BS", "Removal of 40% of largest BS"))

setwd("/Users/DGravel/Desktop/ms_NicheAllometric/allometric_v2")
dev.copy2pdf(file = "Degree.pdf")




