
###### Evolution of rostral scale and mimicry in the genus Xenodon Boie, 1826 (Serpentes, Dipsadidae, Xenodontinae) #####
################################# Hugo Cabral, Pier Cacciali & Diego J. Santana #########################################
# Library
library(caper)
library(phytools)
library(geiger)
library(ade4)
library(picante)
library(phytools)
library(caper)
library(ape)
rm(list=ls()); gc()
setwd("")

# Load table with information of characters
X<-read.csv("",row.names=1, sep=";")
rostral<-setNames(X[,1],rownames(X))
mimic<-setNames(X[,2], rownames(X))
habitat <-setNames(X[,4], rownames(X))

# Load tree
Tree<-read.tree("")
Tree
head(xenodon_out)
plot(Tree, lwd=4)
axisPhylo()

####ANCESTRAL CHARACTER RECONSTRUCTION####
# Perform the ancestral reconstruction of rostral and mimicry with fitDiscrate. Here you can changes between rostral or mimicry. 
# In $data the next numbers is the column in or traits table, would depend on the location of or trait##
# When data are binary, just use ER and ARD model, as the case of rostral, in the case of mimicry use the three models
xen_geiger <- treedata(Tree, X)
xen_geiger$phy
xen_geiger$data

xengeigerer <- fitDiscrete(xenodon, xen_geiger$data [,1], type = "discrete", model = "ER")
xengeigerer

xengeigersym <- fitDiscrete(xenodon, xen_geiger$data [,1], type = "discrete", model = "SYM")
xengeigersym

xengeigerard <- fitDiscrete(xenodon, xen_geiger$data [,1], type = "discrete", model = "ARD")
xengeigerard

# Model Selection
modSel.geiger(xengeigerer, xengeigersym, xengeigerard, type="AICc")

ER_SYM <- abs(xengeigerer$opt$aicc - xengeigersym$opt$aicc)
ER_SYM

ER_ARD <- abs(xengeigerer$opt$aicc - xengeigerard$opt$aicc)
ER_ARD

SYM_ARD <- abs(xengeigersym$opt$aicc - xengeigerard$opt$aicc)
SYM_ARD


dER_ard <- abs(2*(xengeigerer$opt$lnL-xengeigerard$opt$lnL))
dER_ard

pvalueER_ard <- pchisq(dER_ard, 2-1, lower.tail = FALSE)
pvalueER_ard

### MCMC stochastic character mapping. To perform the mimicry mapping, changes rostral to mimicry, following the original table

cols<- c("black", "#DF536B")
mtrees<-make.simmap(Tree,rostral,model="ER")
mtrees
plot(mtrees)
plot(mtrees,type="phylogram",fsize=0.8,ftype="i")

add.simmap.legend(colors=cols2,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(xenodon)),fsize=0.8)

mtrees2<-make.simmap(xenodon,rostral,model="ER",nsim=100)

par(mfrow=c(10,10))
null<-sapply(mtrees2,plot,lwd=1,ftype="off")
pd<-summary(mtrees2,plot=FALSE)
pd
plot(pd,fsize=0.6,ftype="i")

plot(mtrees2[[1]],type="phylogram",fsize=0.8,ftype="i", lwd=4)

nodelabels(pie=pd$ace,piecol=cols,cex=0.5)
add.simmap.legend(colors=cols2,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(xenodon)),fsize=0.8)
axisPhylo()

# Other options of representation

plot(fitER$lik.anc,pd$ace,xlab="marginal ancestral states",
     ylab="posterior probabilities from stochastic mapping")
lines(c(0,1),c(0,1),lty="dashed",col="red",lwd=2)

# Threshold model

# Testing correlation between rostral and habitat
matrix <- read.csv2("binary_matrix.csv", row.names = 1) #load Binary Matrix
hab.ros<-matrix[,c(1,2)]
hab.ros<-as.matrix(hab.ros)
da<- match.phylo.data(xenodon, hab.ros)
da$data

# set some parameters for the MCMC
sample <- 100
ngen <- 1000000
burnin <- 0.2 * ngen


# Correlation between rostral and habitat, This procedure would take at least 8 hours, depending on your computer specifications,
# and the number of generations ngen
hab.ros <- threshBayes(xenodon, da$data, types = c("disc", "disc"),
                       ngen = 1000000,control=list(sample=sample))

plot(hab_ros2)
d <- density(hab_ros2)
dev.off()
print(hab_ros2)
plot(d)
title(main=paste("posterior density of correlation",
                 "coefficient, r,\nfrom threshold model"),
      font.main=3)
# what is our estimate of the correlation?
mean(hab_ros2$par[(burnin/sample + 1):nrow(hab_ros2$par), "r"])

# plot our likelihood profile
plot(hab_ros2$par[, "logL"], type = "l", xlab = "generation", 
     ylab = "logL")

# plot our posterior sample for the correlation
h<-hist(hab_ros2$par[(burnin/sample + 1):nrow(hab_ros2$par), "r"],
        xlab = "r", ylab = "frequency",  main = NULL)
lines(c(0.03,0.03),c(0,max(h$density)),lty="dashed")
plot(density(hab_ros2$par[(burnin/sample+1):nrow(hab_ros2$par),
                          "r"],bw=0.1),xlab="r",main="posterior density for r")
lines(c(r,r),c(0,1000),lty="dashed")


##################################333
# Phylogenetic signal  for discrate characters####
#This script follows this workd:
#Rui Borges, João Paulo Machado, Cidália Gomes, Ana Paula Rocha, Agostinho Antunes; 
#Measuring phylogenetic signal between categorical traits and phylogenies, 
#Bioinformatics, https://doi.org/10.1093/bioinformatics/bty800
#More information can be found here https://github.com/mrborges23/delta_statistic



setwd("C:/Users/Hugodryas/Documents/Documentos/Doctorado/Xenodon/Xenodon nostril evolution/R/Deltasignal")
library(ape)
source("code.R")

dir()

#Load the Tree and traits

#trait
X<-read.csv("",row.names=1, sep=";")
rostral<-setNames(X[,1],rownames(X))
mimicry<-setNames(X[,2], rownames(X))
habitat<-setNames(X[,4], rownames(X))

#tree
setwd("")

Tree<-read.tree("")
Tree
head(xenodon_out)
plot(Tree, lwd=4)
axisPhylo()

#It is important to guarantee that all the branches are positive as this method requires a 
#metric-tree (i.e., branch_lengths > 0). Here, we take 1% of the 1% quantile to fill 
#in the null branches:

xenodon$edge.length[xenodon$edge.length==0] <- quantile(xenodon$edge.length,0.1)*0.1

#confirm that the trait order follow the species order in the tree
xenodon$tip.label

#calculate delta

deltaRostral<- delta(rostral,xenodon,0.1,0.0589,10000,10,100) 

#calculated p-value, compute the probability p(random_delta>deltaA) in the null distribution, which returns the p-value.

random_delta_rostral <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(rostral)
  random_delta_rostral[i] <- delta(rtrait,xenodon,0.1,0.0589,10000,10,100)
}
p_value <- sum(random_delta_rostral>deltaRostral)/length(random_delta_rostral)
boxplot(random_delta_rostral)
abline(h=deltaRostral,col="red")

#The same for mimic

deltaMimic<- delta(mimic,xenodon,0.1,0.0589,10000,10,100) 

#calculated p-value, compute the probability p(random_delta>deltaA) in the null distribution, which returns the p-value.

random_delta_mimic <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(mimic)
  random_delta_mimic[i] <- delta(rtrait,xenodon,0.1,0.0589,10000,10,100)
}
p_value_mimic <- sum(random_delta_mimic>deltaMimic)/length(random_delta_mimic)
boxplot(random_delta_mimic)
abline(h=deltaMimic,col="red")

