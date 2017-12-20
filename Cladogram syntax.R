
#                            Creating Cladograms in R

###############################################################################################

#Install/Load Packages

#ape is a program used in r that has functions for manipulating phylogenetic trees
install.packages("ape",dependencies = TRUE)
library(ape)

#phangorn is a program in r that uses estimation methods for phylogenetic tree building
install.packages("phangorn",dependencies = TRUE)
library(phangorn)


#seqinr is a program in r that helps read in DNA and protien data files 
install.packages("seqinr",dependencies = TRUE)
library(seqinr)

#Natalie's section
################################################################################################

#Import sequence data 

smcc<-read.dna("convertedraw(2).phy",format = "interleaved")
head(smcc)

#convert the alignment to a "phyDat" object in order to use treebuilding operations in phangorn.
smc_phyDat <- phyDat(smcc, type = "DNA", levels = NULL)

#Milo's section
################################################################################################

#Model Testing and Distance Matrice

#using the phangorn program to create a distrance matrix of our data.
#the distance matrix turns sequence alignments into a matrix of pairwise distances. 
#This is the first step for the phangorn program to calculate distance between species.
mt<-modelTest(smc_phyDat)
print(mt)

#this syntax helps compute pairwise distances. the arguments under model can be either "JC69" or "F81"
#the "JC69" assumes equal base frequences between reads, the "F81" model does not. 
dna_dist<-dist.ml(smc_phyDat,model = "JC69")

#Natalie's section
###############################################################################################

#Neighbor joining, UPGMA, and Maximim Parsimony

#UPGMA stands for Unweighted Pair Group Method with Arithmetic Mean it is a agglomerative hierarchical clustering method. 
#we use the phangorn package to run this test, it is built into the package
smc_UPGMA <- upgma(dna_dist)
plot(smc_UPGMA, main="UPGMA")

#Neighbor joining is a tree building method computing the lengths of the branches of this tree.
#we use the phangorn package to run this test, it is already built into the package
smc_NJ<- NJ(dna_dist)
plot(smc_NJ, main="Neighbor Joining")

#Both
###############################################################################################

#Testing which tree is the best 

#the parsimony function returns a parsimony score
#parsimony gives a score to the tree using either the sankoff or fitch algorithm.

# Fitch's algorithm determines the set of possible states for each node using a bottom- up from leaves to root, it also looks from a top- down point of view where it picks the ancestral state for each node from the set of possibilities.
parsimony(smc_UPGMA,smc_phyDat,method="fitch")
parsimony(smc_NJ,smc_phyDat,method="fitch")

#Sankoff's algorithem relies on counting the smallest number of possible (weighted) state changes needed on a given tree.
parsimony(smc_UPGMA,smc_phyDat,method="sankoff")
parsimony(smc_NJ,smc_phyDat,method="sankoff")

#based on these results the UPGMA is the better tree

#Milo's section
###############################################################################################

#Maximum Likelyhood/ Bootstrapping

#this is another method to test the relationships between species
#here we will use both the ape and phangorn package to find the maximum likelyhood of both plots

#key
#pml returns an object which contains the data, the tree, and parameters of the model
#optim.pml computes the likelihood of a phylogenetic tree given a sequence alignment and a model. the model is "JC" beacuse we are assuming equal base frequencies
#bootstrap.pml produces a list of bootstrapped data sets. 

#the following is the process to create the final plots with bootstrap confidence.

#plot for UPGMA
fitU<- pml(smc_UPGMA, smc_phyDat) 
print(fitU)

fitJCU<-optim.pml(fitU, model = "JC", rearrangement = "stochastic")
logLik(fitJCU)

bsU<-bootstrap.pml(fitJCU,bs=100, optNni=TRUE,multicore=TRUE,control=pml.control(trace=0))

plotBS(midpoint(fitJCU$tree),bs, p=50, type="p",title(main="UPGMA"))

#plot for Neighbor joining
fitN <- pml(smc_NJ, smc_phyDat)
print(fitN)

fitJCN <-optim.pml(fitN, model = "JC", rearrangement = "stochastic")
logLik(fitJCN)

bsN<-bootstrap.pml(fitJCN,bs=100, optNni=TRUE,multicore=TRUE,control=pml.control(trace=0))

plotBS(midpoint(fitJCN$tree),bs, p=50, type="p",title(main="Neighbor Joining"))

#Both
###############################################################################################
