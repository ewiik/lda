## script that does LDA for diatom count data from Foy Lake and compares this with
##    previous analysis
## see also https://www.quora.com/Could-latent-Dirichlet-allocation-solved-by-Gibbs-sampling-versus-variational-EM-yield-different-results

library("topicmodels")
library("rioja")
library("vegan")
library("ggplot2")
library("reshape")

## read in data
## openly available with publication Spanbaueretal2014ProlongedIinstabilityBeforeRegimeShift
if (!file.exists('../data/Spanbaueretal2014supp.csv')) {
  stop("get supporting data from Spanbaueretal2014 paper")
}
foy <- read.csv('../data/Spanbaueretal2014supp.csv') 
foyyears <- foy$YB1950

## remove things that aren't diatoms
foy <- foy[,-which(names(foy) %in% c('YB1950', 'Sample'))]

## make % data into 'count' data
foy <- foy*100
foycounts <- as.matrix(foy)
foycounts <- apply(foycounts, 2, round)


#run LDA for multiple settings of group number.. Gibbs vs VEM?
SEED <- 2010
nlist <- 3:8
reptimes <- length(nlist)
creps <- do.call("list", replicate(reptimes, foycounts, simplify = FALSE))

mods <- mapply(LDA, k=nlist, x=creps)
modsgibbs <- mapply(LDA, k=nlist, x=creps, method="Gibbs")

# compare them
aics <- do.call(rbind,lapply(mods, AIC))
bics <- do.call(rbind,lapply(mods, BIC))
#aicsgibbs <- do.call(rbind,lapply(modsgibbs, AIC))

pdf("../data/private/ICs-foy.pdf", onefile = TRUE)
plot(aics ~ nlist, ylab = 'AIC', xlab='n(groups)')
plot(bics ~ nlist, ylab = 'BIC', xlab='n(groups)')
dev.off()

#plot(aics ~ nlist, ylab = 'AIC', xlab='n(groups)', ylim=c(130000,300000))
#points(aicsgibbs ~ nlist, pch = 3)

## select best manually.... though gibbs is different
## NOTE.... Is the rogue community (6) really necessary out of mods? 
ldafoy <- mods[[3]]
#ldafoygibbs <- modsgibbs[[]]

#get parameter estimates
z=posterior(ldafoy)
commun.plot=z$topics
commun.spp=z$terms

#plot relative abundance of component communities for each sampling unit 
#(i.e., along the gradient)
groups <- dim(commun.plot)[2]
plot(NA,NA,xlim=c(0,nrow(foy)),ylim=c(0,1),xlab='Gradient/Sampling units',ylab='Relative abundance')
for (i in 1:groups){
  lines(1:nrow(foy),commun.plot[,i],col=i)
}

#plot relative abundance of species 1:n in component community
species <- dim(commun.spp)[2]
opar <- par()
par(mfrow=c(groups,1))
for (i in 1:groups){
  plot(1:species,commun.spp[i,],type='l',col=i,ylim=c(0,0.5),xlab='Species',ylab='Relative abundance')
}
par(opar)



## ==================================================================================================
## CONISS
## ==================================================================================================
##  Improvising this since no clustering was done before for Foy data
diss <- vegdist(foy, method='euclidean')
clust <- chclust(diss, method = 'coniss')
bclust <- bstick(clust) # --> seems like 6

foyclust <- foy
foyclust$Group <- cutree(clust, k = 6)

## ==================================================================================================
## SAVE OBJECTS
## ==================================================================================================

## save chosen lda model
saveRDS(ldafoy, "../data/lda-foy.rds")
