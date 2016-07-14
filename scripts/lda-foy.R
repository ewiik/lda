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


#run LDA for multiple settings of group number, and for VEM and Gibbs
nlist <- 3:8
reptimes <- length(nlist)
myctrl <- list(seed = 2010)

creps <- do.call("list", replicate(reptimes, foycounts, simplify = FALSE))
ctrls <- do.call("list", replicate(reptimes, myctrl, simplify = FALSE))

mods <- mapply(LDA, k=nlist, x=creps, control = ctrls)
modsgibbs <- mapply(LDA, k=nlist, x=creps, control = ctrls, method="Gibbs")

# compare them with AIC and BIC
aics <- do.call(rbind,lapply(mods, AIC))
bics <- do.call(rbind,lapply(mods, BIC))

aicsgibbs <- do.call(rbind,lapply(modsgibbs, AIC))
bicsgibbs <- do.call(rbind,lapply(modsgibbs, BIC))

pdf("../data/ICs-foy.pdf", onefile = TRUE)
plot(aics ~ nlist, ylab = 'AIC', xlab='n(groups)')
plot(bics ~ nlist, ylab = 'BIC', xlab='n(groups)')
dev.off()

pdf("../data/ICsGibbs-foy.pdf", onefile = TRUE)
plot(aicsgibbs ~ nlist, ylab = 'AIC', xlab='n(groups)')
plot(bicsgibbs ~ nlist, ylab = 'BIC', xlab='n(groups)')
dev.off()

## select best manually.... though gibbs is different
## NOTE.... Is the rogue community (6) really necessary out of mods? 
ldafoy <- mods[[3]]
ldafoygibbs <- modsgibbs[[5]]

#get parameter estimates
z=posterior(ldafoy)
commun.plot=z$topics
commun.spp=z$terms

zgibbs=posterior(ldafoygibbs)
commun.plotgibbs=zgibbs$topics
commun.sppgibbs=zgibbs$terms

#plot relative abundance of component communities for each sampling unit 
#(i.e., along the gradient)
groups <- dim(commun.plot)[2]
plot(NA,NA,xlim=c(0,nrow(foy)),ylim=c(0,1),xlab='Gradient/Sampling units',ylab='Relative abundance')
for (i in 1:groups){
  lines(1:nrow(foy),commun.plot[,i],col=i)
}
groups <- dim(commun.plotgibbs)[2]
plot(NA,NA,xlim=c(0,nrow(foy)),ylim=c(0,1),xlab='Gradient/Sampling units',ylab='Relative abundance')
for (i in 1:groups){
  lines(1:nrow(foy),commun.plotgibbs[,i],col=i)
}


## ==================================================================================================
## CONISS
## ==================================================================================================
##  Improvising this since no clustering was done before for Foy data
diss <- vegdist(foy, method='euclidean')
clust <- chclust(diss, method = 'coniss')
bclust <- bstick(clust) # --> seems like 6

foyclust <- foy
foygroups <- cutree(clust, k = 6)

## ==================================================================================================
## SAVE OBJECTS
## ==================================================================================================

## save chosen lda model
saveRDS(ldafoy, "../data/lda-foy.rds")
saveRDS(ldafoygibbs, "../data/ldagibbs-foy.rds")

## save clustering objects
saveRDS(foygroups, "../data/foy-coniss.rds")
saveRDS(foyyears, "../data/foyyears.rds")

