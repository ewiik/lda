## script that does LDA for count data from Cunswick Tarn, etc. and compares it with previously
##    done analysis (taken code from cunspaper.pcurves.R and followed wiiketal)
## see also https://www.quora.com/Could-latent-Dirichlet-allocation-solved-by-Gibbs-sampling-versus-variational-EM-yield-different-results

library("topicmodels")
library("analogue")
library("ggplot2")
library("reshape")

## read in data
## Cunswick data from my files in static R folder (the 'archive'); multiplied to individuals
##    per 100ml of sediment
if (any(!file.exists("../data/cunscladofinalless.csv",
                     "../data/cladoyears.rds"))) {
  stop("Get cladoceran data from Emma or from open data repo TBC") }
allclados <- read.csv('../data/cunscladofinalless.csv') 
cladoyears <- readRDS('../data/cladoyears.rds') 

## remove things that aren't cladocera or are ambiguous
clados <- allclados[,-which(names(allclados) %in% c('chydEPH', 'chaojaw'))]
clados$Year <- cladoyears
ccounts <- clados[,-which(names(clados) %in% c('Year', 'depth'))]
ccounts <- as.matrix(ccounts)
ccounts <- apply(ccounts, 2, round)

#run LDA for multiple settings of group number
SEED <- 2010
nlist <- 2:9
reptimes <- length(nlist)
creps <- do.call("list", replicate(reptimes, ccounts, simplify = FALSE))
mods <- mapply(LDA, k=nlist, x=creps)

# compare them
aics <- do.call(rbind,lapply(mods, AIC))
bics <- do.call(rbind,lapply(mods, BIC))
pdf("../data/private/ICs-clado.pdf", onefile = TRUE)
plot(aics ~ nlist, ylab = 'AIC', xlab='n(groups)')
plot(bics ~ nlist, ylab = 'BIC', xlab='n(groups)')
dev.off()

# select best manually
ldaclado <- mods[[5]]

#get parameter estimates
z=posterior(ldaclado)
commun.plot=z$topics
commun.spp=z$terms

#plot relative abundance of component communities for each sampling unit 
#(i.e., along the gradient)
groups <- dim(commun.plot)[2]
plot(NA,NA,xlim=c(0,33),ylim=c(0,1),xlab='Gradient/Sampling units',ylab='Relative abundance')
for (i in 1:groups){
  lines(1:32,commun.plot[,i],col=i)
}

#plot relative abundance of species 1:n in component community
opar <- par()
par(mfrow=c(groups,1))
for (i in 1:groups){
  plot(1:22,commun.spp[i,],type='l',col=i,ylim=c(0,0.5),xlab='Species',ylab='Relative abundance')
}
par(opar)


## ==================================================================================================
## PRINCIPAL CURVES AND MVPART RESULTS...
## ==================================================================================================
## make data relative
cladorel <- clados[,-which(names(clados) %in% c('Year', 'depth'))]
cladorel <- decostand(cladorel, method = 'total')

cladocurve <- prcurve(sqrt(cladorel), method = "pca", smoother = smoothSpline, vary = TRUE, 
                      finalCV = FALSE, axis = 1, thresh = 0.001, penalty = 1.4,plotit = TRUE) #
plot(cladocurve$lambda  ~ cladoyears, type = "o")

## mvpart not available but these were the four boundaries established (note that it depended
##    on whether or not relative or not data) : 1427, 1898, 1940, 1990
mvpartsplits <- c(1427, 1898, 1940, 1990)

saveRDS(cladocurve, "../data/cuns-pcurve.rds")
saveRDS(mvpartsplits, "../data/cuns-mvpartsplits.rds")

## ==================================================================================================
## SAVE OBJECTS
## ==================================================================================================

## save chosen model
saveRDS(ldaclado, "../data/lda-cuns.rds")
