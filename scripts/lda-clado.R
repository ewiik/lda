## script that does LDA for count data from Cunswick Tarn, etc.

library("topicmodels")

## read in data
## Cunswick data from my files in static R folder (the 'archive'); multiplied to individuals
##    per 100ml of sediment
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
plot(aics ~ nlist, ylab = 'AIC', xlab='n(groups)')

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

## save chosen model
saveRDS(ldaclado, "../data/lda-clado.rds")
