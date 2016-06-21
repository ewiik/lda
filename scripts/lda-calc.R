## script that does LDA for count data from Cunswick Tarn, etc.

library("topicmodels")

## read in data
## Cunswick data from my files in static R folder (the 'archive'); multiplied to individuals
##    per 100ml of sediment
allclados <- read.csv('~/git/lda/data/private/cunscladofinalless.csv') 
cladoyears <- readRDS('~/git/lda/data/private/cladoyears.rds') 

## remove things that aren't cladocera or are ambiguous
names(clados)
clados <- allclados[,-which(names(clados) %in% c('chydEPH', 'chaojaw'))]
clados$Year <- cladoyears
ccounts <- clados[,-which(names(clados) %in% c('Year', 'depth'))]
ccounts <- as.matrix(ccounts)
ccounts <- apply(ccounts, 2, round)

#run LDA
SEED=2010
VEM=LDA(ccounts,k=5, control = list(seed = SEED))

#get parameter estimates
z=posterior(VEM)
commun.plot=z$topics
commun.spp=z$terms

#plot relative abundance of component communities for each sampling unit 
#(i.e., along the gradient)
plot(NA,NA,xlim=c(0,33),ylim=c(0,1),xlab='Gradient/Sampling units',ylab='Relative abundance')
for (i in 1:5){
  lines(1:32,commun.plot[,i],col=i)
}

#plot relative abundance of species 1,...,2000 in component community
opar <- par()
par(mfrow=c(5,1))
for (i in 1:5){
  plot(1:22,commun.spp[i,],type='l',col=i,ylim=c(0,0.5),xlab='Species',ylab='Relative abundance')
}
par(opar)
