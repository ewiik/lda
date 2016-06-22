## script that does LDA for diatom count data from Gill Lake

library("topicmodels")

## read in data
## Gill data from email Heather sent June 20, 2016 'Gall Lake diatom raw counts'
##    They are as indicated raw counts, with number of diatoms counted per sample as a column
alldiats <- read.csv('../data/private/7.5mCoreRawCounts.csv') 

## remove things that aren't diatoms
diats <- alldiats[,-which(names(alldiats) %in% c('bottom', 'top', 'AGE', 'DI.Depth',
                                                 'Unknowns', 'X..Diatoms', 'Chrysophyte.scales',
                                                 'Chysophyte.cysts', 'C.D.index', 'Microspheres'))]
diats[is.na(diats)] <- 0

discardrows <- which(rowSums(diats, na.rm = TRUE) == 0)
discardcols <- which(colSums(diats, na.rm = TRUE) == 0)
diats <- diats[-discardrows, -discardcols]

dcounts <- as.matrix(diats)
## FIXME: what about rare species? prune for e.g. 5%? and these are all just counts up to ca 400
## individuals, not normalised to an amount of sediment for example

#run LDA for multiple settings of group number
SEED <- 2010
nlist <- 2:6
reptimes <- length(nlist)
creps <- do.call("list", replicate(reptimes, dcounts, simplify = FALSE))
mods <- mapply(LDA, k=nlist, x=creps)

# compare them
aics <- do.call(rbind,lapply(mods, AIC))
plot(aics ~ nlist, ylab = 'AIC', xlab='n(groups)')

# select best manually
ldadiat <- mods[[2]]

#get parameter estimates
z=posterior(ldadiat)
commun.plot=z$topics
commun.spp=z$terms

#plot relative abundance of component communities for each sampling unit 
#(i.e., along the gradient)
groups <- dim(commun.plot)[2]
plot(NA,NA,xlim=c(0,nrow(diats)),ylim=c(0,1),xlab='Gradient/Sampling units',ylab='Relative abundance')
for (i in 1:groups){
  lines(1:nrow(diats),commun.plot[,i],col=i)
}

#plot relative abundance of species 1:n in component community
species <- dim(commun.spp)[2]
opar <- par()
par(mfrow=c(groups,1))
for (i in 1:groups){
  plot(1:species,commun.spp[i,],type='l',col=i,ylim=c(0,0.5),xlab='Species',ylab='Relative abundance')
}
par(opar)

## save chosen model
saveRDS(ldaclado, "../data/lda-clado.rds")
