## script that does LDA for diatom count data from Gill Lake and compares this with
##    CONISS as per previous analysis

library("topicmodels")
library("rioja")
library("vegan")
library("ggplot2")
library("reshape")

papertheme <- theme_bw(base_size=12, base_family = 'Arial') +
  theme(legend.position='top')

## read in data
## Gill data from email Heather sent June 20, 2016 'Gall Lake diatom raw counts'
##    They are as indicated raw counts, with number of diatoms counted per sample as a column
alldiats <- read.csv('../data/private/7.5mCoreRawCounts.csv') 

## remove things that aren't diatoms and change NA to 0
alldiats[is.na(alldiats)] <- 0
notdiats <- which(names(alldiats) %in% c('bottom', 'top', 'AGE', 'DI.Depth',
                                         'Unknowns', 'X..Diatoms', 'Chrysophyte.scales',
                                         'Chysophyte.cysts', 'C.D.index', 'Microspheres'))

if(any(rowSums(alldiats[,-notdiats], na.rm = TRUE) == 0)) {
  discardrows <- which(rowSums(alldiats[,-notdiats]) == 0)
  alldiats <- alldiats[-discardrows, ]
}
if(any(colSums(alldiats, na.rm = TRUE) == 0)) {
  discardcols <- which(colSums(alldiats, na.rm = TRUE) == 0)
  alldiats <- alldiats[, -discardcols]
}
notdiats <- which(names(alldiats) %in% c('bottom', 'top', 'AGE', 'DI.Depth',
                                         'Unknowns', 'X..Diatoms', 'Chrysophyte.scales',
                                         'Chysophyte.cysts', 'C.D.index', 'Microspheres'))

diats <- alldiats[,-notdiats]

dcounts <- as.matrix(diats)
## FIXME: what about rare species? prune for e.g. 5%? 

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

## make into ggplot
agedf <- data.frame(cbind(ldadiat@gamma, topics(ldadiat), alldiats$AGE))
colnames <- c(as.character(seq_len(ldadiat@k)), "Cluster", "Year")
names(agedf) <- colnames
agedf$Cluster <- factor(agedf$Cluster)

agestack <- melt(agedf, id.vars=c('Year','Cluster'), variable_name = 'Group')
agesplit <- with(agestack, split(agestack, list(Group)))
#lapply(agesplit, summary)

ggplot(data = agestack, aes(x=Year, y=value, lty=Group)) +
  papertheme +
  scale_linetype_manual(name='Group', values = c("solid", "longdash","dotdash")) +
  scale_colour_manual(name="Cluster", values = c("#5e3c99", "#e66101","#b2abd2"))+
  geom_line() +
  geom_point(aes(col=Cluster)) +
  theme(legend.box = "horizontal") +
  ylab("Proportion of population")



## ==================================================================================================
## CONISS
## ==================================================================================================
##  Diatom zones in the composite core were determined by a depth-constrained cluster analysis, 
##      using Euclidean distance as a measure of similarity, conducted with Tilia v. 2.02 (Grimm 1987), 
##      and by optimal zonation using the broken-stick method in PSIMPOLL v. 4.10 (Bennett 1996).
## from heather: "The broken stick was just something we ran and i dont think ever used the results... 
##    it was similar enough to coniss that we just left the coniss"
diss <- vegdist(diats, method='euclidean')
clust <- chclust(diss, method = 'coniss')
bclust <- bstick(clust) # --> 4 zones as drawn in publication

diats$Group <- cutree(clust, k = 4)

## ==================================================================================================
## SAVE OBJECTS
## ==================================================================================================

## save chosen lda model
saveRDS(ldaclado, "../data/lda-clado.rds")
