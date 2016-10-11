## script that does LDA for diatom count data from Gull Lake and compares this with
##    CONISS as per previous analysis
## see also https://www.quora.com/Could-latent-Dirichlet-allocation-solved-by-Gibbs-sampling-versus-variational-EM-yield-different-results

library("topicmodels")
library("rioja")
library("vegan")
library("ggplot2")
library("reshape2")
library("extrafont")

## create a theme to save linespace in plots
papertheme <- theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top')

## read in data
## Gull data from email Heather sent June 20, 2016 'Gall Lake diatom raw counts'
##    They are as indicated raw counts, with number of diatoms counted per sample as a column
if (!file.exists('../data/private/7.5mCoreRawCounts.csv')) {
  stop("get Gull data from Heather or Emma")
}
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
nlist <- 2:6
reptimes <- length(nlist)
myctrl <- list(seed = 2010)

creps <- do.call("list", replicate(reptimes, dcounts, simplify = FALSE))
ctrls <- do.call("list", replicate(reptimes, myctrl, simplify = FALSE))

mods <- mapply(LDA, k=nlist, x=creps, control=ctrls)
modsgibbs <- mapply(LDA, k=nlist, x=creps, control = ctrls, method="Gibbs")

## compare them
aics <- do.call(rbind,lapply(mods, AIC))
bics <- do.call(rbind,lapply(mods, BIC))

aicsgibbs <- do.call(rbind,lapply(modsgibbs, AIC))
bicsgibbs <- do.call(rbind,lapply(modsgibbs, BIC))

## make into one data frame for ggplot
aicsall <- as.data.frame(rbind(aics, aicsgibbs))
bicsall <- as.data.frame(rbind(bics, bicsgibbs))

aicsall$IC <- "AIC"
bicsall$IC <- "BIC"
aicsall$Sampling <- rep(c("VEM", "Gibbs"), each=nrow(aicsall)/2)
bicsall$Sampling <- rep(c("VEM", "Gibbs"), each=nrow(bicsall)/2)

icsall <- rbind(aicsall, bicsall)
icsall$Ngroups <- rep(nlist)

## plot output
ICplot <- ggplot(icsall, aes(x=Ngroups, y=V1, group=IC)) +
  papertheme +
  geom_point() +
  facet_grid(Sampling ~ IC, scales='free_y') +
  guides(colour=guide_legend(nrow=1,bycol =TRUE,title.position = 'left')) +
  theme(axis.title.y=element_blank()) +
  xlab("n(groups)")
## FIXME: Gibbs just likes lots of groups...

pdf("../data/private/ICs-gull.pdf", onefile = TRUE)
ICplot
dev.off()

# select best manually
ldagull <- mods[[2]]
ldagullgibbs <- mods[[5]]

#get parameter estimates
z=posterior(ldagull)
commun.plot=z$topics
commun.spp=z$terms

zgibbs=posterior(ldagullgibbs)
commun.plotgibbs=zgibbs$topics
commun.sppgibbs=zgibbs$terms

#plot relative abundance of component communities for each sampling unit 
#(i.e., along the gradient)
groups <- dim(commun.plot)[2]
plot(NA,NA,xlim=c(0,nrow(diats)),ylim=c(0,1),xlab='Gradient/Sampling units',ylab='Relative abundance')
for (i in 1:groups){
  lines(1:nrow(diats),commun.plot[,i],col=i)
}

groups <- dim(commun.plotgibbs)[2]
plot(NA,NA,xlim=c(0,nrow(diats)),ylim=c(0,1),xlab='Gradient/Sampling units',ylab='Relative abundance')
for (i in 1:groups){
  lines(1:nrow(diats),commun.plotgibbs[,i],col=i)
}

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
#(i.e. black line considerably above red line)

diatgroups <- cutree(clust, k = 4)

## ==================================================================================================
## SAVE OBJECTS
## ==================================================================================================

## save chosen lda model
saveRDS(ldagull, "../data/private/lda-gull.rds")
saveRDS(ldagullgibbs, "../data/private/ldagibbs-gull.rds")

## save other strati objects
saveRDS(diatgroups, "../data/private/gull-coniss.rds")
