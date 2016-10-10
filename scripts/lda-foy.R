## script that does LDA for diatom count data from Foy Lake and compares this with
##    previous analysis
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
## openly available with publication Spanbaueretal2014ProlongedIinstabilityBeforeRegimeShift
if (!file.exists('../data/Spanbaueretal2014supp.csv')) {
  stop("get supporting data from Spanbaueretal2014 paper")
}
foy <- read.csv('../data/Spanbaueretal2014supp.csv') 
foyyears <- foy$YB1950
## and here, the less openly available raw counts, which also span a wider time frame
##    Provided by Gavin 4th Oct 2016 (see emails)
if (!file.exists('../data/private/Foy-RawDiatomCounts.csv')) {
  stop("get raw count data from Emma or Gavin")
}
foyraw <- read.csv("../data/private/Foy-RawDiatomCounts.csv", na.strings = c('x','\'','\`'))
foyrawyears <- foyraw$YB1950

## remove things that aren't diatoms
foy <- foy[,-which(names(foy) %in% c('YB1950', 'Sample'))]

foyraw <- foyraw[,-which(names(foyraw) %in% c('YB1950', 'Sample', 'Sample..'))]
# all names after Total Diatoms are not-diatoms....
nomore <- grep("Total", names(foyraw))
foyraw <- foyraw[,-c(nomore:length(names(foyraw)))]

## make % data into 'count' data
foy <- foy*100
foycounts <- as.matrix(foy)
foycounts <- apply(foycounts, 2, round)

## remove aggregate designations (but save separately just in case)
## Note that this still keeps genera-level aggregates like Cyclotella sp.
generics <- c(grep('oid', names(foyraw)), grep('ales', names(foyraw)), grep('Dissolved', names(foyraw)))
foyaggs <- foyraw[,generics]
foyraw <- foyraw[,-generics]

## all 0 are NA in raw...; also 'Each row of the input matrix needs to contain 
##    at least one non-zero entry'
## FIXME: why are there some rows with no entries?
foyraw[is.na(foyraw)] <- 0
losers <- which(rowSums(foyraw) == 0)
foyraw <- foyraw[-losers,]
foyrawyears <- foyrawyears[-losers]

#run LDA for multiple settings of group number, and for VEM and Gibbs
nlist <- 3:8
reptimes <- length(nlist)
myctrl <- list(seed = 2010)

creps <- do.call("list", replicate(reptimes, foycounts, simplify = FALSE))
crepsraw <- do.call("list", replicate(reptimes, foyraw, simplify = FALSE))
ctrls <- do.call("list", replicate(reptimes, myctrl, simplify = FALSE))

mods <- mapply(LDA, k=nlist, x=creps, control = ctrls)
modsgibbs <- mapply(LDA, k=nlist, x=creps, control = ctrls, method="Gibbs")

modsraw <- mapply(LDA, k=nlist, x=crepsraw, control = ctrls)
modsgibbsraw <- mapply(LDA, k=nlist, x=crepsraw, control = ctrls, method="Gibbs")

## compare them with AIC and BIC
aics <- do.call(rbind,lapply(mods, AIC))
bics <- do.call(rbind,lapply(mods, BIC))

aicsraw <- do.call(rbind,lapply(modsraw, AIC))
bicsraw <- do.call(rbind,lapply(modsraw, BIC))

aicsgibbs <- do.call(rbind,lapply(modsgibbs, AIC))
bicsgibbs <- do.call(rbind,lapply(modsgibbs, BIC))

aicsgibbsraw <- do.call(rbind,lapply(modsgibbsraw, AIC))
bicsgibbsraw <- do.call(rbind,lapply(modsgibbsraw, BIC))

## make into one data frame for ggplot
aicsall <- as.data.frame(rbind(aics, aicsraw, aicsgibbs, aicsgibbsraw))
bicsall <- as.data.frame(rbind(bics, bicsraw, bicsgibbs, bicsgibbsraw))

aicsall$IC <- "AIC"
bicsall$IC <- "BIC"
aicsall$Sampling <- rep(c("VEM", "Gibbs"), each=nrow(aicsall)/2)
aicsall$Data <- rep(c("Published", "Raw"), each=nrow(aicsall)/4, times=2)

bicsall$Sampling <- rep(c("VEM", "Gibbs"), each=nrow(bicsall)/2)
bicsall$Data <- rep(c("Published", "Raw"), each=nrow(bicsall)/4, times=2)

icsall <- rbind(aicsall, bicsall)
icsall$Ngroups <- rep(nlist)
icsall$Ngroups[icsall$IC == "AIC"] <- icsall$Ngroups[icsall$IC == "AIC"] + 0.2

## plot all ICs
ICplot <- ggplot(icsall, aes(x=Ngroups, y=V1, col=IC, group=Data)) +
  papertheme +
  geom_point() +
  scale_color_brewer(type="qual", palette = 'Set1') +
  facet_grid(Data ~ Sampling, scales='free_y') +
  guides(colour=guide_legend(nrow=1,bycol =TRUE,title.position = 'left')) +
  theme(axis.title.y=element_blank()) +
  xlab("n(groups) horizontally jittered")


pdf("../data/ICs-foy.pdf", onefile = TRUE)
ICplot
dev.off()

## select best manually.... though gibbs is different
## NOTE.... Is the rogue community (6) really necessary out of mods? 
ldafoy <- modsraw[[4]]
ldafoygibbs <- modsgibbsraw[[6]]

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
plot(NA,NA,xlim=c(0,nrow(foyraw)),ylim=c(0,1),xlab='Gradient/Sampling units',ylab='Relative abundance')
for (i in 1:groups){
  lines(1:nrow(foyraw),commun.plot[,i],col=i)
}
groups <- dim(commun.plotgibbs)[2]
plot(NA,NA,xlim=c(0,nrow(foyraw)),ylim=c(0,1),xlab='Gradient/Sampling units',ylab='Relative abundance')
for (i in 1:groups){
  lines(1:nrow(foyraw),commun.plotgibbs[,i],col=i)
}


## ==================================================================================================
## CONISS
## ==================================================================================================
##  Improvising this since no clustering was done before for Foy data
diss <- vegdist(foyraw, method='euclidean')
clust <- chclust(diss, method = 'coniss')
bclust <- bstick(clust) # --> seems like 6

foygroups <- cutree(clust, k = 6)

## ==================================================================================================
## SAVE OBJECTS
## ==================================================================================================

## save chosen lda model
saveRDS(ldafoy, "../data/lda-foy.rds")
saveRDS(ldafoygibbs, "../data/ldagibbs-foy.rds")

## save clustering objects
saveRDS(foygroups, "../data/foy-coniss.rds")
saveRDS(foyrawyears, "../data/foyyears.rds")

