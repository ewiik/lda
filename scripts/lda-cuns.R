## script that does LDA for count data from Cunswick Tarn, etc. and compares it with previously
##    done analysis (taken code from cunspaper.pcurves.R and followed wiiketal)
## see also https://www.quora.com/Could-latent-Dirichlet-allocation-solved-by-Gibbs-sampling-versus-variational-EM-yield-different-results

library("topicmodels")
library("analogue")
library("ggplot2")
library("reshape")
library("extrafont")

## create a theme to save linespace in plots
papertheme <- theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top')

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
nlist <- 2:9
reptimes <- length(nlist)
myctrl <- list(seed = 2010)

creps <- do.call("list", replicate(reptimes, ccounts, simplify = FALSE))
ctrls <- do.call("list", replicate(reptimes, myctrl, simplify = FALSE))

mods <- mapply(LDA, k=nlist, x=creps, control = ctrls)
modsgibbs <- mapply(LDA, k=nlist, x=creps, control = ctrls, method="Gibbs")

# compare them
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

pdf("../data/private/ICs-cuns.pdf", onefile = TRUE)
ICplot
dev.off()

# select best manually
ldacuns <- mods[[5]]
ldacunsgibbs <- modsgibbs[[8]]
## FIXME: need to increase n for Gibbs seeing as monotonic large decrease for ICs 
##    up to 9 groups!!

#get parameter estimates
z=posterior(ldacuns)
commun.plot=z$topics
commun.spp=z$terms

zgibbs=posterior(ldacunsgibbs)
commun.plotgibbs=zgibbs$topics
commun.sppgibbs=zgibbs$terms

#plot relative abundance of component communities for each sampling unit 
#(i.e., along the gradient)
groups <- dim(commun.plot)[2]
spec <- dim(commun.plot)[1]
plot(NA,NA,xlim=c(0,spec+1),ylim=c(0,1),xlab='Gradient/Sampling units',ylab='Relative abundance')
for (i in 1:groups){
  lines(1:spec,commun.plot[,i],col=i)
}

groups <- dim(commun.plotgibbs)[2]
specs <- dim(commun.plotgibbs)[1]
plot(NA,NA,xlim=c(0,spec+1),ylim=c(0,1),xlab='Gradient/Sampling units',ylab='Relative abundance')
for (i in 1:groups){
  lines(1:spec,commun.plotgibbs[,i],col=i)
}


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

## ==================================================================================================
## SAVE OBJECTS
## ==================================================================================================

## save chosen model
saveRDS(ldacuns, "../data/lda-cuns.rds")
saveRDS(ldacunsgibbs, "../data/ldagibbs-cuns.rds")

## save old cluster data
saveRDS(cladocurve, "../data/cuns-pcurve.rds")
saveRDS(mvpartsplits, "../data/cuns-mvpartsplits.rds")
