## script that does LDA for count data from Cunswick Tarn, etc. and compares it with previously
##    done analysis (taken code from cunspaper.pcurves.R and followed wiiketal)

library("topicmodels")
library("analogue")
library("ggplot2")

papertheme <- theme_bw(base_size=12, base_family = 'Arial') +
  theme(legend.position='top')

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

## make into ggplot
agedf <- data.frame(cbind(ldaclado@gamma, topics(ldaclado), cladoyears))
colnames <- c(as.character(seq_len(ldaclado@k)), "Group", "Year")
names(agedf) <- colnames
agedf$Group <- factor(agedf$Group)

agestack <- melt(agedf, id.vars=c('Year','Group'), variable_name = 'Topic')
agesplit <- with(agestack, split(agestack, list(Group)))
lapply(agesplit, summary)

ggplot(data = agestack, aes(x=Year, y=value, col=Topic, lty=Topic)) +
  papertheme +
  scale_linetype_manual(name='Group', values = c("dotdash", "solid","solid", "longdash", "dotdash", 
                                                            "dotdash")) +
  scale_colour_manual(name="Group", values = c("#b2abd2", "#e66101","#b2abd2", "#5e3c99", "#e66101",
                                              "#5e3c99"))+
  geom_vline(xintercept = c(1390, 1900, 1930, 2000), col='grey50') +
  geom_line() 
  


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


## ==================================================================================================
## SAVE OBJECTS
## ==================================================================================================

## save chosen model
saveRDS(ldaclado, "../data/lda-clado.rds")
