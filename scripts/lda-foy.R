## script that does LDA for diatom count data from Foy Lake and compares this with
##    previous analysis
## see also https://www.quora.com/Could-latent-Dirichlet-allocation-solved-by-Gibbs-sampling-versus-variational-EM-yield-different-results

library("topicmodels")
library("rioja")
library("vegan")
library("ggplot2")
library("reshape")

papertheme <- theme_bw(base_size=12, base_family = 'Arial') +
  theme(legend.position='top')

## read in data
## openly available with publication Spanbaueretal2014ProlongedIinstabilityBeforeRegimeShift
foy <- read.csv('../data/Spanbaueretal2014supp.csv') 
foyyears <- foy$YB1950

## remove things that aren't diatoms
foy <- foy[,-which(names(foy) %in% c('YB1950', 'Sample'))]

## make % data into 'count' data
foy <- foy*100
foycounts <- as.matrix(foy)
foycounts <- apply(foycounts, 2, round)


#run LDA for multiple settings of group number.. Gibbs vs VEM
SEED <- 2010
nlist <- 3:8
reptimes <- length(nlist)
creps <- do.call("list", replicate(reptimes, foycounts, simplify = FALSE))

mods <- mapply(LDA, k=nlist, x=creps)
modsgibbs <- mapply(LDA, k=nlist, x=creps, method="Gibbs")

# compare them
aics <- do.call(rbind,lapply(mods, AIC))
aicsgibbs <- do.call(rbind,lapply(modsgibbs, AIC))

plot(aics ~ nlist, ylab = 'AIC', xlab='n(groups)', ylim=c(130000,300000))
points(aicsgibbs ~ nlist, pch = 3)

# select best manually.... 6 groups? though gibbs is different
ldafoy <- mods[[4]]
ldafoygibbs <- modsgibbs[[5]]

#get parameter estimates
z=posterior(ldafoy)
commun.plot=z$topics
commun.spp=z$terms

#plot relative abundance of component communities for each sampling unit 
#(i.e., along the gradient)
groups <- dim(commun.plot)[2]
plot(NA,NA,xlim=c(0,nrow(foy)),ylim=c(0,1),xlab='Gradient/Sampling units',ylab='Relative abundance')
for (i in 1:groups){
  lines(1:nrow(foy),commun.plot[,i],col=i)
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
agedf <- data.frame(cbind(ldafoy@gamma, topics(ldafoy), foyyears))
colnames <- c(as.character(seq_len(ldafoy@k)), "Cluster", "Year")
names(agedf) <- colnames
agedf$Cluster <- factor(agedf$Cluster)

agestack <- melt(agedf, id.vars=c('Year','Cluster'), variable_name = 'Group')
agesplit <- with(agestack, split(agestack, list(Group)))
#lapply(agesplit, summary)

ggplot(data = agestack, aes(x=Year, y=value, lty=Group)) +
  papertheme +
  scale_linetype_manual(name='Group', values = c("solid", "longdash","dotdash",
                                                 "longdash", "dotdash","solid")) +
  scale_colour_manual(name="Cluster", values = c("#5e3c99", "#e66101","#b2abd2", 
                                                 "yellow", "blue","black"))+
  geom_line() +
  geom_point(aes(col=Cluster)) +
  theme(legend.box = "horizontal") +
  ylab("Proportion of population")

foytermit <- ggplot(data = agestack, aes(y=Year, x=Group, col=factor(Cluster), size=value)) +
  papertheme +
  #scale_color_distiller(name="Value", palette = 'PuOr') +
  scale_size_continuous(name = "Relative abundance") + #guide = FALSE
  scale_color_brewer(name="Community", palette = 'Dark2') +
  #geom_abline(intercept = c(1390, 1900, 1930, 2000), col='grey50', slope = 0) +
  geom_point(alpha=0.4) +
  theme(legend.box = "vertical") +
  xlab("Species group") +
  guides(colour=guide_legend(nrow=1,byrow=TRUE))

## NOTE.... Is the rogue community (6) really necessary? what would happen if we punished it
##    down to 5 groups?

## ==================================================================================================
## CONISS
## ==================================================================================================
##  Improvising this since no clustering was done before for Foy data
diss <- vegdist(foy, method='euclidean')
clust <- chclust(diss, method = 'coniss')
bclust <- bstick(clust) # --> seems like 6

foyclust <- foy
foyclust$Group <- cutree(clust, k = 6)

## ==================================================================================================
## SAVE OBJECTS
## ==================================================================================================

## save chosen lda model
saveRDS(ldafoy, "../data/lda-foy.rds")
