## plot LDA and old models for Gull Lake

## load packages
library("topicmodels")
library("analogue")
library("ggplot2")
library("reshape2")
library("cowplot")
library("gridExtra")
library("extrafont") # Arial innit

## create a theme to save linespace in plots
papertheme <- theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top')

## read in data
if (any(!file.exists(c("../data/private/lda-gull.rds", "../data/private/lda-gull.rds")))) {
  source("lda-gull.R")
}
ldagull <- readRDS("../data/private/lda-gull.rds")
ldagullgibbs <- readRDS("../data/private/ldagibbs-gull.rds")

if (!file.exists('../data/private/7.5mCoreRawCounts.csv')) {
  stop("get Gull data from Heather or Emma")
}
alldiats <- read.csv('../data/private/7.5mCoreRawCounts.csv') 

## clean as per lda-gull.R
## remove things that aren't diatoms and change NA to 0
alldiats[is.na(alldiats)] <- 0
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

diatages <- alldiats$AGE - 1950 #make to BP
alldiats <- alldiats[,-notdiats]

## create data frame to colour backtround by Community
agedf <- data.frame(cbind(ldagull@gamma, topics(ldagull), diatages))
colnames <- c(as.character(seq_len(ldagull@k)), "Cluster", "Year")
names(agedf) <- colnames
agedf$Cluster <- factor(agedf$Cluster)

agestack <- melt(agedf, id.vars=c('Year','Cluster'), variable_name = 'Group')
agesplit <- with(agestack, split(agestack, list('Group')))
#lapply(agesplit, summary)

agedfgibbs <- data.frame(cbind(ldagullgibbs@gamma, topics(ldagullgibbs), diatages))
colnamesgibbs <- c(as.character(seq_len(ldagullgibbs@k)), "Cluster", "Year")
names(agedfgibbs) <- colnamesgibbs
agedfgibbs$Cluster <- factor(agedfgibbs$Cluster)

agestackgibbs <- melt(agedfgibbs, id.vars=c('Year','Cluster'), variable_name = 'Group')
agesplitgibbs <- with(agestackgibbs, split(agestackgibbs, list('Group')))

## create info to colour background
diffs <- which(abs(diff(as.numeric(agedf$Cluster))) > 0)

clusterbins <- agedf$Cluster[diffs] # cluster which changes in the next row

yearstarts <- agedf$Year[diffs] # year which changes cluster in next row
yearstops <- agedf$Year[diffs + 1] # year after cluster change
yearbins <- rowMeans(cbind(yearstarts, yearstops))

rectdf <- data.frame(xmin = c(agedf$Year[1], yearbins), 
                     xmax = c(yearbins, agedf$Year[nrow(agedf)]),
                     Cluster = c(clusterbins, agedf$Cluster[nrow(agedf)]))
##  --> http://stackoverflow.com/questions/26741703/adding-multiple-shadows-rectangles-to-ggplot2-graph

diffsgibbs <- which(abs(diff(as.numeric(agedfgibbs$Cluster))) > 0)

clusterbinsgibbs <- agedfgibbs$Cluster[diffsgibbs] # cluster which changes in the next row

yearstartsgibbs <- agedfgibbs$Year[diffsgibbs] # year which changes cluster in next row
yearstopsgibbs <- agedfgibbs$Year[diffsgibbs + 1] # year after cluster change
yearbinsgibbs <- rowMeans(cbind(yearstartsgibbs, yearstopsgibbs))

rectdfgibbs <- data.frame(xmin = c(agedfgibbs$Year[1], yearbinsgibbs), 
                          xmax = c(yearbinsgibbs, agedfgibbs$Year[nrow(agedfgibbs)]),
                          Cluster = c(clusterbinsgibbs, agedfgibbs$Cluster[nrow(agedfgibbs)]))

## create basic line plot of groups
gullrels <- ggplot(data = agestack, aes(x=Year, y=value, lty=variable)) +
  papertheme +
  geom_rect(data=rectdf, inherit.aes = FALSE, 
            aes(xmin = xmin, xmax= xmax, ymin=-Inf, ymax=Inf, group=factor(Cluster),
                fill=factor(Cluster)), alpha = 0.3) +
  scale_linetype_manual(name='Species Group', values = c("solid", "longdash","dotdash")) +
  scale_fill_manual(name="Community", values = c("#5e3c99", "#e66101","#b2abd2"))+
  geom_line() +
  #geom_point(aes(col=Cluster)) +
  theme(legend.box = "horizontal") +
  ylab("Proportion of population") + xlab("Year (BP)")

gullrelsgibbs <- ggplot(data = agestackgibbs, aes(x=Year, y=value, lty=variable, col=variable)) +
  papertheme +
  geom_rect(data=rectdfgibbs, inherit.aes = FALSE, 
            aes(xmin = xmin, xmax= xmax, ymin=-Inf, ymax=Inf, group=factor(Cluster),
                fill=factor(Cluster)), alpha = 0.3) +
  scale_linetype_manual(name='Species Group', values = c("solid", "longdash","dotdash",
                                                         "solid", "longdash","dotdash",
                                                         "dotdash", "dotdash","solid")) +
  scale_fill_manual(name="Community", values = c("#5e3c99", "#e66101","#b2abd2",
                                                 "grey50", "yellow","blue",
                                                 "#5e3c99", "#e66101","#b2abd2"))+
  scale_colour_manual(name="Species Group", values = c("#5e3c99", "#e66101","#b2abd2",
                                                 "grey50", "yellow","blue",
                                                 "#5e3c99", "#e66101","#b2abd2"))+
  geom_line() +
  #geom_point(aes(col=Cluster)) +
  theme(legend.box = "horizontal") +
  ylab("Proportion of population") +
  xlab("Year (BP)")

## create bubble plot of groups
gulltermit <- ggplot(data = agestack, aes(x=Year, y=variable, col=factor(Cluster), size=value)) +
  papertheme +
  #scale_color_distiller(name="Value", palette = 'PuOr') +
  scale_size_continuous(name = "Relative abundance") + #guide = FALSE
  scale_color_brewer(name="Community", palette = 'Dark2') +
  #geom_abline(intercept = c(1390, 1900, 1930, 2000), col='grey50', slope = 0) +
  geom_point(alpha=0.6) +
  theme(legend.box = "vertical") +
  ylab("Species group") +
  guides(colour=guide_legend(nrow=1,byrow=TRUE)) + xlab("Year (BP)")

## create bubble plot of species
diatspec <- as.data.frame(t(posterior(ldagull)$terms))
diatmax <- do.call(rbind, as.list(apply(FUN = max, MARGIN=1, X=diatspec)))
removenames <- rownames(diatmax)[which(diatmax < 0.01)]
remove <- which(diatmax < 0.01)
diatred <- diatspec[-remove,]

diatstack <- cbind(stack(diatred), rep(rownames(diatred), times = ldagull@k))
names(diatstack)[2] <- "Cluster"
names(diatstack)[3] <- "ind"
diatstack$ind <- as.character(diatstack$ind)

## can we order them?
speclist <- vector(mode="list")
for (i in 1:ldagull@k) {
  speclist[[i]] <- diatstack$ind[which(diatstack$values > 0.05 & diatstack$Cluster == i)]
}
speclist <- unlist(speclist)
remove <- which(duplicated(speclist))
speclist <- as.character(speclist[-remove])
allspec <- c(speclist, unique(diatstack$ind))
allspec <- allspec[-which(duplicated(allspec))]
bubbleorder <- match(allspec, unique(diatstack$ind))
diatstack$ind <-factor(diatstack$ind, levels=unique(diatstack$ind)[c(bubbleorder)])

speciesbubble <- ggplot(data = diatstack, aes(y=ind, x=Cluster, size=values)) +
  papertheme +
  scale_size_continuous(name="Relative abundance", breaks=c(0.01, 0.05, 0.1, 0.3)) + 
  geom_point(alpha=0.6) +
  theme(legend.box = "vertical") +
  xlab("Species group") +
  ylab(NULL)

## ===============================================================================
## SAVE PLOTS
## ===============================================================================
saveRDS(gullrels, "../data/private/gg-gull-line.rds")
saveRDS(gulltermit, "../data/private/gg-gull-bubble.rds")
saveRDS(speciesbubble, "../data/private/gg-gull-spbubble.rds")

