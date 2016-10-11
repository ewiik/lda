## plot LDA and old models for Foy Lake

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
if (any(!file.exists(c("../data/lda-foy.rds", "../data/ldagibbs-foy.rds")))) {
  source("lda-foy.R")
}
ldafoy <- readRDS("../data/lda-foy.rds")
ldafoygibbs <- readRDS("../data/ldagibbs-foy.rds")
foyyears <- readRDS("../data/foyyears.rds")

## create suitable data frames
agedf <- data.frame(cbind(ldafoy@gamma, topics(ldafoy), foyyears))
colnames <- c(as.character(seq_len(ldafoy@k)), "Cluster", "Year")
names(agedf) <- colnames
agedf$Cluster <- factor(agedf$Cluster)

agestack <- melt(agedf, id.vars=c('Year','Cluster'), variable_name = 'Group')
agesplit <- with(agestack, split(agestack, list('Group')))
#lapply(agesplit, summary)

agedfgibbs <- data.frame(cbind(ldafoygibbs@gamma, topics(ldafoygibbs), foyyears))
colnamesgibbs <- c(as.character(seq_len(ldafoygibbs@k)), "Cluster", "Year")
names(agedfgibbs) <- colnamesgibbs
agedfgibbs$Cluster <- factor(agedfgibbs$Cluster)

agestackgibbs <- melt(agedfgibbs, id.vars=c('Year','Cluster'), variable_name = 'Group')
agesplitgibbs <- with(agestackgibbs, split(agestackgibbs, list('Group')))

## create data frame to colour backtround by Community
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

## create basic point/line plot of groups
foyrels <- ggplot(data = agestack, aes(x=Year, y=value)) + #, lty=Group
  papertheme +
  geom_rect(data=rectdf, inherit.aes = FALSE, 
            aes(xmin = xmin, xmax= xmax, ymin=-Inf, ymax=Inf, group=factor(Cluster),
                fill=factor(Cluster)), alpha = 0.3) +
  #scale_linetype_manual(name='Species group', values = c("solid", "longdash","dotdash",
  #                                              "longdash", "dotdash")) +
  #geom_point(aes(col=variable), alpha=0.5, size=0.7) +
  scale_colour_manual(name="Species group", values = c("#5e3c99", "#e66101","#b2abd2", 
                                                       "black", "blue", "#A63603"))+
  scale_fill_manual(name="Community", values = c("#5e3c99", "#e66101","#b2abd2",
                                                 "#F7F7F7", "#3182BD", "#A63603"))+
  stat_smooth(method='loess', mapping=aes(col=variable), size=0.5, show.legend = TRUE, span=0.01, se=FALSE) + 
  theme(legend.box = "horizontal") +
  ylab("Proportion of population \n (LOESS fit)") + xlab("Year (BP)") +
  guides(fill=guide_legend(nrow=1,bycol =TRUE,title.position = 'left'),
         col=guide_legend(nrow=1, bycol =TRUE,title.position = 'left')) 
  #guides(colour = guide_legend(override.aes = list(alpha = 1)))
## FIXME: how to take out lines from the fill legend?

## FIXME: not functional, awaiting resolution on n(groups) for gibbs
foyrelsgibbs <- ggplot(data = agestackgibbs, aes(x=Year, y=value, lty=variable)) +
  papertheme +
  geom_rect(data=rectdfgibbs, inherit.aes = FALSE, 
            aes(xmin = xmin, xmax= xmax, ymin=-Inf, ymax=Inf, group=factor(Cluster),
                fill=factor(Cluster)), alpha = 0.3) +
  scale_linetype_manual(name='Species Group', values = c("solid", "longdash","dotdash")) +
  scale_fill_manual(name="Community", values = c("#5e3c99", "#e66101","#b2abd2"))+
  geom_line() +
  #geom_point(aes(col=Cluster)) +
  theme(legend.box = "horizontal") +
  ylab("Proportion of population") + xlab("Year (BP)")

## create bubble plot of groups
foytermit <- ggplot(data = agestack, aes(y=variable, x=Year, col=factor(Cluster), size=value)) +
  papertheme +
  #scale_color_distiller(name="Value", palette = 'PuOr') +
  scale_size_continuous(name = "Relative abundance") + #guide = FALSE
  scale_color_brewer(name="Community", palette = 'Dark2') +
  #geom_abline(intercept = c(1390, 1900, 1930, 2000), col='grey50', slope = 0) +
  geom_point(alpha=0.3) +
  theme(legend.box = "vertical") +
  ylab("Species group") + xlab("Year (BP)") +
  guides(colour=guide_legend(nrow=1,byrow=TRUE, override.aes = list(alpha = 1)))

## bubble plot of species
diatspec <- as.data.frame(t(posterior(ldafoy)$terms))
diatmax <- do.call(rbind, as.list(apply(FUN = max, MARGIN=1, X=diatspec)))
removenames <- rownames(diatmax)[which(diatmax < 0.01)]
remove <- which(diatmax < 0.01)
diatred <- diatspec[-remove,]

diatstack <- cbind(stack(diatred), rep(rownames(diatred), times = ldafoy@k))
names(diatstack)[2] <- "Cluster"
names(diatstack)[3] <- "ind"

## can we order them?
speclist <- vector(mode="list")
for (i in 1:ldafoy@k) {
  speclist[[i]] <- diatstack$ind[which(diatstack$values > 0.18 & diatstack$Cluster == i)]
}
speclist <- unlist(speclist)
remove <- which(duplicated(speclist))
speclist <- as.character(speclist[-remove])
allspec <- c(speclist, ldafoy@terms)
allspec <- allspec[-which(duplicated(allspec))]
bubbleorder <- match(allspec, unique(diatstack$ind))
diatstack$ind <-factor(diatstack$ind, levels=unique(diatstack$ind)[c(bubbleorder)])

speciesbubble <- ggplot(data = diatstack, aes(y=ind, x=Cluster, size=values)) +
  papertheme +
  scale_size_continuous(name="Relative abundance") + #guide = FALSE
  geom_point(alpha=0.6) +
  theme(legend.box = "vertical") +
  xlab("Species group") +
  ylab(NULL) +
  guides(colour=guide_legend(nrow=1,byrow=TRUE))


## ============================================================================================
## SAVE PLOTS
## ============================================================================================
saveRDS(foyrels, "../data/private/gg-foy-line.rds")
saveRDS(foytermit, "../data/private/gg-foy-bubble.rds")
saveRDS(speciesbubble, "../data/private/gg-foy-spbubble.rds")
