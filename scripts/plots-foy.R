## plot LDA and old models for Foy Lake

## load packages
library("topicmodels")
library("analogue")
library("ggplot2")
library("reshape")
library("cowplot")
library("gridExtra")
library("extrafont") # Arial innit


papertheme <- theme_bw(base_size=12, base_family = 'Arial') +
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
agesplit <- with(agestack, split(agestack, list(Group)))
#lapply(agesplit, summary)

agedfgibbs <- data.frame(cbind(ldafoygibbs@gamma, topics(ldafoygibbs), cladoyears))
colnamesgibbs <- c(as.character(seq_len(ldafoygibbs@k)), "Cluster", "Year")
names(agedfgibbs) <- colnamesgibbs
agedfgibbs$Cluster <- factor(agedfgibbs$Cluster)

agestackgibbs <- melt(agedfgibbs, id.vars=c('Year','Cluster'), variable_name = 'Group')
agesplitgibbs <- with(agestackgibbs, split(agestackgibbs, list(Group)))

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

## create basic line plot of groups
foyrels <- ggplot(data = agestack, aes(x=Year, y=value)) + #, lty=Group
  papertheme +
  geom_rect(data=rectdf, inherit.aes = FALSE, 
            aes(xmin = xmin, xmax= xmax, ymin=-Inf, ymax=Inf, group=factor(Cluster),
                fill=factor(Cluster)), alpha = 0.3) +
  #scale_linetype_manual(name='Species group', values = c("solid", "longdash","dotdash",
  #                                              "longdash", "dotdash")) +
  scale_colour_manual(name="Species group", values = c("#5e3c99", "#e66101","#b2abd2", 
                                                       "black", "blue"))+
  scale_fill_manual(name="Community", values = c("#5e3c99", "#e66101","#b2abd2",
                                                 "#F7F7F7", "#FDB863"))+
  #geom_line() +
  geom_point(aes(col=Group), alpha=0.5) +
  theme(legend.box = "horizontal") +
  ylab("Proportion of population")

foyrelsgibbs <- ggplot(data = agestackgibbs, aes(x=Year, y=value, lty=Group)) +
  papertheme +
  geom_rect(data=rectdfgibbs, inherit.aes = FALSE, 
            aes(xmin = xmin, xmax= xmax, ymin=-Inf, ymax=Inf, group=factor(Cluster),
                fill=factor(Cluster)), alpha = 0.3) +
  scale_linetype_manual(name='Species Group', values = c("solid", "longdash","dotdash")) +
  scale_fill_manual(name="Community", values = c("#5e3c99", "#e66101","#b2abd2"))+
  geom_line() +
  #geom_point(aes(col=Cluster)) +
  theme(legend.box = "horizontal") +
  ylab("Proportion of population")

## create bubble plot of groups
foytermit <- ggplot(data = agestack, aes(y=Year, x=Group, col=factor(Cluster), size=value)) +
  papertheme +
  #scale_color_distiller(name="Value", palette = 'PuOr') +
  scale_size_continuous(name = "Relative abundance") + #guide = FALSE
  scale_color_brewer(name="Community", palette = 'Dark2') +
  #geom_abline(intercept = c(1390, 1900, 1930, 2000), col='grey50', slope = 0) +
  geom_point(alpha=0.3) +
  theme(legend.box = "vertical") +
  xlab("Species group") +
  guides(colour=guide_legend(nrow=1,byrow=TRUE))

## create bubble plot of species


## ============================================================================================
## SAVE PLOTS
## ============================================================================================
saveRDS(foyrels, "../data/gg-foy-line.rds")
saveRDS(foytermit, "../data/gg-foy-bubble.rds")

