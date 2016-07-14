## plot LDA and old models for Cunswick Tarn

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
if (any(!file.exists(c("../data/lda-cuns.rds", "../data/ldagibbs-cuns.rds")))) {
  source("lda-cuns.R")
  }
ldacuns <- readRDS("../data/lda-cuns.rds")
ldacunsgibbs <- readRDS("../data/ldagibbs-cuns.rds")


## create suitable data frames
agedf <- data.frame(cbind(ldacuns@gamma, topics(ldacuns), cladoyears))
colnames <- c(as.character(seq_len(ldacuns@k)), "Cluster", "Year")
names(agedf) <- colnames
agedf$Cluster <- factor(agedf$Cluster)

agestack <- melt(agedf, id.vars=c('Year','Cluster'), variable_name = 'Group')
agesplit <- with(agestack, split(agestack, list(Group)))

agedfgibbs <- data.frame(cbind(ldacunsgibbs@gamma, topics(ldacunsgibbs), cladoyears))
colnamesgibbs <- c(as.character(seq_len(ldacunsgibbs@k)), "Cluster", "Year")
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
cunsrels <- ggplot(data = agestack, aes(x=Year, y=value, col=Group, lty=Group)) +
  papertheme +
  geom_rect(data=rectdf, inherit.aes = FALSE, 
            aes(xmin = xmin, xmax= xmax, ymin=-Inf, ymax=Inf, group=factor(Cluster),
                fill=factor(Cluster)), alpha = 0.3) +
  scale_linetype_manual(name='Species group', values = c("solid", "dotdash","solid", 
                                                         "solid", "longdash", "dotdash")) +
  
  scale_colour_manual(name="Species group", values = c("#b2abd2", "#5e3c99","#5e3c99", 
                                                       "#e66101", "#e66101", "#b2abd2"))+
  
  scale_fill_manual(name="Community", values = c("#D8DAEB", "#998EC3","#542788",
                                                 "#B35806", "#F1A340", "#FEE0B6"))+
  #geom_vline(xintercept = c(1390, 1900, 1930, 2000), col='grey50') +
  geom_line() +
  theme(legend.box = "horizontal") +
  guides(fill=guide_legend(nrow=1,bycol =TRUE,title.position = 'left'),
         lty=guide_legend(nrow=1, bycol =TRUE,title.position = 'left')) +
  ylab("Proportion of population")

cunsrelsgibbs <- ggplot(data = agestackgibbs, aes(x=Year, y=value, col=Group, lty=Group)) +
  papertheme +
  geom_rect(data=rectdfgibbs, inherit.aes = FALSE, 
            aes(xmin = xmin, xmax= xmax, ymin=-Inf, ymax=Inf, group=factor(Cluster),
                fill=factor(Cluster)), alpha = 0.3) +
  scale_linetype_manual(name='Species group', values = c("solid", "dotdash","solid", 
                                                         "solid", "longdash", "dotdash",
                                                         "solid", "longdash", "dotdash")) +
  
  scale_colour_manual(name="Species group", values = c("#b2abd2", "#5e3c99","#5e3c99", 
                                                       "#e66101", "#e66101", "#b2abd2",
                                                       "grey50", "grey50", "yellow"))+
  
  scale_fill_manual(name="Community", values = c("#D8DAEB", "#998EC3","#542788",
                                                 "#B35806", "#F1A340", "#FEE0B6",
                                                 "grey50", "grey50", "yellow"))+
  #geom_vline(xintercept = c(1390, 1900, 1930, 2000), col='grey50') +
  geom_line() +
  theme(legend.box = "horizontal") +
  guides(fill=guide_legend(nrow=1,bycol =TRUE,title.position = 'left'),
         lty=guide_legend(nrow=1, bycol =TRUE,title.position = 'left')) +
  ylab("Proportion of population")

## create bubble plot of groups
cunstermit <- ggplot(data = agestack, aes(y=Year, x=Group, col=Cluster, size=value)) +
  papertheme +
  #scale_color_distiller(name="Value", palette = 'PuOr') +
  scale_size_continuous(name="Relative abundance") + #guide = FALSE
  scale_color_brewer(name="Community", palette = 'Dark2') +
  #geom_abline(intercept = c(1390, 1900, 1930, 2000), col='grey50', slope = 0) +
  geom_point(alpha=0.6) +
  theme(legend.box = "vertical") +
  xlab("Species group") +
  guides(colour=guide_legend(nrow=1,byrow=TRUE))

## create bubble plot of species
cladospec <- as.data.frame(posterior(ldacuns)$terms)
cladostack <- cbind(stack(cladospec), rep(factor(1:6), times = ncol(cladospec)))
names(cladostack)[3] <- "Cluster"
ggplot(data = cladostack, aes(y=ind, x=Cluster, size=values)) +
  papertheme +
  #scale_color_distiller(name="Value", palette = 'PuOr') +
  scale_size_continuous(name="Relative abundance") + #guide = FALSE
  scale_color_brewer(name="Community", palette = 'Dark2') +
  #geom_abline(intercept = c(1390, 1900, 1930, 2000), col='grey50', slope = 0) +
  geom_point(alpha=0.6) +
  theme(legend.box = "vertical") +
  xlab("Species group") +
  ylab(NULL) +
  guides(colour=guide_legend(nrow=1,byrow=TRUE))

## ============================================================================================
## SAVE PLOTS
## ============================================================================================
saveRDS(cunsrels, "../data/gg-cuns-line.rds")
saveRDS(cunstermit, "../data/gg-cuns-bubble.rds")
