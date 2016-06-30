## plot all the models

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
if (any(!file.exists("../data/lda-foy.rds", "../data/lda-clado.rds",
                     "../data/private/lda-diat.rds"))) {
  source("lda-clado.R")
  source("lda-foy.R")
  source("lda-diat.R")
}
ldafoy <- readRDS("../data/lda-foy.rds")
ldaclado <- readRDS("../data/lda-clado.rds")
ldadiat <- readRDS("../data/private/lda-diat.rds")

## ============================================================================================
## Foy data
## ============================================================================================
agedf <- data.frame(cbind(ldafoy@gamma, topics(ldafoy), foyyears))
colnames <- c(as.character(seq_len(ldafoy@k)), "Cluster", "Year")
names(agedf) <- colnames
agedf$Cluster <- factor(agedf$Cluster)

agestack <- melt(agedf, id.vars=c('Year','Cluster'), variable_name = 'Group')
agesplit <- with(agestack, split(agestack, list(Group)))
#lapply(agesplit, summary)

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

## ============================================================================================
## Gull data
## ============================================================================================
agedf <- data.frame(cbind(ldadiat@gamma, topics(ldadiat), alldiats$AGE))
colnames <- c(as.character(seq_len(ldadiat@k)), "Cluster", "Year")
names(agedf) <- colnames
agedf$Cluster <- factor(agedf$Cluster)

agestack <- melt(agedf, id.vars=c('Year','Cluster'), variable_name = 'Group')
agesplit <- with(agestack, split(agestack, list(Group)))
#lapply(agesplit, summary)

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

diatrels <- ggplot(data = agestack, aes(x=Year, y=value, lty=Group)) +
  papertheme +
  geom_rect(data=rectdf, inherit.aes = FALSE, 
            aes(xmin = xmin, xmax= xmax, ymin=-Inf, ymax=Inf, group=factor(Cluster),
                fill=factor(Cluster)), alpha = 0.3) +
  scale_linetype_manual(name='Species Group', values = c("solid", "longdash","dotdash")) +
  scale_fill_manual(name="Community", values = c("#5e3c99", "#e66101","#b2abd2"))+
  geom_line() +
  #geom_point(aes(col=Cluster)) +
  theme(legend.box = "horizontal") +
  ylab("Proportion of population")

diatermit <- ggplot(data = agestack, aes(y=Year, x=Group, col=factor(Cluster), size=value)) +
  papertheme +
  #scale_color_distiller(name="Value", palette = 'PuOr') +
  scale_size_continuous(name = "Relative abundance") + #guide = FALSE
  scale_color_brewer(name="Community", palette = 'Dark2') +
  #geom_abline(intercept = c(1390, 1900, 1930, 2000), col='grey50', slope = 0) +
  geom_point(alpha=0.6) +
  theme(legend.box = "vertical") +
  xlab("Species group") +
  guides(colour=guide_legend(nrow=1,byrow=TRUE))

diatspec <- as.data.frame(t(posterior(ldadiat)$terms))
diatmax <- do.call(rbind, as.list(apply(FUN = max, MARGIN=1, X=diatspec)))
removenames <- rownames(diatmax)[which(diatmax < 0.0001)]
remove <- which(diatmax < 0.0001)
diatred <- diatspec[-remove,]

diatstack <- cbind(stack(diatred), rep(rownames(diatred), times = 3))
names(diatstack)[3] <- "Species"
ggplot(data = diatstack, aes(y=Species, x=ind, size=values)) +
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

## ==========================================================================================
## Cunswick data
## ==========================================================================================
agedf <- data.frame(cbind(ldaclado@gamma, topics(ldaclado), cladoyears))
colnames <- c(as.character(seq_len(ldaclado@k)), "Cluster", "Year")
names(agedf) <- colnames
agedf$Cluster <- factor(agedf$Cluster)

agestack <- melt(agedf, id.vars=c('Year','Cluster'), variable_name = 'Group')
agesplit <- with(agestack, split(agestack, list(Group)))
#lapply(agesplit, summary)

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

cladorels <- ggplot(data = agestack, aes(x=Year, y=value, col=Group, lty=Group)) +
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
  

cladotermit <- ggplot(data = agestack, aes(y=Year, x=Group, col=Cluster, size=value)) +
  papertheme +
  #scale_color_distiller(name="Value", palette = 'PuOr') +
  scale_size_continuous(name="Relative abundance") + #guide = FALSE
  scale_color_brewer(name="Community", palette = 'Dark2') +
  #geom_abline(intercept = c(1390, 1900, 1930, 2000), col='grey50', slope = 0) +
  geom_point(alpha=0.6) +
  theme(legend.box = "vertical") +
  xlab("Species group") +
  guides(colour=guide_legend(nrow=1,byrow=TRUE))

cladospec <- as.data.frame(posterior(ldaclado)$terms)
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

## =========================================================================================
## Arrange all plots
## =========================================================================================
allrels <- plot_grid(cladorels, diatrels, foyrels, ncol = 1, 
                    labels=c('Cunswick Tarn', "Gull Lake", "Foy Lake"))

termits <- plot_grid(cladotermit, diatermit, foytermit, ncol = 3, 
                        labels = c('Cunswick Tarn', "Gull Lake", "Foy Lake"))

ggsave("../docs/private/allrels.pdf", allrels, scale=1.5) #width=28, height=18, units = 'cm'
ggsave("../docs/private/termits.pdf", termits, width=28, height=13, units = 'cm',
       scale=1.5)
