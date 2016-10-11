## load packages
library('ggplot2')
library('gridExtra')

## read in data
if (any(!file.exists("../data/gg-cuns-line.rds", "../data/gg-foy-line.rds",
                     "../data/private/gg-gull-line.rds"))) {
  source("plots-cuns.R")
  source("plots-foy.R")
  source("plots-gull.R")
}
cunsrels <- readRDS("../data/private/gg-cuns-line.rds")
cunstermit <- readRDS("../data/private/gg-cuns-bubble.rds")
cunsbubble <- readRDS("../data/private/gg-cuns-spbubble.rds")

foyrels <- readRDS("../data/private/gg-foy-line.rds")
foytermit <- readRDS("../data/private/gg-foy-bubble.rds")
foybubble <- readRDS("../data/private/gg-foy-spbubble.rds")

gullrels <- readRDS("../data/private/gg-gull-line.rds")
gulltermit <- readRDS("../data/private/gg-gull-bubble.rds")
gullbubble <- readRDS("../data/private/gg-gull-spbubble.rds")

## =========================================================================================
## Arrange all plots
## =========================================================================================
allrels <- plot_grid(cunsrels, gullrels, foyrels, ncol = 1, 
                     labels=c('Cunswick Tarn', "Gull Lake", "Foy Lake"))

termits <- plot_grid(cunstermit, gulltermit, foytermit, ncol = 3, 
                     labels = c('Cunswick Tarn', "Gull Lake", "Foy Lake"))
## FIXME: how to disable the fill in the size legend?

ggsave("../docs/private/allrels.pdf", allrels, scale=1.5) #width=28, height=18, units = 'cm'
ggsave("../docs/private/termits.pdf", termits, width=28, height=13, units = 'cm',
       scale=1.5)

## =========================================================================================
## Species bubbles
## =========================================================================================
plots <- list(cunsbubble, gullbubble, foybubble)
plotnames <- c('cunsbubble', 'gullbubble', 'foybubble')

invisible( # this means I don't get the list [[1:3]] returned on screen
  lapply(
    seq_along(plots), 
    function(x) ggsave(filename=paste0("../docs/private/gg-bubbles-", plotnames[x], ".pdf"), 
                       plot=plots[[x]], scale=1.6) # width=7, height=5, units = 'in'
  ) )

