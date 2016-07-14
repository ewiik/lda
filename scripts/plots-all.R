## read in data
if (any(!file.exists("../data/gg-cuns-line.rds", "../data/gg-foy-line.rds",
                     "../data/private/gg-gull-line.rds"))) {
  source("plots-cuns.R")
  source("plots-foy.R")
  source("plots-gull.R")
}
cunsrels <- readRDS("../data/gg-cuns-line.rds")
cunstermit <- readRDS("../data/gg-cuns-bubble.rds")

foyrels <- readRDS("../data/gg-foy-line.rds")
foytermit <- readRDS("../data/gg-foy-bubble.rds")

gullrels <- readRDS("../data/private/gg-gull-line.rds")
gulltermit <- readRDS("../data/private/gg-gull-bubble.rds")

## =========================================================================================
## Arrange all plots
## =========================================================================================
allrels <- plot_grid(cunsrels, gullrels, foyrels, ncol = 1, 
                     labels=c('Cunswick Tarn', "Gull Lake", "Foy Lake"))

termits <- plot_grid(cunstermit, gulltermit, foytermit, ncol = 3, 
                     labels = c('Cunswick Tarn', "Gull Lake", "Foy Lake"))

ggsave("../docs/private/allrels.pdf", allrels, scale=1.5) #width=28, height=18, units = 'cm'
ggsave("../docs/private/termits.pdf", termits, width=28, height=13, units = 'cm',
       scale=1.5)
