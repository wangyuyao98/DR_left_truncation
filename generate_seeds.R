# Generate random seeds for simulating data sets

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

set.seed(123)
n.seed = 500
seeds = round(runif(n.seed, 0, 100000))
save(seeds, file = "seeds.rda")
