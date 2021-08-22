rm(list = ls())

library(xtable)
library(data.table)

filenames = list.files("./experiments/results/results_non_const", pattern="*", full.names=TRUE)

param.names = c("setup", "n", "p")

setup.values = c('A','B','C', "D")
raw = data.frame(t(sapply(filenames, function(fnm) {
  output = read.csv(fnm)[,-1]
  params = strsplit(fnm, "-")[[1]][2:4]
  if (params[[1]] != "D") {
    output = t(matrix(as.vector(as.matrix(output)), 5, 200))
  }

  mse.mean = colMeans(as.matrix(output))

  c(params,
    mse=sprintf("%.8f", round(mse.mean, 8)))
})))


rownames(raw) = 1:nrow(raw)
names(raw) = c(param.names,
               "did", "rt", "rs", "t", "ols")
raw = data.frame(raw)

options(stringsAsFactors = FALSE)
raw = data.frame(apply(raw, 1:2, as.character))

raw = raw[order(as.numeric(raw$p)),]
raw = raw[order(as.numeric(raw$n)),]
raw = raw[order(as.character(raw$setup)),]
rownames(raw) = 1:nrow(raw)

raw.round = raw
for (col in 4:8){
  raw.round[,col] <- round(as.numeric(unlist(raw[,col])),2)
}
raw = data.frame(apply(raw, 1:2, as.character))
raw.round = data.frame(apply(raw.round, 1:2, as.character))

# write raw csv output file
write.csv(raw.round, file="experiments/simulations/output_hte_full.csv")

# write latex tables
raw <- read.csv(file="experiments/simulations/output_hte_full.csv", header=TRUE, sep=",")
raw = raw[,-1]
#raw = round(raw,3)
tab.all = cbind("", raw)
rmse.idx = c(4:9)
for(iter in 1:nrow(tab.all)) {
  best.idx = rmse.idx[which(as.numeric(tab.all[iter,rmse.idx]) == min(as.numeric(tab.all[iter,rmse.idx])))]
  tab.all[iter,best.idx] = paste("\\bf", tab.all[iter,best.idx])
}

tab.all=tab.all[,-1]
xtab.all = xtable(tab.all, align=rep("c",9))
print(xtab.all, include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = identity
      , file = "simulation_results_hte.tex")

