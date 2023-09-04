args = commandArgs(TRUE)
ematpath = as.character(args[1])
condpath = as.character(args[2])
grn_structure = as.character(args[3])
outtpath = as.character(args[4])


# setwd('/Users/fernando/Documents/Research/DraCooN/evaluation/simulated_data')
# ematpath =  '/Users/fernando/Documents/Research/DraCooN/evaluation/simulated_data/treelike/simulations/DR/perturbation_noise/pertnoise_0.5_run2_perts_0,0,10,_biomdata.csv'
# condpath =  '/Users/fernando/Documents/Research/DraCooN/evaluation/simulated_data/treelike/simulations/DR/perturbation_noise/pertnoise_0.5_run2_perts_0,0,10,_conddata.csv'
# grn_structure = '/Users/fernando/Documents/Research/DraCooN/evaluation/simulated_data/treelike/simulations/DR/perturbation_noise/pertnoise_0.5_run2_perts_0,0,10,_structure.csv'
# outtpath ='/Users/fernando/Documents/Research/DraCooN/evaluation/simulated_data/treelike/results/DR/perturbation_noise/pertnoise_0.5_run2_perts_0,0,10'

#grn_structure = 'simulations/DR/gene_number/GeneNum_N20_DR_run2_structure.csv'
print('loading necessary packages')
suppressMessages(library(igraph))
suppressMessages(library(dcanr))
suppressMessages(library(EBcoexpress))
suppressMessages(library(dplyr))

#library(COSINE)

print('loading datasets')
emat <- read.csv(file = ematpath, row.names = 1)
condition <- read.csv(file = condpath)
grn <- read.csv(file = grn_structure, row.names = 1)
grn <- grn[,c('source', 'target')]
grn <- graph_from_data_frame(grn, directed = F)
#visNetwork::visIgraph(grn)

list <- as.list(condition)
cond <- condition$condition
names(cond) <- condition$sample

#print(emat)
#print(cond)
print('running methods')
for (run in 0:4){
  for (method in dcMethods()) {
    if (method %in% c("ldgm", "mindy", "ecf")) next
    print(method)
    #if (method %in% c("ldgm", "mindy")) next
    z_scores <-dcScore(emat, cond, dc.method = method)
    raw_p <- dcTest(z_scores, emat, cond)
    adj_p <- dcAdjust(raw_p, f = p.adjust, method = 'fdr')
    # print(adj_p[1:5, 1:5])
    dcnet <- dcNetwork(z_scores,adj_p, thresh =1.1 )
    drnet <- dcnet %s% grn
    # visNetwork::visIgraph(drnet)
    #plot(dcnet, vertex.label = '')
    
    edgedf <- get.data.frame(drnet, what = 'edges')
    edgedf <- rename(edgedf, source = from, target = to)
      
    # Reshape df2 into a long format
    df2_melt <- reshape2::melt(adj_p)
    names(df2_melt) <- c("source", "target", "adj_pval")
    
    # Merge df1 and df2_melt
    edgedf <- merge(edgedf, df2_melt, by = c("source", "target"), all.x = TRUE)

    # If the 'adj_pval' is NA after merge, replace it with 0
    edgedf$adj_pval[is.na(edgedf$adj_pval)] <- 1

    file <- paste(outtpath,'_', method, '_algr', run, ".tsv")
    # print(file)
    write.csv(edgedf, file, row.names = F, quote = F, sep="\t")
  }
}


