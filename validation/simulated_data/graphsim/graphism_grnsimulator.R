

library("igraph")
library("gplots")
library("graphsim")
library("scales")

#grn_structure_path <- '/Users/fernando/Documents/Research/DraCooN/checkups/graphsim/graphsim_structure.csv'
#outpath <- '/Users/fernando/Documents/Research/DraCooN/checkups/graphsim/graphsim_simdata.csv'
#nsamples <- 200
#correlation_level <- 1
#mean_expression <- 5
#simdata_noise <- 1

args <- commandArgs(TRUE)
grn_structure_path <- as.character(args[1])
nsamples <- as.integer(args[2])
outpath <- as.character(args[3])
correlation_level <- as.numeric(args[4])
mean_expression <- as.numeric(args[4])
simdata_noise <- as.numeric(args[5])


structuredf <- read.csv(grn_structure_path)
#structuredf <- structuredf[,c('source', 'target')]

#library(visNetwork)
#visNetwork::visIgraph(igraph::graph_from_data_frame(structuredf))
graph_structure <- graph_from_data_frame(structuredf, directed = T)
E(graph_structure)$state <- structuredf$state

expr <- generate_expression(n = nsamples,
                            graph_structure,
                            cor = correlation_level,
                            mean = mean_expression,
                            sd = simdata_noise,
                            comm = FALSE,
                            dist = TRUE,
                            absolute = FALSE,
                            laplacian = FALSE)
write.csv(expr, outpath, quote = FALSE)





# heatmap.2(expr, scale = "none", trace = "none", col = bluered(50), colsep = 1:4, rowsep = 1:4)
corrmat <- cor(t(expr))
correlating_pairs <- NULL
correlating_pairs_genes <- NULL
edge_states <- structuredf$state
write.csv(expr, "expr_output.csv", row.names = FALSE)
write.csv(corrmat, "corrmat_output.csv", row.names = TRUE)

#
# for(i in 1:nrow(structuredf)) {
#   print(i)
#   pair <- c(structuredf[i, 'source'], structuredf[i, 'target'])
#   # print(pair)
#   correlating_pairs <- c(correlating_pairs, corrmat[pair[1],pair[2]])
#   correlating_pairs_genes <- c(correlating_pairs_genes, pair)
# }
#   plot_directed(graph_structure, layout = layout.kamada.kawai, state = edge_states, col.arrow = c("red", 'darkgreen')[edge_states / 2 + 1.5])
#
# not_correlating_pairs <- corrmat[corrmat %in% correlating_pairs == FALSE]
# not_correlating_pairs <- not_correlating_pairs[!duplicated(not_correlating_pairs)]
#
# df <- data.frame(rbind(cbind(not_correlating_pairs, rep('notcorr', length(not_correlating_pairs))),
#       cbind(correlating_pairs, rep('corr', length(correlating_pairs)))))
# colnames(df) <- c('correlation_value', 'type')
# df$correlation_value <- as.numeric(df$correlation_value)
# df$type <- as.factor(df$type)
# library(ggplot2)
# ggplot(df) +
#  aes(x = type, y = correlation_value, fill = type) +
#  geom_boxplot() +
#  scale_fill_hue(direction = 1) +
#  theme_minimal()

#cor(expr['g8',], expr['g12',])
#adj_mat <- make_adjmatrix_graph(graph_structure)
# laplacian_mat <- make_laplacian_graph(graph_structure)
# heatmap.2(make_laplacian_graph(graph_structure),
#           scale = "none", trace = "none",
#           col = bluered(50),colsep = 1:4, rowsep = 1:4)
#
# heatmap.2(make_adjmatrix_graph(graph_structure), scale = "none", trace = "none",
#           col = colorpanel(3, "grey75", "white", "blue"),
#           colsep = 1:4, rowsep = 1:4)
#
# heatmap.2(cor(t(expr)), scale = "none", trace = "none",
#            col = bluered(50), colsep = 1:4, rowsep = 1:4)
# }
