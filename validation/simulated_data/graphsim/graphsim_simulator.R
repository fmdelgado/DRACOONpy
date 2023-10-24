generate_expression_analysis <- function(grn_structure,
                                         nsamples,
                                         mean_expression,
                                         simdata_noise) {
    library("igraph")
    library("gplots")
    library("graphsim")
    library("scales")
    library("ggplot2")

    graph_structure <- graph_from_data_frame(grn_structure, directed = T)

    expr <- generate_expression(n = nsamples,
                                graph_structure,
                                cor = 1,
                                mean = mean_expression,
                                sd = simdata_noise,
                                comm = FALSE,
                                dist = TRUE,
                                absolute = FALSE,
                                laplacian = FALSE)
    return(expr)
}
