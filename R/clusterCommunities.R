#' Title
#'
#' @param Object
#' @param pcRange
#' @param resolutionRange
#' @param identPrefix
#'
#' @return
#' @export
#'
#' @import Seurat
#'
#' @examples
gridFindClusters <- function(Object,pcRange,resolutionRange, identPrefix) {
    for (pc in pcRange) {
        Object <- FindNeighbors(Object, dims = 1:pc)
        for (resX in resolutionRange) {
            Object <- FindClusters(Object, resolution = resX, verbose = FALSE)
            Object[[paste(identPrefix, pc, resX, sep = "_")]]<- Idents(Object)
        }
    }
    Object@misc$data_NClust <- res_Nclust(Object,pcRange,resolutionRange, identPrefix)
    return(Object)
}


#' Title
#'
#' @param Object
#' @param pcRange
#' @param resolutionRange
#' @param identPrefix
#'
#' @return
#' @export
#'
#' @import Seurat
#'
#' @examples
res_Nclust <- function(Object,pcRange,resolutionRange, identPrefix) {
    result <-data.frame(  cluster_id = character(),
                          pc = numeric(), res = numeric(), Num_Clus = numeric()
        )
    for (pc in pcRange) {
        for (resX in resolutionRange) {
            colX = paste(identPrefix, pc, resX, sep = "_")
            # print(max(levels( scaled_cluster[[]][,colX])))
            result <- rbind(result,
                                data.frame(
                                    cluster_id = colX,
                                    pc = pc,
                                    res = resX,
                                    Num_Clus = (1 + max(
                                        as.numeric(levels(Object[[]][, colX]))
                                    ) )
                      ))
        }
    }
    return(result)
}


#' Title
#'
#' @return
#' @export
#'
#' @import clues igraph
#'
#' @examples
findSimilarClusterSolution <- function(Object, identPrefix, similarityCut){
    cluster_results <- t( Object[[]][, grep(identPrefix,colnames(Object[[]]))])
    invert_adjustedRand <- function(v1, v2){  clues::adjustedRand(v1,v2)["HA"] }
    cluster_dist <- usedist::dist_make(cluster_results,invert_adjustedRand,"Adjusted rand index")


    df.dist_lou <- as.matrix(cluster_dist, labels=TRUE)
    df.dist_lou[df.dist_lou < similarityCut] =0
    # pheatmap::pheatmap(df.dist_lou,cluster_rows = F, cluster_cols = F, show_colnames = TRUE)

    Object@misc$cluster_G1 <- graph.adjacency(df.dist_lou, mode = "undirected", weighted = TRUE, diag = TRUE)

    Object@misc$clusterlouvain <- cluster_louvain(G1)
    return(Object)
}
