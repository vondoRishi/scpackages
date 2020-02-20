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
gridFindClusters <- function(Object,pcRange,resolutionRange, identPrefix, ...) {
    for (pc in pcRange) {
        Object <- FindNeighbors(Object, dims = 1:pc, ...)
        for (resX in resolutionRange) {
            Object <- FindClusters(Object, resolution = resX,...)
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
                          pc = numeric(), resolution = numeric(), Cell_types = numeric()
        )
    for (pc in pcRange) {
        for (resX in resolutionRange) {
            colX = paste(identPrefix, pc, resX, sep = "_")
            result <- rbind(result,
                                data.frame(
                                    cluster_id = colX,
                                    pc = pc,
                                    resolution = resX,
                                    Cell_types = length(levels(Object[[]][, colX]))

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
#' @import clues igraph usedist
#'
#' @examples
findSimilarClusterSolution <- function(Object, identPrefix, similarityCut){
    cluster_results <- t( Object[[]][, grep(identPrefix,colnames(Object[[]]))])
    invert_adjustedRand <- function(v1, v2){  clues::adjustedRand(as.numeric(factor(v1)),
                                                                  as.numeric(factor(v2)))["HA"] }
    cluster_dist <- usedist::dist_make(cluster_results,invert_adjustedRand,"Adjusted rand index")

    df.dist_lou <- as.matrix(cluster_dist, labels=TRUE)
    df.dist_lou[df.dist_lou < similarityCut] =0

    Object@misc$cluster_G1 <- graph.adjacency(df.dist_lou, mode = "undirected", weighted = TRUE, diag = TRUE)
    Object@misc$clusterlouvain <- cluster_louvain(Object@misc$cluster_G1)

    meta_cluster <- data.frame(
            cluster_id = Object@misc$clusterlouvain$names,
            membership = Object@misc$clusterlouvain$membership
        )

    meta_cluster <- meta_cluster %>% separate(cluster_id,
                                  c("id_x", "pca", "resolution"),
                                  sep = "_",
                                  remove = FALSE) %>% select(-id_x)


    Object@misc$meta_cluster <-
        meta_cluster %>% left_join(Object@misc$data_NClust[, c("cluster_id", "Cell_types")])
    return(Object)
}
