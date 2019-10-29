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
    print("loop ends")

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
            res_Nclust <- rbind(res_Nclust,
                                data.frame(
                                    cluster_id = colX,
                                    pc = pc,
                                    res = resX,
                                    Num_Clus = (1 + as.numeric(max(
                                        as.numeric(levels(Object[[]][, colX]))
                                    ) ) )
                      ))
        }
    }
    return(result)
}
