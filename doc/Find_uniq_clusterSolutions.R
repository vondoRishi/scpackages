## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup,message=FALSE,warning=FALSE-----------------------------------
library(scpackages)

## ----readSeurat,warning=FALSE,fig.width=8,fig.height=6-------------------
scaled_cluster <- readRDS("../scaled_cluster.13_February_2020.rds")

## ----gene_p,warning=FALSE,fig.width=8,fig.height=6-----------------------
tmp_hvf <-HVFInfo(scaled_cluster)
tmp_hvf$Symbol <- rownames(tmp_hvf)
# write.csv(tmp_hvf, file = "HVF_info.csv",row.names = FALSE)

## gene PCA Explore ##
scaled_cluster@misc$loading <- scpackages::pca_gene_loading(scaled_cluster)
scaled_cluster@misc$loading <-
  scaled_cluster@misc$loading %>% left_join(tmp_hvf, by = "Symbol")

library(ggrepel)

ggplot(scaled_cluster@misc$loading, 
       aes(x=factor(pca_fix),y=variance.standardized, label = Symbol)) +
  geom_point(color = ifelse(scaled_cluster@misc$loading$variance.standardized > 4 
                            , "red", "black")) +
  geom_text_repel(
    data = subset(scaled_cluster@misc$loading, variance.standardized > 6 )
  )+ggtitle(label = "Visualize the gene loading in each PCA ")

## ----param,warning=FALSE,fig.width=8,fig.height=6------------------------

scaled_cluster@misc$elbowPlot+ggtitle("Where is the elbow")
## ??Where is the elbow does it matters ??

## ----optimal_param,warning=FALSE,fig.width=8,fig.height=6----------------
pcRange <- seq(15,25,by = 2)
resolutionRange <- seq(0.4,1.2,by = 0.2)

scaled_cluster <- scpackages::gridFindClusters( Object = scaled_cluster, 
                                                pcRange, resolutionRange,
                                                identPrefix = "solution", 
                                                verbose =FALSE)

ggplot(scaled_cluster@misc$data_NClust,
       aes(x = resolution, y = Cell_types)) + geom_line() + facet_wrap(. ~ pc)

## ----vis_solution,warning=FALSE,fig.width=8,fig.height=6-----------------
library(easyalluvial)
alluvial_wide(dplyr::select(
  scaled_cluster[[]],
  solution_17_0.4,
  solution_17_0.6,
  solution_17_0.8
)) + theme_bw()

## ----group_clusters,warning=FALSE,fig.width=8,fig.height=6---------------
scaled_cluster <-  scpackages::findSimilarClusterSolution(scaled_cluster,
                                         identPrefix = "solution",
                                         similarityCut = 0.8)
plot(scaled_cluster@misc$cluster_G1, vertex.color = 
       rainbow(10, alpha = 0.8)[scaled_cluster@misc$clusterlouvain$membership])


## ----warning=FALSE,fig.width=8,fig.height=6------------------------------
ggplot(  scaled_cluster@misc$meta_cluster,  aes(
                y = resolution,
                x = pca,
                # shape = factor(Cell_types),
                label = Cell_types,#membership,
                color = factor(membership)
              )
)  + geom_text()
#+  scale_color_manual(values =  rainbow(8)) #+  theme_dark()

## ----warning=FALSE,fig.width=8,fig.height=6------------------------------
cluster_representatives <- as_tibble( scaled_cluster@misc$meta_cluster %>% 
                                        group_by(membership,Cell_types) %>% 
                                        sample_n(1)) %>% arrange(pca,resolution)


x <- as.list(cluster_representatives %>% 
               filter( Cell_types==9 ) %>% 
               select(cluster_id))
alluvial_wide(scaled_cluster[[]][, as.character(x$cluster_id)]) + 
  theme(axis.text.x =  element_text(angle = 90)) +
  theme_bw() + 
  ggtitle(paste(    length(x$cluster_id),
    " different solutions with equal number of clusters"
  ))

