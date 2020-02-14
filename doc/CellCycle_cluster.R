## ----setup,message=FALSE,warning=FALSE-----------------------------------
library(scpackages)

## ----readSeurat,warning=FALSE,fig.width=8,fig.height=6-------------------
seuratObj <- readRDS("/home/dasroy/data/Project/scpackages/1k_Brain_E18_Mm.12_February_2020.rds")

## ----normScaleHVG,warning=FALSE,fig.width=8,fig.height=6-----------------
seuratObj <- normScaleHVG(seuratObj,seuratVerbose = FALSE,nfeatures=1000)

## ----hvg_info,eval=FALSE-------------------------------------------------
#  tmp_hvf <-HVFInfo(seuratObj)
#  tmp_hvf$Symbol <- rownames(tmp_hvf)
#  write.csv(tmp_hvf, file = "HVF_info.csv",row.names = FALSE)

## ----CellCycle,message=FALSE,warning=FALSE,fig.width=8,fig.height=6------
library(readr)
cell_cycle_genes <- read_csv(system.file("extdata", "cell_cycle_genes.csv", package = "scpackages"))
cell_cycle_genes

## ----CellCycle_regress,message=FALSE,warning=FALSE,fig.width=8,fig.height=6----
seuratObj <-qc_regress_CellCycle(seuratObj, cell_cycle_genes)

## ----Elbowplot,message=FALSE,warning=FALSE,fig.width=8,fig.height=6------
seuratObj <- pcaProcess(seuratObj,features = VariableFeatures(seuratObj),jackStraw = FALSE)
# seuratObj@misc$elbowPlot

## ----Clustering,message=FALSE,warning=FALSE,fig.width=8,fig.height=6-----
seuratObj <- makeClusterSeurat(seuratObj,maxDims = 20, res = 0.5)
# seuratObj@misc$umapPlot

## ----Clustering_markers,message=FALSE,warning=FALSE,fig.width=8,fig.height=6----
markers_info <-  FindAllMarkers( object = seuratObj,
                                 only.pos = TRUE, min.pct = 0.25,
                                 thresh.use = 0.25, verbose = FALSE )
head(markers_info)
seuratObj@misc$markers_info <- markers_info

## ----echo=TRUE,eval=FALSE------------------------------------------------
#  write.csv(markers_info, file = "cluster_markers.csv",row.names = FALSE)

## ----Dummy,message=FALSE,warning=FALSE,fig.width=8,fig.height=6----------
marker_gene <- read_csv("/home/dasroy/data/Project/GeneExpression/scrna_workflow/marker_genes.csv")
marker_gene <- marker_gene[!duplicated(marker_gene$Symbol),]
head(marker_gene)

## ----Annotation,message=FALSE,warning=FALSE,fig.width=8,fig.height=6-----
library(dplyr)
markers_info <- left_join(markers_info,marker_gene, 
                          by = c("gene" = "Symbol")) %>% 
    replace_na(list(Gene_type = "Non skin"))

# Storing the info 
seuratObj@misc$markers_info <- markers_info


## ----Rename,message=FALSE,warning=FALSE,fig.width=8,fig.height=6---------
seuratObj <- renameCluster(ClusterInfo = markers_info, Object = seuratObj)

## ----report,message=FALSE,warning=FALSE,fig.width=8,fig.height=6---------
report_Cluster(obj_scRNA = seuratObj,title = "../scaled_cluster")

