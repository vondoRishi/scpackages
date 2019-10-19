#' Title
#'
#' @param tenxPath
#' @param pMin.cells
#' @param pMin.features
#' @param markerGenes
#' @param mitoPattern
#' @param projectName
#' @param cellRangerAggregated
#'
#' @return
#' @export
#'
#' @import Seurat stringr
#'
#' @examples
read10XwithMarkergenes <-  function(tenxPath,
           pMin.cells = 3, pMin.features = 20,
           markerGenes, mitoPattern = "^mt-",
           projectName, cellRangerAggregated =TRUE) {

  if( dir.exists(path)){
    AGG.data <- Read10X(data.dir =  path )
  }else{
    AGG.data <- Read10X_h5(filename = path )
  }
  AGG <- CreateSeuratObject( counts = AGG.data,
                             min.cells = pMin.cells,
                             min.features = pMin.features,
                             project = projectName
                             )
  AGG[["percent.mito"]] <- PercentageFeatureSet(AGG, pattern = mitoPattern)
  for (genes in markerGenes) {
    AGG[[paste("percent", genes, sep = ".")]] <-
      PercentageFeatureSet(AGG, pattern = paste("^", genes, "$", sep = ""))
  }

  AGG@misc$messages <-list()
  AGG@misc$messages <- str_c(AGG@misc$messages, paste("Data loaded from",path),collapse = " ")
  AGG@misc$messages <- c(AGG@misc$messages, paste(capture.output(show(AGG)),collapse = "") )
  group_by <- NULL
  if(cellRangerAggregated){
    AGG <- ExtractSamplesFromAggCellRanger(Object = AGG)
    AGG@misc$cellRangerAggregated <- cellRangerAggregated
    group_by <- "sample"
  }
  AGG <- QC_plot(AGG,group.by = group_by)
  return(AGG)
}

ExtractSamplesFromAggCellRanger <- function(Object, splitBy="-") {
  Object[["old.ident"]] <- Idents(object = Object)
  Object[["sample"]]<- str_split(colnames(Object),pattern = splitBy,simplify = TRUE)[,2]
  Idents(object = Object) <- "sample"
  return(Object)
}

#' Title
#'
#' @param Object
#' @param group.by
#'
#' @return
#' @export
#'
#' @import ggplot2 cowplot ggpubr
#' @examples
QC_plot <- function(Object, group.by = NULL) {

  QC_umi_mito <-
    ggplot(Object[[]],
           aes(x = nCount_RNA, color = nFeature_RNA, y = percent.mito)) + geom_point() +
    theme(legend.position = "top")+theme_bw()

  QC_umi_gene <-
    ggplot(Object[[]],
           aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mito)) + geom_point() +
    theme(legend.position = "top") +theme_bw()
  xdens <- axis_canvas(QC_umi_gene, axis = "x")+
    geom_density(data = Object[[]], aes(x = nCount_RNA, fill = !!sym(group.by)),
                 alpha = 0.7, size = 0.2)+
    ggpubr::fill_palette("jco")
  ydens <- axis_canvas(QC_umi_gene, axis = "y", coord_flip = TRUE)+
    geom_density(data = Object[[]], aes(x = nFeature_RNA, fill = !!sym(group.by)),
                 alpha = 0.7, size = 0.2)+
    coord_flip()+
    ggpubr::fill_palette("jco")

  p1 <- insert_xaxis_grob(QC_umi_gene, xdens, grid::unit(.2, "null"), position = "top")
  p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
  QC_umi_gene<- ggdraw(p2)


  Object@misc$QC_plot <-
    CombinePlots(plots = list(QC_umi_mito, QC_umi_gene))

  if (is.null(group.by)) {
    Object@misc$QC_vln_plot <-
      VlnPlot(
        Object,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
        ncol = 3
      )
  } else{
    Object@misc$QC_vln_plot <-
      VlnPlot(
        Object,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
        ncol = 3,
        group.by = group.by
      )
  }

  return(Object)
}

#' Title
#'
#' @param Object
#' @param mito.range
#' @param gene_range
#' @param umi_range
#'
#' @return
#' @export
#'
#' @import Seurat
#'
#' @examples
filterSeurat <-
  function(Object,
           mito.range = c(0, 100),
           gene_range = c(0, Inf),
           umi_range = c(0, Inf)) {
    filter_message = paste(
      "Keeping cells with mitochndrial genes between",paste(mito.range, collapse = ","),"%",
      " gene count within",paste(gene_range, collapse = ","),
      " and umi count within",paste(umi_range, collapse = ","),sep = ""
    )
    Object@misc$messages <- c(Object@misc$messages, filter_message )
    x<-  subset(Object[[]],
                percent.mito > min(mito.range) &
                  percent.mito < max(mito.range) &
                  nCount_RNA > min(gene_range) &
                  nCount_RNA < max(gene_range) &
                  nFeature_RNA > min(umi_range) &
                  nFeature_RNA < max(umi_range)
    )
    Object <- subset( Object, cells= rownames(x))
    Object@misc$messages <- c(Object@misc$messages, paste("After filter ",capture.output(show(Object))[2],collapse = "") )
    return(Object)
  }


topFeatures <- function(Object, selection.method, topN) {
  hvg_info <- Object[["RNA"]]@meta.features
  if(selection.method == "vst"){
    return(rownames(head(
      hvg_info[order(-hvg_info$vst.variance.standardized),] %>% subset(vst.variable ==
                                                                         TRUE), n = topN
    )))
  }else if(selection.method == "mvp"){
    return(rownames(head(
      hvg_info[order(-hvg_info$mvp.dispersion.scaled), ] %>% subset(mvp.variable ==
                                                                      TRUE), n = topN
    )))
  }
}

#' Title
#'
#' @param Object
#' @param features
#' @param jackStraw
#'
#' @return
#' @export
#'
#' @import Seurat
#'
#' @examples
pcaProcess <- function(Object, features,jackStraw= FALSE ){
  pDims = if(length(features) < 50 ) length(features) else 50
  Object <-
    RunPCA(
      Object,npcs = pDims,
      features = features
    )
  if(jackStraw){
      Object <- JackStraw(Object, num.replicate = 100, dims = pDims, verbose = FALSE)
      Object <- ScoreJackStraw(Object, dims = 1:pDims)
  }
  Object@misc$elbowPlot <- ElbowPlot(placodeObj_fltred,ndims = pDims)
  return(Object)
}

#' Title
#'
#' @param Object
#' @param maxDims
#' @param res
#' @param identName
#'
#' @return
#' @export
#'
#' @import Seurat
#' @examples
makeClusterSeurat <- function(Object, maxDims, res, identName=NULL){
  Object <- FindNeighbors(Object, dims = 1:maxDims)
  Object <- FindClusters(Object, resolution = res, verbose = FALSE)
  Object <- RunUMAP(Object, dims = 1:maxDims, verbose = FALSE)
  if(is.null(identName)){
    identName<-paste("Ident",maxDims,res,sep = "_")
  }
  Object[[identName]]<- Idents(Object)
  return(Object)
}

#' Title
#'
#' @param Object
#' @param marker_genes
#' @param cell_cycle_genes
#'
#' @return
#' @export
#'
#' @import Seurat dplyr
#'
#' @examples
clusterMarkers <-  function(Object, marker_genes = NULL,
                            cell_cycle_genes = NULL) {
  markers_info <-  FindAllMarkers( object = Object,
      only.pos = TRUE, min.pct = 0.25,
      thresh.use = 0.25, verbose = FALSE
    )

  if (!is.null(marker_genes)) {
    markers_info <- left_join(markers_info,
                              marker_genes, by = c("gene" = "Symbol")) %>%
      replace_na(list(Gene_type = "Non marker"))
  }

  if (!is.null(cell_cycle_genes)) {
    markers_info <-
      markers_info %>% left_join(cell_cycle_genes, by = c("gene" = "Symbol")) %>%
      replace_na(list(Phase = "Not CC", Gene_type = "Non marker"))
  }
  return(markers_info)
}

#' Title
#'
#' @param ClusterInfo
#' @param Object
#'
#' @return
#' @export
#'
#' @import Seurat dplyr
#'
#' @examples
renameCluster <- function(ClusterInfo, Object) {
  varCluster_top <- ClusterInfo %>% group_by(cluster) %>% top_n(10, avg_logFC)
  cluster_rename <- as.data.frame(table(varCluster_top$Gene_type,varCluster_top$cluster)) %>% group_by(Var2) %>% filter(Freq== max(Freq))
  new.name <- paste(cluster_rename$Var1, cluster_rename$Var2, sep="-")
  names(new.name) <- levels(placodeObj_fltred)
  Object <- RenameIdents(Object,new.name)
  return(Object)
}


pca_gene_loading <- function(Object) {
  raw_pca_loading <- Loadings(Object[["pca"]])
  max_pc<-colnames(raw_pca_loading)[apply(abs(raw_pca_loading),1,which.max)]
  raw_gene_pca_loading <- as.data.frame(cbind(Symbol=rownames(raw_pca_loading),PCA=max_pc))
  raw_gene_pca_loading <- raw_gene_pca_loading %>% mutate(pca_fix=as.numeric(substring(PCA,4)))
  return(raw_gene_pca_loading)
}


#' Title
#'
#' @param Object
#' @param CellCycle_genes
#'
#' @return
#' @export
#'
#' @import Seurat dplyr
#'
#' @examples
qc_regress_CellCycle<- function(Object,CellCycle_genes) {
  ### Check CC gene loadings
  Object <- RunPCA(Object, features = VariableFeatures(object = Object),verbose = FALSE)

  raw_gene_pca_loading <- pca_gene_loading(Object)

  raw_gene_pca_loading <-  left_join(raw_gene_pca_loading,CellCycle_genes,  by = "Symbol") %>%
    replace_na(list( Phase= "Not CC"))

  Object@misc$CC_pca_bar <- ggplot(subset(raw_gene_pca_loading, pca_fix <6),aes(x=factor(pca_fix)))+
    geom_bar(aes(fill=Phase))+
    ggtitle("CC genes on PCAs with default HVGs")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme_bw()
  ### Check CC gene loadings
  ### CellCycleScoring
  Object <-CellCycleScoring(Object,
      s.features = subset(CellCycle_genes, Phase == "S")$Symbol,
      g2m.features = subset(CellCycle_genes, Phase == "M")$Symbol,
      set.ident = TRUE
    )

  Object <- RunPCA(Object, verbose = FALSE,
                   features = subset(CellCycle_genes, Phase %in% c("S", "M"))$Symbol)
  Object@misc$before_cc_pca <- DimPlot(Object)+theme_bw() + ggtitle("Before regressing out")
  ### CellCycleScoring
  ### regress out CC
  Object <- ScaleData(Object,
    vars.to.regress = c("S.Score", "G2M.Score", "percent.mito", "nCount_RNA"),
    features = mvp_features,
    verbose = FALSE
  )

  Object <- RunPCA( Object, verbose = FALSE,
      features = VariableFeatures(Object)
    )

  Object@misc$after_cc_pca <- DimPlot(Object)+theme_bw() + ggtitle("After regressing out")
  ### regress out CC
  return(Object)
}


#' Title
#'
#' @param Object
#'
#' @return
#' @export
#'
#' import Seurat
#'
#' @examples
normScaleHVG <- function(Object,verbose =TRUE) {
  print("Normalizing")
  Object <- NormalizeData(    Object,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  print("Finding Variable Features")
  Object <- FindVariableFeatures(Object, selection.method = "vst", verbose = TRUE)
  vst_features <-VariableFeatures(Object,selection.method = "vst")
  Object@misc$hvgPlot <- topVariableFeaturePlot(Object,10)
  print(Object@misc$hvgPlot)

  print("scaling")
  Object <-  ScaleData(Object,features = vst_features,
                       vars.to.regress = c("percent.mito","nCount_RNA"),
                       verbose = TRUE
  )
}

topVariableFeaturePlot <- function(Object, topN){
  topGenes <- head(VariableFeatures(Object), topN)

  plot1 <- VariableFeaturePlot(Object)
  return(LabelPoints(plot = plot1, points = topGenes, repel = TRUE))
}
