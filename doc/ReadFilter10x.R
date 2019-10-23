## ----setup,message=FALSE,warning=FALSE-----------------------------------
library(scpackages)

## ----warning=FALSE,fig.width=8,fig.height=6------------------------------
seurat_Obj <- read10XwithMarkergenes(
    tenxPath =  "/tmp/filtered_gene_bc_matrices/mm10/",
    pMin.cells = 3,
    pMin.features = 20,
    markerGenes = NULL,
    projectName = "1k_Brain_E18_Mm", 
    cellRangerAggregated = FALSE # If 10x data are aggrgated
)

## ----filter,warning=FALSE,fig.width=8,fig.height=6-----------------------
seurat_Obj <- filterSeurat(seurat_Obj, 
                           mito.range = c(0, 5), # Cells with mitochondrial genes 0-5%
                           gene_range = c(1000,Inf) # Cells with more than 1000 genes 
                           )

## ----report,message=FALSE,warning=FALSE,fig.width=8,fig.height=6---------
report_QC(obj_scRNA = seurat_Obj,title = "1k_Brain_E18_Mm")

