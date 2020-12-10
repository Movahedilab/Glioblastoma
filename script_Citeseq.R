suppressPackageStartupMessages({
 library(Seurat)
 library(dplyr)
 library(VennDiagram)
 library(harmony)
 library(clustree)
})

################################################################################
########## RNA data preprocessing
###############################################################################
#### Same as described in script_scRNAseq.R


################################################################################
########## ADT data preprocessing
###############################################################################

########## 1. Read the raw ADT data into a sparse matrix  ##########

#### Read the raw (unfiltered) antibody-barcode matrix from the
#### ADT Cellranger output of the CITEseq sample into a sparse matrix
ADT_raw <- Read10X("outs/raw_feature_bc_matrix/")[["Antibody Capture"]]


#### Remove barcodes with zero counts for all antibodies
dim(ADT_raw)
##### Calculate number of UMI counts per barcode
cell.nUMI<-colSums(ADT_raw)
ADT_raw<-ADT_raw[, cell.nUMI>0]
dim(ADT_raw)

#### Remove antibodies with zero counts for all cells
gene.total.nUMI<-rowSums(ADT_raw)
print("ADT features with 0 counts:")
print(rownames(ADT_raw[gene.total.nUMI==0,]))
### Remove antobodies with 0 counts
print(dim(ADT_raw))
ADT_raw<-ADT_raw[ADT.total.nUMI>0,]
print(dim(ADT_raw))


########## 2. Perform ASINH_GEOM normalization ##########

ASINH_GEOM_normalisation<-function(ADT_raw_matrix){
  ##### Find the total number of antibodies
  num.Abs <- length(rownames(ADT_raw_matrix))
  ##### Calculate modified geometric mean counts per antibody in asinh scale
  geo.means <- unlist(lapply(1:num.Abs,function(x){
    sinh(mean(asinh(as.matrix(ADT_raw_matrix[x,]))))
  }))
  ##### Perform ASINH_GEOM normalisation
  ADT_norm_matrix <- ADT_raw_matrix/geo.means
  ADT_norm_matrix<- asinh(ADT_norm_matrix)
  ADT_norm_matrix
}

ADT_norm<-ASINH_GEOM_normalisation(ADT_raw)


########## 3. Make some QC plots (optional) ##########

#### Calculate for each cell the normalised UMI counts of the antibody with
#### highest expression for this cell (max_UMI_counts)
max_UMI_counts<-  as.vector(apply(ADT_norm,2,max,na.rm=TRUE))

#### Plot a sample profile plot: the distribution of the maximum UMI count per barcode
#### The y axis show the number of cells with maxUMI larger than the respective x-value,
####  and the x -axis show the maxUMI values. The sample profile plot can be compared between
#### samples to reveal those with abnormal profile
A.tab.x <- sort(unique(max_UMI_counts))
A.tab.y  <- unlist(lapply(A.tab.x ,function(x){
  length(which(max_UMI_counts  >= x))
        }))
plot(A.tab.x,A.tab.y,xlab="max(UMI)",ylab="cell count >= UMI")


#### Compare the cells remaining after applying filtering with different thresholds on the ADT data with
#### the cells that passed the quality filtering of the RNA data for the same sample

#### Read the seurat object, based on the RNA data of the same CITEseq sample
seuratObj<-readRDS("path/to/seuratObj.rds")
seuratObj$barcode<-colnames(seur)

#### Define a function that filters the ADT matrix below a specified number of cells, based on the maximum UMI count per barcode
get_input <- function(numcells,tab.y=A.tab.y,max_UMI=max_UMI_counts,tab.x=A.tab.x,data_norm=ADT_norm){
  index_K_cells <- which(tab.y <= numcells)[1]
  input_clean <- which(max_UMI >= tab.x[index_K_cells])
  colnames(data_norm)[input_clean]
}

#### Define a function that plots a Venn Diagram of the overlap between cells passing the filtering of the ADT data
####vs those that pass the filteirng of the RNA data
RNA_ADT_venndiagram <- function(cells.input,seur.obj=seur,tab.y=A.tab.y,max_UMI=max_UMI_counts,tab.x=A.tab.x,data_norm=ADT_norm){
  ### cells.input - threshold for maximum number of cells after filtering
  ADT_barcodes_filt<-get_input(numcells=cells.input)
  perc<-(length(intersect(ADT_barcodes_filt,colnames(seur.obj)))/length(colnames(seur.obj)))*100
  v<-venn.diagram(
      list (ADT=ADT_barcodes_filt,
          RNA=colnames(seur.obj)
        ),
    main = paste0("Threshold: <=", cells.input, " cells (max cutoff) - ",round(perc,1),"% overlap with RNA 'true' cells"),
    filename=NULL,alpha = c( 0.5,0.5),cat.pos=c(270,100),margin=0.05,ext.text=T, fill = c("lightblue","orange")
   )
  grid.newpage()
  grid.draw(v)
}

#### Plot the Venn Diagram
RNA_ADT_venndiagram(cells.input=10000)
RNA_ADT_venndiagram(cells.input=20000)


########## 4. Filter cells ##########

#### Filter out the ADT data that don't overlap with RNA high quality cells
ADT_raw_filtered<-ADT_raw[,intersect(seur$barcode,colnames(ADT_raw))]
cat("raw ADT matrix:",dim(ADT_raw))
cat("\n intersect with RNA QC filtered matrix: ",dim(ADT_raw_filtered))


#### Rerun the ASINH_GEOM normalization on the filtered data
ADT_norm_filtered<-ASINH_GEOM_normalisation(ADT_raw_filtered)



################################################################################
########## Analyze jointly the ADT and RNA data
################################################################################

########## 1.Merge the ADT and RNA data into single seurat object ##########

##### Add the ADT data to the RNA seurat object
ADT_raw_filtered<-ADT_raw_filtered[,seur$barcode]
seur[["ADT"]]<-CreateAssayObject(counts =as.matrix(ADT_raw_filtered))
rownames( ADT_norm_filtered)<-rownames(seur[["ADT"]])
ADT_norm_filtered<-ADT_norm_filtered[,seur$barcode]
seur<-SetAssayData(object =seur,slot="data",new.data =ADT_norm_filtered, assay="ADT")
seur$max_ADT_counts<- as.vector(apply(as.matrix(GetAssayData(seur, assay="ADT", slot="counts")),2,max,na.rm=TRUE))

#### Merge seurat objects from several samples (optional)
seur$sample="sample1"
seur2=readRDS("path/to/sample2_SeuratObject.rds")
seur2$sample="sample2"
seur_merged<-merge(seur,seur2)


########## 2. ADT UMAP dimensionality reduction and clustering ##########

#### ADT feature selection
## Visualise the protein expression through density plots after filtering
ADT_data_plotting<-as.data.frame(t( as.matrix(GetAssayData(seur_merged, slot="data", assay="ADT"))))
ADT_data_plotting$sample<-seur_merged$sample[rownames(ADT_data_plotting)]
ADT_data_plotting<-reshape2::melt(ADT_data_plotting)
colnames(ADT_data_plotting)[colnames(ADT_data_plotting)=="variable"]<-"protein"
ggplot(ADT_data_plotting,  aes(x=value,color=sample, fill=sample))+
  geom_density(alpha=0.2)+
  facet_wrap(~protein)+xlab("Arsinh geom counts")+
  geom_rug()

## Visualise the protein expression on the RNA-based UMAP plot, with a 20-95 quantile cutoff for each feature
p=DimPlot(seur_merged,label=T, repel=T, group.by="cluster")
for (i in rownames(seur_merged[["ADT"]])){
  p1<-FeaturePlot(object = seur_merged, paste0("adt_",i), slot = "data" ,
     min.cutoff = "q20", max.cutoff = "q95",cols = c("yellow", "red"))+ggtitle(paste0(i, " arcsinh geom norm counts"))
  p2<-FeaturePlot(object = seur_merged, paste0("adt_",i) , slot = "counts",
    min.cutoff = "q20", max.cutoff = "q95",cols = c("yellow", "red"))+ggtitle(paste0(i, " counts"))
  print(plot_grid(p1,p2, ncol=2))
  print(scater::multiplot(p+NoLegend(),"", cols=2))
  plot.new()
  dev.off()
  rm(p1)
  rm(p2)
}
## Tag antibodies with low and non-specific expression in the studied sample, based on inspecting the density plots and feature plots
lowly.expr.ADT<-c("CD80",  "CD275-A0009" , "CD137L","CD70",  "CD30",  "CD277", "CD154", "CD34", "CD138-A0055" , "CD269", "CD56-A0084" ,  "TIGIT", "IgG2a-Mouse-k",
                  "IgG2b-Mouse-k" , "IgG2b-Rat-k" , "CD294","CD45R-B220" , "CD326", "CD133-A0126" ,  "CD140a" , "CD140b" ,  "Cadherin11", "CD340", "CD324" ,
                   "IgM",   "TCRgd", "CD183", "CD185", "CD197", "CD152", "CD223", "CD96", "CD178", "CD79b", "CD235a-CD235b" , "CD370", "XCR1",  "Notch1",
                  "CD268", "CD42b", "Notch3" , "IgG1-Rat-k", "IgG2a-Rat-k" , "IgG-Hamster", "CD122", "CD62E", "KLRG1-A0250","CD135", "CD254", "CD357", "IgG",
                  "B7-H4", "CD193", "CD301", "CD207-h",  "CD284", "CD243", "CD171", "CD230", "DopamineReceptorD4", "GABRB3",  "CD207-mh" , "CD338", "CD235a" ,
                  "CD79a", "mastCellTryptase",   "TCR-Vd2", "TCR-Vg9",  "CD202b",   "CD133-A0594", "CD110", "CD158f" ,   "CD337", "CD253", "CD186", "CD226-A0805",
                  "CD11a-CD18" ,  "CD307c" , "CD307d", "CD319", "CD138-A0831", "CD199", "CD45RB", "CD99",  "CD257", "CLEC1B", "IgE",   "CD365","Ig-lightChainK" ,
                  "Mac2",  "Ig-lightChainλ", "GARP",  "CD4-A0922" , "NKp80" )
## Select the remaining antibodies
features.use<-rownames(seur_merged[["ADT"]])[!rownames(seur_merged[["ADT"]]) %in% lowly.expr.ADT]

#### Scale ADT data
seur_merged<-ScaleData(seur_merged, assay = "ADT", verbose = F)

#### Run UMAP, using the scaled ARSINH GEOM transformed values (excluding lowly expressed proteins)
umap<- RunUMAP(seur_merged, assay = "ADT", features = features.use, reduction.key = "adtUMAP_", slot="scale.data")
seur_merged [["adtUMAP"]] <- CreateDimReducObject(embeddings = Embeddings(umap, reduction = "umap"), key = "adtUMAP_", assay ="ADT")
rm(umap)


#### Visualise the ADT UMAP plot, coloured by sample to check for batch effects
DimPlot(seur_merged, reduction = "adtUMAP",group.by = "sample")+
  ggtitle(paste("Samples: ADT UMAP"))+
  xlab("UMAP_1")+ylab("UMAP_2")

#### Visualise the ADT UMAP plot, coloured by RNA clusters
DimPlot(seur_merged, reduction = "adtUMAP",group.by = "cluster", label=T, repel=T)+
    ggtitle(paste("RNA clusters: ADT UMAP"))+
    xlab("UMAP_1")+ylab("UMAP_2")

#### Louvain clustering based on the ADT distance matrix
adt.data <- GetAssayData(seur_merged, slot = "data", assay="ADT")[features.use,]
adt.dist <- dist(t(as.matrix(adt.data)))
seur_merged[["snn_adt"]] <- FindNeighbors(adt.dist)$snn
for ( i in seq(0,1.5, 0.25)){  # Run FindClusters for a range of resolutions
  seur_merged <- FindClusters(seur_merged, resolution = i, graph.name="snn_adt")
}

#### Determine the final resolution of clustering, based on the stability of clustering, visualised in the clustering tree plot
clustree(seur_merged, prefix = "snn_adt_res.")+
  ggtitle("ADT")

#### Visualise the ADT UMAP plot, coloured by ADT clusters
res=0.5
DimPlot(seur_merged, reduction = "adtUMAP",group.by = paste0("snn_adt_res.",res), label=T, repel=T)+
      ggtitle(paste("ADT clusters: ADT UMAP"))+
      xlab("UMAP_1")+ylab("UMAP_2")

########## 3.Harmony batch correction -if needed, apply batch correction

#### Select a value for theta - diversity clustering penalty parameter.
#### Default theta=2. Larger values of theta result in stronger integration, but can lead to over-correction)
theta.use<-0

#### Run Harmony on the normalized ADT count matrix
harmony_embeddings <- HarmonyMatrix(t(GetAssayData(seur_merged, assay = "ADT",slot="data")),
                      meta_data=seur_merged$sample,  do_pca=FALSE,theta = theta.use)

#### Run UMAP on the harmony-corrected count matrix
umap<-CreateSeuratObject(counts=GetAssayData(seur_merged, assay = "ADT",slot="counts"), assay="ADT")
umap<-SetAssayData(object =umap,slot="data",new.data =t(harmony_embeddings), assay="ADT")
umap<-ScaleData(umap, assay="ADT")
umap<- RunUMAP(umap, assay = "ADT", features = features.use, reduction.key = paste0("adtHarmonyTheta",theta.use,"Umap_"), slot="scale.data")
seur_merged[[paste0("adtHarmonyTheta",theta.use,"Umap")]] <- CreateDimReducObject(embeddings = Embeddings(umap, reduction = "umap"),
                                                                           key = paste0("adtHarmonyTheta",theta.use,"Umap_"),
                                                                           assay ="ADT")
rm(umap)

#### Visualise the ADT UMAP plot, coloured by sample to check if the  batch effects were resolved
DimPlot(seur_merged, reduction = paste0("adtHarmonyTheta",theta.use,"Umap"),group.by = "sample")+
  ggtitle(paste("Samples: ADT UMAP Harmony theta=",theta.use))+
  xlab("UMAP_1")+ylab("UMAP_2")

#### Visualise the ADT UMAP plot, coloured by clusters
DimPlot(seur_merged, reduction =  paste0("adtHarmonyTheta",theta.use,"Umap"),group.by = "cluster",
        label=T, repel=T)+
    ggtitle(paste("RNA clusters: ADT UMAP Harmony theta=",theta.use))+
    xlab("UMAP_1")+ylab("UMAP_2")

#### Louvain clustering based on the harmony-corrected ADT distance matrix
adt.data.harmony <- harmony_embeddings [,colnames(harmony_embeddings) %in% features.use]
adt.dist.harmony <- dist(as.matrix(adt.data.harmony))
seur_merged[[paste0("snn_adt_harmony_theta",theta.use)]] <- FindNeighbors(adt.dist.harmony)$snn
for ( i in seq(0,1.5, 0.25)){
  seur_merged <- FindClusters(seur_merged, resolution = i, graph.name=paste0("snn_adt_harmony_theta",theta.use))
}

#### Determine the final resolution of clustering, based on the stability of clustering, visualised in the clustering tree plot
clustree(seur_merged, prefix = paste0("snn_adt_harmony_theta",theta.use,"_res."))+
  ggtitle(paste("ADT: Harmony theta =",theta.use))

#### Visualise the ADT UMAP plot, coloured by ADT clusters
res=0.5
DimPlot(seur_merged, reduction =  paste0("adtHarmonyTheta",theta.use,"Umap"),
        group.by = paste0("snn_adt_harmony_theta",theta.use,"_res.",res), label=T, repel=T)+
    ggtitle(paste0("ADT clusters: ADT UMAP Harmony theta=",theta.use))+
    xlab("UMAP_1")+ylab("UMAP_2")

sessionInfo()
'''
R version 4.0.2 (2020-06-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
LAPACK: /home/daliya/anaconda3/lib/libmkl_rt.so

locale:
  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=de_BE.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=de_BE.UTF-8    LC_MESSAGES=en_US.UTF-8
[7] LC_PAPER=de_BE.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=de_BE.UTF-8 LC_IDENTIFICATION=C

attached base packages:
  [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
  [1] scater_1.16.2               cowplot_1.0.0               ggrepel_0.8.2               clustree_0.4.3              ggraph_2.0.3                SingleR_1.2.4
[7] VennDiagram_1.6.20          futile.logger_1.4.3         DropletUtils_1.8.0          SingleCellExperiment_1.10.1 SummarizedExperiment_1.18.2 DelayedArray_0.14.1
[13] matrixStats_0.56.0          Biobase_2.48.0              GenomicRanges_1.40.0        GenomeInfoDb_1.24.2         IRanges_2.22.2              S4Vectors_0.26.1
[19] BiocGenerics_0.34.0         Matrix_1.2-18               dplyr_1.0.0                 Seurat_3.2.2                ggplot2_3.3.2

loaded via a namespace (and not attached):
  [1] utf8_1.1.4                    reticulate_1.16               R.utils_2.9.2                 tidyselect_1.1.0              RSQLite_2.2.0
[6] AnnotationDbi_1.50.3          htmlwidgets_1.5.1             BiocParallel_1.22.0           Rtsne_0.15                    munsell_0.5.0
[11] codetools_0.2-16              ica_1.0-2                     future_1.18.0                 miniUI_0.1.1.1                withr_2.2.0
[16] colorspace_1.4-1              knitr_1.29                    rstudioapi_0.11               ROCR_1.0-11                   tensor_1.5
[21] listenv_0.8.0                 labeling_0.3                  GenomeInfoDbData_1.2.3        polyclip_1.10-0               bit64_4.0.5
[26] farver_2.0.3                  rhdf5_2.32.2                  vctrs_0.3.2                   generics_0.0.2                lambda.r_1.2.4
[31] xfun_0.16                     BiocFileCache_1.12.0          R6_2.4.1                      ggbeeswarm_0.6.0              graphlayouts_0.7.0
[36] rsvd_1.0.3                    locfit_1.5-9.4                bitops_1.0-6                  spatstat.utils_1.17-0         assertthat_0.2.1
[41] promises_1.1.1                scales_1.1.1                  beeswarm_0.2.3                gtable_0.3.0                  globals_0.12.5
[46] goftest_1.2-2                 tidygraph_1.2.0               rlang_0.4.7                   splines_4.0.2                 textTinyR_1.1.3
[51] lazyeval_0.2.2                BiocManager_1.30.10           yaml_2.2.1                    reshape2_1.4.4                abind_1.4-5
[56] httpuv_1.5.4                  tools_4.0.2                   ellipsis_0.3.1                RColorBrewer_1.1-2            ggridges_0.5.2
[61] Rcpp_1.0.5                    plyr_1.8.6                    base64enc_0.1-3               zlibbioc_1.34.0               purrr_0.3.4
[66] RCurl_1.98-1.2                rpart_4.1-15                  deldir_0.1-28                 pbapply_1.4-3                 viridis_0.5.1
[71] zoo_1.8-8                     cluster_2.1.0                 magrittr_1.5                  data.table_1.13.0             futile.options_1.0.1
[76] lmtest_0.9-37                 RANN_2.6.1                    fitdistrplus_1.1-1            patchwork_1.0.1.9000          mime_0.9
[81] evaluate_0.14                 xtable_1.8-4                  gridExtra_2.3                 compiler_4.0.2                tibble_3.0.3
[86] KernSmooth_2.23-17            crayon_1.3.4                  R.oo_1.24.0                   htmltools_0.5.0               mgcv_1.8-33
[91] later_1.1.0.1                 tidyr_1.1.0                   DBI_1.1.0                     tweenr_1.0.1                  formatR_1.7
[96] ExperimentHub_1.14.0          dbplyr_1.4.4                  MASS_7.3-53                   rappdirs_0.3.1                cli_2.0.2
[101] R.methodsS3_1.8.1             igraph_1.2.5                  pkgconfig_2.0.3               plotly_4.9.2.1                vipor_0.4.5
[106] dqrng_0.2.1                   XVector_0.28.0                stringr_1.4.0                 digest_0.6.25                 sctransform_0.3
[111] RcppAnnoy_0.0.16              spatstat.data_1.4-3           rmarkdown_2.3                 leiden_0.3.3                  uwot_0.1.8
[116] edgeR_3.30.3                  DelayedMatrixStats_1.10.1     curl_4.3                      shiny_1.5.0                   lifecycle_0.2.0
[121] nlme_3.1-149                  jsonlite_1.7.0                Rhdf5lib_1.10.1               BiocNeighbors_1.6.0           viridisLite_0.3.0
[126] limma_3.44.3                  fansi_0.4.1                   pillar_1.4.6                  lattice_0.20-41               fastmap_1.0.1
[131] httr_1.4.2                    survival_3.2-7                interactiveDisplayBase_1.26.3 glue_1.4.2                    spatstat_1.64-1
[136] png_0.1-7                     BiocVersion_3.11.1            bit_4.0.4                     ggforce_0.3.2                 stringi_1.4.6
[141] HDF5Array_1.16.1              blob_1.2.1                    BiocSingular_1.4.0            AnnotationHub_2.20.0          memoise_1.1.0
[146] irlba_2.3.3                   future.apply_1.6.0

'''
