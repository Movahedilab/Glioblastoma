suppressPackageStartupMessages({
 library(Seurat)
 library(dplyr)
 library(VennDiagram)
 library(harmony)
 library(clustree)
})

path_data<-"~/path/to/data/"
path_output<-"~/path/to/output_folder/"

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
ADT_raw <- Read10X(paste0(path_data,"raw_feature_bc_matrix/"))[["Antibody Capture"]]


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

##### Save object
saveRDS(seur_merged, file=paste0(path_output,"Robjects/ADT.RNA.seuratObj.rds"))
