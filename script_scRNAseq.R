suppressPackageStartupMessages({
  library(Seurat)
  library(scater)
  library(dplyr)
  library(mvoutlier)
  library(VennDiagram)
  library(harmony)
  library(clustree)
})


path_data<-"~/path/to/data/"
path_output<-"~/path/to/output_folder/"


################################################################################
########## QC: CELLS
################################################################################


##########  1 .Read the raw filtered RNA count data ##########

expr.mat <-Read10X(paste0(path_data,"filtered_feature_bc_matrix/"))
## For CITE-seq data:
# expr.mat <-Read10X("outs/filtered_feature_bc_matrix/")[["Gene Expression"]]

##### Convert to SingleCellExperiment object
sce<-SingleCellExperiment(assays=list(counts=expr.mat))
rm(expr.mat)

## When loaded from a Cellranger aggregated matrix, find which cells come from which sample:
#sce$sample<-sapply(strsplit(colnames(sce), "-"), "[[",2)
#table(sce$sample)


##########  2. Calculate QC metrics ##########  

##### Get mitochondrial genes
mito.genes <- grep("^mt-", rownames(sce))
## For human data:
# mito.genes <- grep("^MT-", rownames(sce))
length(mito.genes)
# 13

##### Calculate QC metrics per cell
sce$percent.mito <- (Matrix::colSums(counts(sce)[mito.genes, ])*100)/Matrix::colSums(counts(sce))  ## percentage expressed mitochondrial genes 
sce$nGene<-apply(counts(sce),  2,  function(x) length(x[x > 0])) # number of expressed genes
sce$nUMI<-apply(counts(sce),  2,  sum) # total UMI counts (library size)
sce$staticNr<-1
dim(colData(sce))
colnames(colData(sce))


########## 3. Get outliers ##########
## Outliers are determined based on the median absolute deviation (MAD). It is possible to modify the nmads parameter 
## (minimum number of MADs away from median required for a value to be called an outlier)

##### Aim: remove outlier cells for low library size, low number of expressed genes and high mitochondrial gene proportions
## UMI counts per cell
sce$nUMI.outlier.low <- isOutlier(sce$nUMI, nmads=3, type="lower", log=TRUE)
sum(sce$nUMI.outlier.low )

## Number of expressed genes
sce$nGene.outlier.low <- isOutlier(sce$nGene, nmads=3, type="lower", log=TRUE)
sum(sce$nGene.outlier.low)

## Mitochondrial gene proportions 
sce$mito.outlier.high <- isOutlier(sce$percent.mito, nmads=3, type="higher", log=TRUE)
sum(sce$mito.outlier.high)


##### Create histograms
hist(sce$nUMI,
     breaks = 100,
     main="Histogram of total UMI counts per cell")
  abline(v = max(sce$nUMI[sce$nUMI.outlier.low]), col = "red")
hist(sce$nGene,
     breaks = 100, 
     main="Histogram of number of genes per cell")
  abline(v = max(sce$nGene[sce$nGene.outlier.low], na.rm=T), col = "red")
hist(sce$percent.mito,
     breaks = 100,
     main="Histogram of % mito genes per cell")
  abline(v = min(sce$percent.mito[sce$mito.outlier.high]), col = "red")


###### Create violin plots 

### Before filtering
metaData=as.data.frame(colData(sce))
ggplot(metaData, aes(staticNr, nUMI)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=nUMI.outlier.low)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle("Total UMI counts per cell")
ggplot(metaData, aes(staticNr, nGene)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=nGene.outlier.low)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle("Number of genes per cell")
ggplot(metaData, aes(staticNr, percent.mito)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=mito.outlier.high)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle("% mito genes per cell")

### After filtering
metaData.filtered<-metaData[! (metaData$nUMI.outlier.low | metaData$nGene.outlier.low | metaData$mito.outlier.high) ,]
ggplot(metaData.filtered, aes(staticNr, nUMI)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=nUMI.outlier.low)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle("Total UMI counts per cell after filtering")
ggplot(metaData.filtered, aes(staticNr, nGene)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=nGene.outlier.low)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle("Number of genes per cell after filtering")
ggplot(metaData.filtered, aes(staticNr, percent.mito)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=mito.outlier.high)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle("% mito genes per cell after filtering")


### Create venn diagram
v <-venn.diagram(
  list (UMI=rownames(colData(sce)[sce$nUMI.outlier.low,]),
        gene=rownames(colData(sce)[sce$nGene.outlier.low,]),
        mito=rownames(colData(sce)[sce$mito.outlier.high,])),
  filename=NULL,
  alpha = c( 0.5,0.5,0.5),
  fill = c("green","orange","blue")
)
grid.newpage()
grid.draw(v)
rm(v)


###### Remove outliers
sce_clean <- sce[,!(sce$nUMI.outlier.low | sce$nGene.outlier.low | sce$mito.outlier.high)]
dim(sce_clean)

### Number of cells removed
ncol(sce_clean)-ncol(sce)



########## 4. Get PCA outliers ########## 
### Get additional outliers, using multivariate outlier detection function of scater, 
### based on PCA computed from QC metric data (uses the mvoutlier package)

### Calculate additional QC metrics with scater
is.mito <-grepl("^mt-",rownames(sce_clean))
sum(is.mito)
sce_clean <-calculateQCMetrics(sce_clean, feature_controls=list(Mt=is.mito))
### In scater v. 1.16.2
### sce_clean <- addPerCellQC(sce_clean, subsets=list(Mt=is.mito), flatten=T)

### Choose QC variables for authomatic outlier detection
selected_variables <- c("pct_counts_in_top_100_features", 
                        "total_features_by_counts", "pct_counts_feature_control", 
                        "total_features_by_counts_feature_control", "log10_total_counts_endogenous", 
                        "log10_total_counts_feature_control")
setdiff(selected_variables, colnames(colData(sce_clean)))
### Detect outliers
sce_clean<-runPCA(sce_clean, use_coldata = TRUE , detect_outliers=TRUE,variables=selected_variables)
### In scater v. 1.16.2
### sce_clean=runColDataPCA(sce_clean, outliers=T, variables=selected_variables)
table(sce_clean$outlier)

###### Create scatter plots 
ggplot(as.data.frame(colData(sce_clean)), aes(nUMI, percent.mito)) + 
  geom_point(fill="gray80", aes(col=outlier)) + 
  scale_color_manual(values=c("#00BFC4", "#F8766D"))

ggplot(as.data.frame(colData(sce_clean)), aes(nUMI, nGene)) + 
  geom_point(fill="gray80", aes(col=outlier)) + 
  scale_color_manual(values=c("#00BFC4", "#F8766D"))

###### Create violin plots 
ggplot(as.data.frame(colData(sce_clean)), aes(staticNr, nUMI)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=outlier)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle("Total UMI counts per cell")
ggplot(as.data.frame(colData(sce_clean)), aes(staticNr, total_features_by_counts)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=outlier)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle("Number of genes per cell")
ggplot(as.data.frame(colData(sce_clean)), aes(staticNr, percent.mito)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=outlier)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle("% mito genes per cell")

##### Color outlier cells on PCA plot calculated on QC metrics 
plotReducedDim(sce_clean, use_dimred = "PCA_coldata", colour_by='outlier',shape_by='outlier') + labs(title="PCA with outliers colored")


##### Remove outlier cells ####
sce_clean <- sce_clean[,!sce_clean$outlier]
dim(sce_clean)

################################################################################
########## QC: GENES
################################################################################
### To define a cutoff of lowly-abundant genes, we plot the distribution of log-means across all genes. 
### The cutoff is placed in middle of the rectangular component of the graph before the peak
ave.counts <- rowMeans(as.matrix(counts(sce_clean)))

thresh<-0.005  #cutoff can be modified here
hist(log10(ave.counts), breaks=100, col="grey80",
     xlab=expression(Log[10]~"mean count"))
abline(v=log10(thresh), col="blue", lwd=2, lty=2)

rowData(sce_clean)$usegenes<-ave.counts>thresh
table(rowData(sce_clean)$usegenes)

### Filter out the lowly-abundant genes
sce_clean<-sce_clean[rowData(sce_clean)$usegenes, ]

################################################################################
########## CONVERT TO SEURAT OBJECT
################################################################################
counts<-counts(sce_clean)
rownames(counts)<-rownames(sce_clean)
colnames(counts)<-colnames(sce_clean)
seuratObj <- CreateSeuratObject(counts =counts )
meta.data.keep<-colData(sce_clean)
meta.data.keep<-meta.data.keep[,c("percent.mito","sample")]
seuratObj<-AddMetaData(seuratObj, as.data.frame(meta.data.keep))

################################################################################
##########  MERGE SEURAT OBJECT FROM SEVERAL SAMPLES (IF NEEDED)
################################################################################

seuratObj$sample="sample1"
seuratObj2=readRDS("path/to/sample2_SeuratObject.rds")
seuratObj2$sample="sample2"
seuratObj<-merge(seuratObj,seuratObj22)

################################################################################
########## NORMALIZATION, HVG DETECTION and PCA
################################################################################

##### Normalize data
seuratObj <- NormalizeData(seuratObj,verbose = F)

##### HVG detection 
seuratObj <- FindVariableFeatures(seuratObj,verbose=F)

##### Scale data per gene
seuratObj <- ScaleData(seuratObj,verbose=F)

##### PCA
seuratObj <- RunPCA(seuratObj, features =VariableFeatures(seuratObj))

#### Select PCs for downstream analysis of the dataset
## Seurat provides a heuristic method to help us select PC components. It generates a ranking of PCs based on 
## the percentage of variance explained by each of them 
ElbowPlot(object = seuratObj,ndims =50)

########## Heatmap of the genes that drive each PC
DimHeatmap(seuratObj, dims = 1:30, cells = 5000, balanced = TRUE)

########## Determine statistically significant PCs - time consuming step
# seuratObj <- JackStraw(seuratObj, num.replicate = 100, dims=30)
# JackStrawPlot(seuratObj, dims = 1:30)


################################################################################
########## CLUSTERING AND UMAP DIMENSIONALITY REDUCTION 
################################################################################

##Final PC selection
dims.use<-30

### Run UMAP non-linear dimensional reduction for visualisation of the data
seuratObj <- RunUMAP(seuratObj, dims = 1:dims.use, verbose=F)

#### Visualise the UMAP plot, coloured by sample to check for batch effects
DimPlot(object = seuratObj, group.by="sample")

#### Louvain clustering
seuratObj <- FindNeighbors(seuratObj, dims = 1:dims.use, verbose=F)
for ( i in seq(0,2, 0.25))  # Run FindClusters for a range of resolutions
  seuratObj <- FindClusters(seuratObj, resolution = i, verbose=F)

#### Plot of a clustering tree showing the relationship between clusterings at different resolutions (using the clustree package)
### Each cluster forms a node in the tree and edges are constructed by considering how the cells move from lower to higher resolution
clustree(seuratObj, prefix = "RNA_snn_res.")+
  ggtitle("Clustering tree")

#### Visualise the clustering for selected resolutions on a UMAP plot
plot<-list()
for ( res in c(0.5, 0.75,1,1.25))
  plot[[as.character(res)]]<-DimPlot(seuratObj, pt.size = 1,label=T,repel=T, group.by = paste0("RNA_snn_res.",res)) +
  ggtitle(paste("res=",res))
plot_grid(plotlist=plot)



################################################################################
########## HARMONY BATCH CORRECTION (IF NEEDED)
################################################################################
#### Select a value for theta - diversity clustering penalty parameter.
#### Default theta=2. Larger values of theta result in stronger integration, but can lead to over-correction)
theta.use<-1
seuratObj<-RunHarmony(seuratObj, group.by.vars="sample",theta =theta.use,  
                      plot_convergence = TRUE, reduction.save =paste0("harmony_theta",theta.use), verbose=F) 

#### Run UMAP on the harmony-corrected PCA
seuratObj <- RunUMAP(seuratObj,reduction =paste0("harmony_theta",theta.use), dims = 1:dims.use,verbose =F, 
                reduction.name = paste0("umapHarmonyTheta",theta.use, "PC",dims.use),
                reduction.key = paste0("umapHarmonyTheta",theta.use, "PC",dims.use,"_"))  

#### Visualise the harmony-corrected UMAP plot, coloured by sample to check if the  batch effects were resolved
DimPlot(seuratObj, reduction = paste0("umapHarmonyTheta",theta.use, "PC",dims.use),group.by = "sample")+
  ggtitle(paste("Samples: UMAP Harmony theta=",theta.use, " PC=",dims.use))+
  xlab("UMAP_1")+ylab("UMAP_2")

#### Louvain clustering on the harmony-corrected PCA
seuratObj <- FindNeighbors(seuratObj, dims = 1:dims.use, reduction =paste0("harmony_theta",theta.use), 
                      graph.name = paste0("RNA_snn_harmony_theta",theta.use, ".PC",dims.use), verbose=F)
for ( i in seq(0,2, 0.25))
  seuratObj <- FindClusters(seuratObj, resolution = i, graph.name=paste0("RNA_snn_harmony_theta",theta.use, ".PC",dims.use), verbose=F)

#### Plot of a clustering tree 
clustree(seuratObj, prefix = paste0("RNA_snn_harmony_theta",theta.use, ".PC",dims.use,"_res."))+
  ggtitle(paste("Harmony theta =",theta.use,"; PC =",dims.use))

plot<-list()
for ( res in c(0.25,0.5, 0.75,1))
  plot[[as.character(res)]]<-DimPlot(seur, pt.size = 1,label=T,repel=T, group.by = paste0("RNA_snn_harmony_theta",theta.use, ".PC",dims.use,"_res.",res),
                                     reduction=paste0("umapHarmonyTheta",theta.use, "PC",dims.use)) +
  ggtitle(paste("Theta =",theta.use,"PC =",dims.use,"res=",res))
plot_grid(plotlist=plot)



##### Save object
saveRDS(seuratObj, file=paste0(path_output,"Robjects/seuratObj.rds"))


