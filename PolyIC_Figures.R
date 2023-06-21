########################################
##                                    ##
##          Install Packages          ##
##                                    ##
########################################

#install.packages('Seurat')

# Needs python
#library(reticulate)
#py_config()
# create a new environment 
#conda_create("r-reticulate")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("biomaRt")
#BiocManager::install(c("edgeR", "DESeq2", "limma"))

#devtools::install_github("neurorestore/Libra")

#install.packages("glmmTMB")
#install.packages('TMB', type = 'source')
#BiocManager::install("EnrichmentBrowser")

#install.packages("cli")
#install.packages('rlang')

### You will need this to update R ###
#install.packages("installr")
#install.packages("stringr")
#library(installr)
#install.packages("dplyr")
#install.packages("cowplot")
#updateR()
#install.packages("harmony")
#install.packages('Rcpp')

########################################
##                                    ##
##               Prep                 ##
##                                    ##
########################################

# Set a working Directory

setwd("XXXX/PolyICFigures")
getwd()

#Example
#setwd("C:/Users/adam/Documents/PolyICFigures")


# Load packages
library(Seurat)
library(cowplot)
library(dplyr)
library(sctransform)
library(ggplot2)
library(Matrix)
library(Rcpp)
library(patchwork)
library(EnrichmentBrowser)


########################################
##                                    ##
##        Load Data and Labels        ##
##                                    ##
########################################

# This R object contains all data in a Seurat object
# Contact me if you'd like to discuss ways to share this object
# Otherwise, it can be recreated using the BAM files uploaded to SRA.
# But I do not have any code right now to do this, because I've never had to. 
# I am working on it.

load("XXX/PolyICFigures/PBMC_20220711.Robj")

# Example
#load("C:/Users/adam/Documents/PolyICFigures/PBMC_20220711.Robj")

# I actually can't remember if I need all of these metadata labels, but they might be useful
pbmc.integrated$disease_weight <- paste(pbmc.integrated$disease,
                                        pbmc.integrated$weight, sep = "_")
pbmc.integrated$cell_disease_treatments <- paste(pbmc.integrated$celltype, 
                                                 pbmc.integrated$disease,
                                                 pbmc.integrated$first, 
                                                 pbmc.integrated$second, sep = "_")
pbmc.integrated$cell_disease_weight_treatments <- paste(pbmc.integrated$celltype, 
                                                        pbmc.integrated$disease,
                                                        pbmc.integrated$weight,
                                                        pbmc.integrated$first, 
                                                        pbmc.integrated$second, sep = "_")
pbmc.integrated$patient_cell_disease_weight_treatments <- paste(pbmc.integrated$celltype,
                                                                pbmc.integrated$disease,
                                                                pbmc.integrated$weight,
                                                                pbmc.integrated$orig.ident,
                                                                pbmc.integrated$first, 
                                                                pbmc.integrated$second, sep = "_")
pbmc.integrated$cell_disease <- paste(pbmc.integrated$celltype, 
                                      pbmc.integrated$disease, sep = "_")
pbmc.integrated$cell_disease_weight <- paste(pbmc.integrated$celltype, 
                                             pbmc.integrated$disease, 
                                             pbmc.integrated$weight,
                                             sep = "_")
pbmc.integrated$treatments <- paste(pbmc.integrated$first, 
                                    pbmc.integrated$second, sep = "_")
pbmc.integrated$weight_treatments <- paste(pbmc.integrated$weight,
                                           pbmc.integrated$first, 
                                           pbmc.integrated$second, sep = "_")


table(Idents(pbmc.integrated))
Idents(pbmc.integrated) <- pbmc.integrated$celltype

########################################
##                                    ##
##                UMAPS               ##
##                                    ##
########################################

# An option for Supplemental Figure

pdf("UMAP_allCond.pdf", height = 10, width = 12, useDingbats=FALSE)
DimPlot(pbmc.integrated, reduction = 'umap', split.by = 'ultra', ncol = 4)
dev.off()

# Figure 2A
pdf("UMAP_all.pdf", height = 10, width = 12, useDingbats=FALSE)
DimPlot(pbmc.integrated, reduction = 'umap', label = TRUE, label.size = 8)
dev.off()

DimPlot(pbmc.integrated, reduction = 'umap', split.by = 'ultra', ncol = 4)

Idents(Mono.integrated) <- Mono.integrated$seurat_clusters


########################################
##                                    ##
##     Load BioMart for GO terms      ##
##                                    ##
########################################

# This code is not used, but can be helpful for illustrative purposes

# Original code came from:
# https://www.biostars.org/p/52101/
library(GO.db)
library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations

# How to use this:

# If you have a pathway you are interested in, and a specific GO term for it
# You can imput in that GO terms, and this code will get all gene names for that GO term,
# and even everything withing subcategories associated.

# Using biomart, pull in all genes included within a specified GO term

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations

# get gene symbol, transcript_id and go_id for all genes annotated with  a Go term
# EXAMPLE: GO:0042110 - Tcell Activation
gene.data <- getBM(attributes=c('hgnc_symbol'),
                   filters = 'go', values = 'GO:0034340', mart = ensembl)
x <- c(gene.data$hgnc_symbol) # x is now a vector with all the gene names

# Then get the percentage of counts per cell associated with those genes
# This is similar to the percentage mitochondria code commonly used
# This creates a meta label with the percentage for that gene set

# Continuing the example with that tcell GO term
DefaultAssay(pbmc.integrated) <- "RNA"
pbmc.integrated <- PercentageFeatureSet(pbmc.integrated, 
                                             features = x, 
                                             col.name = "tcell", 
                                             assay = "RNA")
# This will recompute the total counts associated with those genes
pbmc.integrated_temp$interferon_tot <- pbmc.integrated_temp$tcell*pbmc.integrated_temp$nCount_RNA


# Make a violin plot
DefaultAssay(pbmc.integrated) <- "SCT"
VlnPlot(pbmc.integrated, 
        features = c("tcell"), 
        split.by = "ultra", 
        idents = c("CD4-T-cell1", "CD4-T-cell2", 
                   "CD4-T-cell3", "CD4-T-cell4", 
                   "CD4-T-cell5", "CD4-T-cell6", 
                   "CD8-T-cell1", "CD8-T-cell2", 
                   "B-cell",
                   "CD14_Monocyte1", "CD14_Monocyte2", "CD16_Monocyte",
                   "NK-cell",
                   "DC"),
        pt.size = 0, ncol = 1)

# NOTE - SOMETIMES THIS DOESN'T WORK
# And I've never figured it out. It's usually due to a single gene in the list
# For the IFNG pathways I use in the paper, I had to break the list into smaller and smaller groups
# to narrow down the bad gene. And yes, it was a single gene that messed the entire code up

########################################
##                                    ##
##      Number of cells/cluster       ##
##                                    ##
########################################

DefaultAssay(pbmc.integrated) <- "SCT"

#Calculate number of cells per cluster from object@ident
table(Idents(pbmc.integrated))
table(pbmc.integrated$orig.ident)

Idents(pbmc.integrated) <- pbmc.integrated$cell_disease_treatments
table(Idents(pbmc.integrated))

cell.num <- as.data.frame(table(Idents(pbmc.integrated)))
write.table(cell.num, file="CellTypes_PBMCscRNA_v2.txt")

cell.num <- as.data.frame(table(pbmc.integrated$patient_cell_disease_weight_treatments))
write.table(cell.num, file="CellTypes_AddedLabels_PBMCscRNA.txt")

########################################
##                                    ##
##            Violin Plots            ##
##                                    ##
########################################

# All violin plots in the manuscript were made using this code.
# This is fairly adjustable, just replace gene names and the list of cell clusters you want to see
# The confusing part is the metadata, which cells do you want to see and how is it split
# In most cases, you want to split the data by the metadata label, "treatments"
# This includes Basal, PolyIC, LPS, and LPS_PolyIC
# Note about the metadata:
# All drug challenges are labeled under First or Second
# First is either Basal or PolyIC
# Second is either Basal or LPS
# So from there you can get combinations
# The metadata label "treatments" is always a combination of First_Second

# Split by treatments (PolyIC and LPS)
DefaultAssay(pbmc.integrated) <- "SCT"
table(Idents(pbmc.integrated))
Idents(pbmc.integrated)<-pbmc.integrated$cell_disease #format: "CD14_Monocyte2_HD"

# For all cell,s very simple, but messy
VlnPlot(pbmc.integrated, features = c("IL1B"), 
                 pt.size = 0, combine = FALSE
                 )

# Reorder the labels for figures
Idents(pbmc.integrated)<-factor(Idents(pbmc.integrated), levels = c("CD4-T-cell1_HC","CD4-T-cell1_HD", 
                                                                    "CD4-T-cell2_HC","CD4-T-cell2_HD", 
                                                                    "CD4-T-cell3_HC","CD4-T-cell3_HD", 
                                                                    "CD4-T-cell4_HC","CD4-T-cell4_HD", 
                                                                    "CD4-T-cell5_HC","CD4-T-cell5_HD", 
                                                                    "CD4-T-cell6_HC","CD4-T-cell6_HD", 
                                                                    "CD8-T-cell1_HC","CD8-T-cell1_HD", 
                                                                    "CD8-T-cell2_HC","CD8-T-cell2_HD", 
                                                                    "B-cell_HC","B-cell_HD",
                                                                    "CD14_Monocyte1_HC","CD14_Monocyte1_HD",
                                                                    "CD14_Monocyte2_HC","CD14_Monocyte2_HD",
                                                                    "CD16_Monocyte_HC","CD16_Monocyte_HD",
                                                                    "NK-cell_HC","NK-cell_HD",
                                                                    "DC_HC","DC_HD"))

pdf("ViolinPlot_Treatments.pdf", height = 6, width = 12, useDingbats=FALSE)
plots <- VlnPlot(pbmc.integrated, features = c("IFNGR1", "IFNGR2"), split.by = "treatments",
                 pt.size = 0, combine = FALSE,
                 idents = c("CD14_Monocyte1_HC","CD14_Monocyte1_HD",
                            "CD14_Monocyte2_HC","CD14_Monocyte2_HD",
                            "CD16_Monocyte_HC","CD16_Monocyte_HD"
                            ))
plots <- lapply(X = plots,
                FUN = function(p) p 
                + ggplot2::theme(axis.title.y = element_blank())
                + ggplot2::scale_fill_manual(values = c('black','red','grey','pink'))
                )
CombinePlots(plots = plots, ncol = 1)
dev.off()


# Split by treatments (PolyIC and LPS) and weight (Lean and Obese)

# First way: (I think I prefer this way)
DefaultAssay(pbmc.integrated) <- "SCT"
table(Idents(pbmc.integrated))
Idents(pbmc.integrated)<-pbmc.integrated$cell_disease #format: "CD14_Monocyte2_HD"

# Reorder the labels for figures - Two ways, user choice
Idents(pbmc.integrated)<-factor(Idents(pbmc.integrated), levels = c("CD4-T-cell1_HC","CD4-T-cell1_HD", 
                                                                    "CD4-T-cell2_HC","CD4-T-cell2_HD", 
                                                                    "CD4-T-cell3_HC","CD4-T-cell3_HD", 
                                                                    "CD4-T-cell4_HC","CD4-T-cell4_HD", 
                                                                    "CD4-T-cell5_HC","CD4-T-cell5_HD", 
                                                                    "CD4-T-cell6_HC","CD4-T-cell6_HD", 
                                                                    "CD8-T-cell1_HC","CD8-T-cell1_HD", 
                                                                    "CD8-T-cell2_HC","CD8-T-cell2_HD", 
                                                                    "B-cell_HC","B-cell_HD",
                                                                    "CD14_Monocyte1_HC","CD14_Monocyte1_HD",
                                                                    "CD14_Monocyte2_HC","CD14_Monocyte2_HD",
                                                                    "CD16_Monocyte_HC","CD16_Monocyte_HD",
                                                                    "NK-cell_HC","NK-cell_HD",
                                                                    "DC_HC","DC_HD"))
Idents(pbmc.integrated)<-factor(Idents(pbmc.integrated), levels = c("CD4-T-cell1_HC","CD4-T-cell2_HC","CD4-T-cell3_HC",
                                                                    "CD4-T-cell4_HC","CD4-T-cell5_HC","CD4-T-cell6_HC",
                                                                    "CD8-T-cell1_HC","CD8-T-cell2_HC","B-cell_HC",
                                                                    "CD14_Monocyte1_HC", "CD14_Monocyte2_HC","CD16_Monocyte_HC",
                                                                    "NK-cell_HC","DC_HC",
                                                                    "CD4-T-cell1_HD","CD4-T-cell2_HD","CD4-T-cell3_HD",
                                                                    "CD4-T-cell4_HD","CD4-T-cell5_HD","CD4-T-cell6_HD", 
                                                                    "CD8-T-cell1_HD","CD8-T-cell2_HD","B-cell_HD",
                                                                    "CD14_Monocyte1_HD","CD14_Monocyte2_HD","CD16_Monocyte_HD",
                                                                    "NK-cell_HD","DC_HD"))

pdf("ViolinPlot_WeightTreatments_IFNGR.pdf", height = 6, width = 12, useDingbats=FALSE)
plots <- VlnPlot(pbmc.integrated, features = c("IFNGR1", "IFNGR2"), split.by = "weight_treatments",
                 pt.size = 0, combine = FALSE,
                 idents = c("CD14_Monocyte1_HC","CD14_Monocyte1_HD",
                            "CD14_Monocyte2_HC","CD14_Monocyte2_HD",
                            "CD16_Monocyte_HC","CD16_Monocyte_HD"
                 ))
plots <- lapply(X = plots,
                FUN = function(p) p 
                + ggplot2::theme(axis.title.y = element_blank())
                + ggplot2::scale_fill_manual(values = c('black','red','grey','pink',
                                                        'black','red','grey','pink'))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()



# Second way:
DefaultAssay(pbmc.integrated) <- "SCT"
table(Idents(pbmc.integrated))
Idents(pbmc.integrated)<-pbmc.integrated$cell_disease_weight #format: "CD14_Monocyte2_HD_Lean"

# Reorder the labels for figures
Idents(pbmc.integrated)<-factor(Idents(pbmc.integrated), 
                                levels = c("CD4-T-cell1_HC_Lean","CD4-T-cell1_HC_Obese","CD4-T-cell1_HD_Lean","CD4-T-cell1_HD_Obese",
                                           "CD4-T-cell2_HC_Lean","CD4-T-cell2_HC_Obese","CD4-T-cell2_HD_Lean","CD4-T-cell2_HD_Obese", 
                                           "CD4-T-cell3_HC_Lean","CD4-T-cell3_HC_Obese","CD4-T-cell3_HD_Lean","CD4-T-cell3_HD_Obese",
                                           "CD4-T-cell4_HC_Lean","CD4-T-cell4_HC_Obese","CD4-T-cell4_HD_Lean","CD4-T-cell4_HD_Obese",
                                           "CD4-T-cell5_HC_Lean","CD4-T-cell5_HC_Obese","CD4-T-cell5_HD_Lean","CD4-T-cell5_HD_Obese",
                                           "CD4-T-cell6_HC_Lean","CD4-T-cell6_HC_Obese","CD4-T-cell6_HD_Lean","CD4-T-cell6_HD_Obese",
                                           "CD8-T-cell1_HC_Lean","CD8-T-cell1_HC_Obese","CD8-T-cell1_HD_Lean","CD8-T-cell1_HD_Obese",
                                           "CD8-T-cell2_HC_Lean","CD8-T-cell2_HC_Obese","CD8-T-cell2_HD_Lean","CD8-T-cell2_HD_Obese",
                                           "B-cell_HC_Lean","B-cell_HC_Obese","B-cell_HD_Lean","B-cell_HD_Obese",
                                           "CD14_Monocyte1_HC_Lean","CD14_Monocyte1_HC_Obese","CD14_Monocyte1_HD_Lean","CD14_Monocyte1_HD_Obese",
                                           "CD14_Monocyte2_HC_Lean","CD14_Monocyte2_HC_Obese","CD14_Monocyte2_HD_Lean","CD14_Monocyte2_HD_Obese",
                                           "CD16_Monocyte_HC_Lean","CD16_Monocyte_HC_Obese","CD16_Monocyte_HD_Lean","CD16_Monocyte_HD_Obese",
                                           "NK-cell_HC_Lean","NK-cell_HC_Obese","NK-cell_HD_Lean","NK-cell_HD_Obese",
                                           "DC_HC_Lean","DC_HC_Obese","DC_HD_Lean","DC_HD_Obese"))

plots <- VlnPlot(pbmc.integrated, features = c("IFNGR1", "IFNGR2"), split.by = "treatments",
                 pt.size = 0, combine = FALSE,
                 idents = c("CD14_Monocyte1_HC_Lean","CD14_Monocyte1_HC_Obese",
                            "CD14_Monocyte1_HD_Lean","CD14_Monocyte1_HD_Obese",
                            "CD14_Monocyte2_HC_Lean","CD14_Monocyte2_HC_Obese",
                            "CD14_Monocyte2_HD_Lean","CD14_Monocyte2_HD_Obese",
                            "CD16_Monocyte_HC_Lean","CD16_Monocyte_HC_Obese",
                            "CD16_Monocyte_HD_Lean","CD16_Monocyte_HD_Obese"
                            ))
plots <- lapply(X = plots,
                FUN = function(p) p 
                + ggplot2::theme(axis.title.y = element_blank()) 
                + ggplot2::scale_fill_manual(values = c('black', 'red','grey','pink'))
                )
CombinePlots(plots = plots, ncol = 1)



########################################
##                                    ##
##             Dot Plots              ##
##        for Monocyte markers        ##
########################################

########################################
#     Cluster Markers for Monocytes    #
########################################
# Subset only the Monocytes and Basal 
# Because the other treatments would influence the baseline markers

# Subset Basal
DefaultAssay(pbmc.integrated) <- "integrated"
Idents(pbmc.integrated) <- pbmc.integrated$treatments # Categorize by treatments
table(Idents(pbmc.integrated))
pbmc.integrated_subset <- subset(pbmc.integrated, idents = c("Basal_Basal")) # Only Basal

#Subset Monocytes
Idents(pbmc.integrated_subset)<-pbmc.integrated_subset$celltype # Categorize by celltype
table(Idents(pbmc.integrated))
MonoSubBasal <- subset(pbmc.integrated_subset, idents = c("CD14_Monocyte1","CD14_Monocyte2","CD16_Monocyte")) # Only monocytes

Idents(MonoSubBasal)  <- MonoSubBasal$celltype
table(Idents(MonoSubBasal))
# find markers for every cluster compared to all remaining cells, report only the positive ones
Mono.markers <- FindAllMarkers(MonoSubBasal, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(Mono.markers, file="ClusterMarkers_Monocytes_Basal_scRNA.txt")

# Dot plot used in Figure 2B
# These genes come from the file we just made of Monocyte Markers
markers.to.plot <- c(
  "H2AC20","CXCL16","FCGR3A",
  "TMSB4X","LST1","ITGAL",
  "PRDX1","ASAH1","TXN","LGALS3","CTSD",
  "CSTB","CD63","SDCBP","CREG1","PLA2G7",
  "MARCKSL1","LUC7L3","STK24","RESF1","ATF7IP",
  "BTG1","THBS1","ELK3","SLC4A7","PABPC1"
)

DefaultAssay(MonoSubBasal) <- "SCT"
pdf("DotPlot_Mono_all.pdf", height = 4, width = 10, useDingbats=FALSE)
DotPlot(MonoSubBasal, features = rev(markers.to.plot), 
        #cols = c("blue","red", "lightblue", "pink"), 
        dot.scale = 8, 
        #split.by = "both",
        assay = 'SCT'
) + RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")
dev.off()


########################################
##                                    ##
##            Split Figures           ##
##              by Treatment          ##
########################################

# For many of the figures in the paper, it is easier to split the Seurat object by treatment
# That way you can look at Basal vs PolyIC with no LPS treated samples in the object
# This also prevents possible down the line confusion, and cell labeling issues.

# POLY IC ONLY
Idents(pbmc.integrated) <- pbmc.integrated$treatments
pbmc.integrated_Poly <- subset(pbmc.integrated, idents = c("Basal_Basal","Poly_Basal"))
table(Idents(pbmc.integrated_Poly))

DefaultAssay(pbmc.integrated_Poly) <- "SCT"
Idents(pbmc.integrated_Poly) <- pbmc.integrated_Poly$cell_disease_weight
table(Idents(pbmc.integrated_Poly))
# NOTE - Re-ordering some of the celltypes, but not all can mess things up, like DotPlots
# So the following line will reorder these cells, but leave the rest unlabels, and ruin some figures
# The better way to do this is to create another subset of just the cells you want
Idents(pbmc.integrated_Poly) <- factor(Idents(pbmc.integrated_Poly),
                                       levels = c("CD14_Monocyte1_HC_Lean","CD14_Monocyte1_HC_Obese",
                                                  "CD14_Monocyte1_HD_Lean","CD14_Monocyte1_HD_Obese",
                                                  "CD14_Monocyte2_HC_Lean","CD14_Monocyte2_HC_Obese",
                                                  "CD14_Monocyte2_HD_Lean","CD14_Monocyte2_HD_Obese",
                                                  "CD16_Monocyte_HC_Lean","CD16_Monocyte_HC_Obese",
                                                  "CD16_Monocyte_HD_Lean","CD16_Monocyte_HD_Obese"
                                      ))

plots <- VlnPlot(pbmc.integrated_Poly, features = c("GBP1","GBP2"), split.by = "first",
                 pt.size = 0, combine = FALSE, 
                 idents = c("CD14_Monocyte1_HC_Lean","CD14_Monocyte1_HC_Obese",
                            "CD14_Monocyte1_HD_Lean","CD14_Monocyte1_HD_Obese"
                            ))
plots <- lapply(X = plots,
                FUN = function(p) p
                + ggplot2::theme(axis.title.y = element_blank())
                + ggplot2::scale_fill_manual(values = c('black', 'red'))
                )
CombinePlots(plots = plots, ncol = 1)


# Dot plot
Idents(pbmc.integrated_Poly)<-pbmc.integrated_Poly$celltype # Categorize by celltype
Mono.integrated_Poly <- subset(pbmc.integrated_Poly, idents = c("CD14_Monocyte1","CD14_Monocyte2")) # only certain monocytes
DefaultAssay(Mono.integrated_Poly) <- "SCT"
Idents(Mono.integrated_Poly)<-Mono.integrated_Poly$cell_disease_weight
Idents(Mono.integrated_Poly) <- factor(Idents(Mono.integrated_Poly),
                                       levels = c("CD14_Monocyte2_HD_Obese","CD14_Monocyte2_HD_Lean",
                                                  "CD14_Monocyte1_HD_Obese","CD14_Monocyte1_HD_Lean",
                                                  "CD14_Monocyte2_HC_Obese","CD14_Monocyte2_HC_Lean",
                                                  "CD14_Monocyte1_HC_Obese","CD14_Monocyte1_HC_Lean"
                                       ))
Idents(Mono.integrated_Poly)<-Mono.integrated_Poly$cell_disease_weight_treatments
Idents(Mono.integrated_Poly) <- factor(Idents(Mono.integrated_Poly),
                                       levels = c("CD14_Monocyte2_HD_Obese_Poly_Basal","CD14_Monocyte2_HD_Obese_Basal_Basal",
                                                  "CD14_Monocyte2_HD_Lean_Poly_Basal","CD14_Monocyte2_HD_Lean_Basal_Basal",
                                                  "CD14_Monocyte1_HD_Obese_Poly_Basal","CD14_Monocyte1_HD_Obese_Basal_Basal",
                                                  "CD14_Monocyte1_HD_Lean_Poly_Basal","CD14_Monocyte1_HD_Lean_Basal_Basal",
                                                  "CD14_Monocyte2_HC_Obese_Poly_Basal","CD14_Monocyte2_HC_Obese_Basal_Basal",
                                                  "CD14_Monocyte2_HC_Lean_Poly_Basal","CD14_Monocyte2_HC_Lean_Basal_Basal",
                                                  "CD14_Monocyte1_HC_Obese_Poly_Basal","CD14_Monocyte1_HC_Obese_Basal_Basal",
                                                  "CD14_Monocyte1_HC_Lean_Poly_Basal","CD14_Monocyte1_HC_Lean_Basal_Basal"
                                       ))

# This is Figure 4C
#Common upregulated genes for all conditions
markers.to.plot <- c(
"BIRC3","CCL20","CCL3","CCL3L3",
"CCL4","CCL4L2","CD83","CXCL1",
"CXCL2","CXCL3","CXCL8","EGR3",
"EHD1","F3","G0S2","GCH1",
"HEY1","HS3ST3B1","ICAM1","IER3",
"IL1A","IL1B","IL6","IRAK2",
"KANK1","KCNJ2","NFKBIA","NFKBIZ",
"NLRP3","PDE4B","PNRC1","PTGS2",
"REL","RNF144B","SERPINB9","SLC2A6",
"SOCS3","TNF","TNFAIP2","TNFAIP3",
"TNFAIP6","TNFAIP8","TNIP3","TRAF1",
"WTAP","ZC3H12A","ZC3H12C"
)
pdf("DotPlot_CD14Mono_Poly_TopCommon.pdf", height = 8, width = 14, useDingbats=FALSE)
DotPlot(Mono.integrated_Poly, features = (markers.to.plot), 
        dot.scale = 8, 
#        split.by = "first",
        assay = 'SCT'
) + RotatedAxis() + theme(axis.text.x = element_text(angle = 90))
dev.off()

# This is Supplemental Figure 4A
# Genes only upregulated in HD
markers.to.plot <- c(
"AC010980.1","ACOD1","ACSL1","ACVR2A",
"ADORA2A","AK4","ATP2B1","B4GALT5",
"BCL2A1","C11orf96","C22orf42","CCL15",
"CCRL2","CD274","CD40","CD83",
"CLEC4E","CLIC4","CSF1","CSF3",
"CXCL10","DDIT4","DRAM1","DUSP5",
"ELOVL7","ETV3","FLT1","GCH1",
"GEM","GJB2","GPR132","HCAR2",
"HCAR3","IFIT2","IL18","IL6",
"INHBA","JAG1","KCNJ2","KMO",
"LAMB3","LAMP3","LCP2","MAP3K8",
"MARCKS","MFSD2A","MSANTD3","MTF1",
"NBN","NFKB2","NR4A3","PALM2AKAP2",
"PELI1","PIM3","PLAU","PLD1",
"PLEKHF2","PNPLA1","PSMA6","PSTPIP2",
"PTGER2","PTX3","RGS16","RIN2",
"RIPK2","RND1","SAV1","SDC4",
"SERPINB2","SERPINB9","SGPP2","SLC2A6",
"SNX10","STAT5A","STX11","TFPI2",
"TJP2","TNFRSF4","TNIP1","TNIP3",
"TSLP","UPB1","USP12","ZBTB10"
)

pdf("DotPlot_CD14Mono_Poly_Top_HDonly.pdf", height = 8, width = 17, useDingbats=FALSE)
DotPlot(Mono.integrated_Poly, features = (markers.to.plot), 
        dot.scale = 8, 
        #split.by = "first",
        assay = 'SCT'
) + RotatedAxis() + theme(axis.text.x = element_text(angle = 90))
dev.off()


# LPS ONLY
Idents(pbmc.integrated) <- pbmc.integrated$treatments
pbmc.integrated_LPS <- subset(pbmc.integrated, idents = c("Basal_Basal","Basal_LPS"))
table(Idents(pbmc.integrated_LPS))

DefaultAssay(pbmc.integrated_LPS) <- "SCT"
Idents(pbmc.integrated_LPS) <- pbmc.integrated_LPS$cell_disease_weight
table(Idents(pbmc.integrated_LPS))
Idents(pbmc.integrated_LPS) <- factor(Idents(pbmc.integrated_LPS),
                                       levels = c("CD14_Monocyte1_HC_Lean","CD14_Monocyte1_HC_Obese",
                                                  "CD14_Monocyte1_HD_Lean","CD14_Monocyte1_HD_Obese",
                                                  "CD14_Monocyte2_HC_Lean","CD14_Monocyte2_HC_Obese",
                                                  "CD14_Monocyte2_HD_Lean","CD14_Monocyte2_HD_Obese",
                                                  "CD16_Monocyte_HC_Lean","CD16_Monocyte_HC_Obese",
                                                  "CD16_Monocyte_HD_Lean","CD16_Monocyte_HD_Obese"
                                       ))
plots <- VlnPlot(pbmc.integrated_LPS, features = c("TNF","IL1B"), split.by = "second",
                 pt.size = 0, combine = FALSE, 
                 idents = c("CD14_Monocyte1_HC_Lean","CD14_Monocyte1_HC_Obese",
                            "CD14_Monocyte1_HD_Lean","CD14_Monocyte1_HD_Obese"
                 ))
plots <- lapply(X = plots,
                FUN = function(p) p
                + ggplot2::theme(axis.title.y = element_blank())
                + ggplot2::scale_fill_manual(values = c('black', 'red'))
                )
CombinePlots(plots = plots, ncol = 1)



########################################
##                                    ##
##              Dot Plots             ##
##              for ISGs              ##
########################################

# List of all ISGs, collected from the literature
markers.to.plot <- c(
  "IFITM1","IFITM2","IFITM3","NCOA7","TRIM5",
  "MX1","MX2",
  "APOBEC1","APOBEC2","APOBEC3A","APOBEC3B","APOBEC3C",
  "APOBEC3D","APOBEC3F","APOBEC3G","APOBEC3H","APOBEC4","AICDA",
  "IFI16",
  "EIF2AK2","ZC3HAV1","PARP12",
  "IFI6","RSAD2",
  "IFIT1","IFIT2","IFIT3","IFIT5",
  "ISG20","OAS1","OAS2","OAS3",
  "BST2","CNP","GBP5",
  "IFI44L","ISG15","OASL"
)

markers.to.plot <- c(
  #Various

  # Egress
  "RSAD2","GBP5","CNP","BST2",
  
  # Replication
  "OAS3","OAS2","OAS1",
  "APOBEC3H","APOBEC3G","APOBEC3F",
  "APOBEC3D","APOBEC3C","APOBEC3A","ADAR",
  
  # Replication and Translation
  "DDX58","IFIH1",
  "IFIT5","IFIT3","IFIT2","IFIT1",
  "IFI6",
  
  # Translation
  "SSBP3","PARP12","OASL","NT5C3A",
  "MS4A4A","MAP3K14","IFI44L","EIF2AK2",
  "DDX60","DDIT4","CGAS",

  # Tranlsation an post Entry
  "ZC3HAV1",
  
  # RNA syn
  "TRIM5","MX1",
  "ISG20","IFI16",
  
  # post entry
  "MOV10",
  
  # Entry
  "NCOA7",
  "IFITM3","IFITM2","IFITM1"
  
  # IFN
  #"IRF1",
  #"IRF7", "TRIM25"
)

# Genes removed with low expression: CH25H,SFLN11,SAT1,APOBEC4,AICDA,TREX1,APOBEC1,APOBEC2,APOBEC3B
# Genes removed with unknown role: NAMPT MX2 HPSE P2RY6 JADE2 RTP4 SLC15A3 SLC25A28 SUN2 
# Genes removed associated with IFNG:  "GBP1","GBP2"
# Genes removed because they make no sense here  #"CD74",  #"ISG15","PML",

# PolyIC DotPlots
# Again, removing cells we don't need for the DotPlot
Idents(pbmc.integrated_Poly) <- pbmc.integrated_Poly$celltype
table(Idents(pbmc.integrated_Poly))
Mono.integrated_Poly <- subset(pbmc.integrated_Poly, idents = c("CD14_Monocyte1","CD14_Monocyte2","CD16_Monocyte", "DC"))
table(Idents(Mono.integrated_Poly))

# There are many ways to make these figures, and split the data in different ways.
# I will keep some options here,but I don't necessarily like them all

# Option 1 - Have DotPlot split the data using the metadata label treatments
# This section is commented out because I don't like it like this

#DefaultAssay(Mono.integrated_Poly) <- "SCT"
#Idents(Mono.integrated_Poly) <- Mono.integrated_Poly$cell_disease_weight #Format: "CD4-T-cell1_HC_Lean"
#table(Idents(Mono.integrated_Poly))

#Idents(Mono.integrated_Poly)<-factor(Idents(Mono.integrated_Poly), 
#                                levels = c("CD14_Monocyte1_HC_Lean","CD14_Monocyte1_HC_Obese",
#                                           "CD14_Monocyte1_HD_Lean","CD14_Monocyte1_HD_Obese",
#                                           "CD14_Monocyte2_HC_Lean","CD14_Monocyte2_HC_Obese",
#                                           "CD14_Monocyte2_HD_Lean","CD14_Monocyte2_HD_Obese",
#                                           "CD16_Monocyte_HC_Lean","CD16_Monocyte_HC_Obese",
#                                           "CD16_Monocyte_HD_Lean","CD16_Monocyte_HD_Obese",
#                                           "DC_HC_Lean","DC_HC_Obese",
#                                           "DC_HD_Lean","DC_HD_Obese"
#                                           ))
#table(Idents(Mono.integrated_Poly))

#pdf("DotPlot_CD14Mono1_HCHD_ISG_v3.pdf", height = 4, width = 16, useDingbats=FALSE)
#DotPlot(Mono.integrated_Poly, features = rev(markers.to.plot), 
#        idents = c("CD14_Monocyte1_HC_Lean","CD14_Monocyte1_HC_Obese",
#                   "CD14_Monocyte1_HD_Lean","CD14_Monocyte1_HD_Obese"
#                   ),
        #cols = c("black","blue"), 
#        dot.scale = 8, 
#        split.by = "treatments",
#        assay = 'SCT'
#) + RotatedAxis()
#dev.off()

# Option 2 - Manually order the cells and treatments they way you want it
DefaultAssay(Mono.integrated_Poly) <- "SCT"
Idents(Mono.integrated_Poly) <- Mono.integrated_Poly$cell_disease_weight_treatments #Format: "CD4-T-cell1_HC_Lean_Poly_Basal"
Idents(Mono.integrated_Poly)<-factor(Idents(Mono.integrated_Poly), 
                                     levels = c("CD14_Monocyte1_HC_Lean_Basal_Basal","CD14_Monocyte1_HC_Lean_Poly_Basal",
                                                "CD14_Monocyte1_HC_Obese_Basal_Basal","CD14_Monocyte1_HC_Obese_Poly_Basal",
                                                "CD14_Monocyte1_HD_Lean_Basal_Basal","CD14_Monocyte1_HD_Lean_Poly_Basal",
                                                "CD14_Monocyte1_HD_Obese_Basal_Basal","CD14_Monocyte1_HD_Obese_Poly_Basal",
                                                "CD14_Monocyte2_HC_Lean_Basal_Basal","CD14_Monocyte2_HC_Lean_Poly_Basal",
                                                "CD14_Monocyte2_HC_Obese_Basal_Basal","CD14_Monocyte2_HC_Obese_Poly_Basal",
                                                "CD14_Monocyte2_HD_Lean_Basal_Basal","CD14_Monocyte2_HD_Lean_Poly_Basal",
                                                "CD14_Monocyte2_HD_Obese_Basal_Basal","CD14_Monocyte2_HD_Obese_Poly_Basal",
                                                "CD16_Monocyte_HC_Lean_Basal_Basal","CD16_Monocyte_HC_Lean_Poly_Basal",
                                                "CD16_Monocyte_HC_Obese_Basal_Basal","CD16_Monocyte_HC_Obese_Poly_Basal",
                                                "CD16_Monocyte_HD_Lean_Basal_Basal","CD16_Monocyte_HD_Lean_Poly_Basal",
                                                "CD16_Monocyte_HD_Obese_Basal_Basal","CD16_Monocyte_HD_Obese_Poly_Basal",
                                                "DC_HC_Lean_Basal_Basal","DC_HC_Lean_Poly_Basal",
                                                "DC_HC_Obese_Basal_Basal","DC_HC_Obese_Poly_Basal",
                                                "DC_HD_Lean_Basal_Basal","DC_HD_Lean_Poly_Basal",
                                                "DC_HD_Obese_Basal_Basal","DC_HD_Obese_Poly_Basal"))
table(Idents(Mono.integrated_Poly))

DefaultAssay(Mono.integrated_Poly) <- "SCT"
Idents(Mono.integrated_Poly) <- Mono.integrated_Poly$cell_disease_weight_treatments #Format: "CD4-T-cell1_HC_Lean_Poly_Basal"
Idents(Mono.integrated_Poly)<-factor(Idents(Mono.integrated_Poly), 
                                     levels = c("DC_HD_Obese_Poly_Basal","DC_HD_Obese_Basal_Basal",
                                                "DC_HD_Lean_Poly_Basal","DC_HD_Lean_Basal_Basal",
                                                "DC_HC_Obese_Poly_Basal","DC_HC_Obese_Basal_Basal",
                                                "DC_HC_Lean_Poly_Basal","DC_HC_Lean_Basal_Basal",
                                              
                                                "CD14_Monocyte1_HD_Obese_Poly_Basal","CD14_Monocyte1_HD_Obese_Basal_Basal",
                                                "CD14_Monocyte1_HD_Lean_Poly_Basal","CD14_Monocyte1_HD_Lean_Basal_Basal",
                                                "CD14_Monocyte1_HC_Obese_Poly_Basal","CD14_Monocyte1_HC_Obese_Basal_Basal",
                                                "CD14_Monocyte1_HC_Lean_Poly_Basal","CD14_Monocyte1_HC_Lean_Basal_Basal",
                                                
                                                "CD14_Monocyte2_HD_Obese_Poly_Basal","CD14_Monocyte2_HD_Obese_Basal_Basal",
                                                "CD14_Monocyte2_HD_Lean_Poly_Basal","CD14_Monocyte2_HD_Lean_Basal_Basal",
                                                "CD14_Monocyte2_HC_Obese_Poly_Basal","CD14_Monocyte2_HC_Obese_Basal_Basal",
                                                "CD14_Monocyte2_HC_Lean_Poly_Basal","CD14_Monocyte2_HC_Lean_Basal_Basal",
                                                
                                                "CD16_Monocyte_HD_Obese_Poly_Basal","CD16_Monocyte_HD_Obese_Basal_Basal",
                                                "CD16_Monocyte_HD_Lean_Poly_Basal","CD16_Monocyte_HD_Lean_Basal_Basal",
                                                "CD16_Monocyte_HC_Obese_Poly_Basal","CD16_Monocyte_HC_Obese_Basal_Basal",
                                                "CD16_Monocyte_HC_Lean_Poly_Basal","CD16_Monocyte_HC_Lean_Basal_Basal"
                                                ))
table(Idents(Mono.integrated_Poly))

# This is Figure 2C
pdf("DotPlot_CD14Mono1_HCHD_ISG_v5.pdf", height = 4, width = 16, useDingbats=FALSE)
DotPlot(Mono.integrated_Poly, features = rev(markers.to.plot), 
        idents = c("CD14_Monocyte1_HC_Lean_Basal_Basal","CD14_Monocyte1_HC_Lean_Poly_Basal",
                   "CD14_Monocyte1_HC_Obese_Basal_Basal","CD14_Monocyte1_HC_Obese_Poly_Basal",
                   "CD14_Monocyte1_HD_Lean_Basal_Basal","CD14_Monocyte1_HD_Lean_Poly_Basal",
                   "CD14_Monocyte1_HD_Obese_Basal_Basal","CD14_Monocyte1_HD_Obese_Poly_Basal"),
        dot.scale = 8, 
        assay = 'SCT'
) + RotatedAxis() + theme(axis.text.x = element_text(angle = 90))
dev.off()

# This is Figure 2C
pdf("DotPlot_CD14Mono2_HCHD_ISG_v5.pdf", height = 4, width = 16, useDingbats=FALSE)
DotPlot(Mono.integrated_Poly, features = rev(markers.to.plot), 
        idents = c("CD14_Monocyte2_HC_Lean_Basal_Basal","CD14_Monocyte2_HC_Lean_Poly_Basal",
                   "CD14_Monocyte2_HC_Obese_Basal_Basal","CD14_Monocyte2_HC_Obese_Poly_Basal",
                   "CD14_Monocyte2_HD_Lean_Basal_Basal","CD14_Monocyte2_HD_Lean_Poly_Basal",
                   "CD14_Monocyte2_HD_Obese_Basal_Basal","CD14_Monocyte2_HD_Obese_Poly_Basal"),
        dot.scale = 8, 
        assay = 'SCT'
) + RotatedAxis() + theme(axis.text.x = element_text(angle = 90))
dev.off()

# This is Figure 2C
pdf("DotPlot_CD16Mono_HCHD_ISG_v5.pdf", height = 4, width = 16, useDingbats=FALSE)
DotPlot(Mono.integrated_Poly, features = rev(markers.to.plot), 
        idents = c("CD16_Monocyte_HC_Lean_Basal_Basal","CD16_Monocyte_HC_Lean_Poly_Basal",
                   "CD16_Monocyte_HC_Obese_Basal_Basal","CD16_Monocyte_HC_Obese_Poly_Basal",
                   "CD16_Monocyte_HD_Lean_Basal_Basal","CD16_Monocyte_HD_Lean_Poly_Basal",
                   "CD16_Monocyte_HD_Obese_Basal_Basal","CD16_Monocyte_HD_Obese_Poly_Basal"),
        dot.scale = 8, 
        assay = 'SCT'
) + RotatedAxis() + theme(axis.text.x = element_text(angle = 90))
dev.off()

# This is an Unused Figure
pdf("DotPlot_DC_HCHD_ISG_v5.pdf", height = 4, width = 16, useDingbats=FALSE)
DotPlot(Mono.integrated_Poly, features = rev(markers.to.plot), 
        idents = c("DC_HC_Lean_Basal_Basal","DC_HC_Lean_Poly_Basal",
                   "DC_HC_Obese_Basal_Basal","DC_HC_Obese_Poly_Basal",
                   "DC_HD_Lean_Basal_Basal","DC_HD_Lean_Poly_Basal",
                   "DC_HD_Obese_Basal_Basal","DC_HD_Obese_Poly_Basal"),
        dot.scale = 8, 
        assay = 'SCT'
) + RotatedAxis()   + theme(axis.text.x = element_text(angle = 90)) 
dev.off()



# LPS DotPLots
# Again, removing cells we don't need for the DotPlot
Idents(pbmc.integrated_LPS) <- pbmc.integrated_LPS$celltype
table(Idents(pbmc.integrated_LPS))
Mono.integrated_LPS <- subset(pbmc.integrated_LPS, idents = c("CD14_Monocyte1","CD14_Monocyte2","CD16_Monocyte","DC"))
table(Idents(Mono.integrated_LPS))

# Option 1 - Have DotPlot split the data using the metadata label treatments
# This section is commented out because I don't like it like this

#DefaultAssay(Mono.integrated_LPS) <- "SCT"
#Idents(Mono.integrated_LPS) <- Mono.integrated_LPS$cell_disease_weight #Format: "CD4-T-cell1_HC_Lean"
#table(Idents(Mono.integrated_LPS))

#Idents(Mono.integrated_LPS)<-factor(Idents(Mono.integrated_LPS), 
#                                     levels = c("CD14_Monocyte1_HC_Lean","CD14_Monocyte1_HC_Obese",
#                                                "CD14_Monocyte1_HD_Lean","CD14_Monocyte1_HD_Obese",
#                                                "CD14_Monocyte2_HC_Lean","CD14_Monocyte2_HC_Obese",
#                                                "CD14_Monocyte2_HD_Lean","CD14_Monocyte2_HD_Obese",
#                                                "CD16_Monocyte_HC_Lean","CD16_Monocyte_HC_Obese",
#                                                "CD16_Monocyte_HD_Lean","CD16_Monocyte_HD_Obese",
#                                                "DC_HC_Lean","DC_HC_Obese",
#                                                "DC_HD_Lean","DC_HD_Obese"
#                                     ))
#table(Idents(Mono.integrated_LPS))

#pdf("DotPlot_CD14Mono1_HCHD_ISG_LPS_v3.pdf", height = 4, width = 16, useDingbats=FALSE)
#DotPlot(Mono.integrated_LPS, features = rev(markers.to.plot), 
#        idents = c("CD14_Monocyte1_HC_Lean","CD14_Monocyte1_HC_Obese",
#                   "CD14_Monocyte1_HD_Lean","CD14_Monocyte1_HD_Obese"
#        ),
        #cols = c("black","blue"), 
#        dot.scale = 8, 
#        split.by = "treatments",
#        assay = 'SCT'
#) + RotatedAxis()
#dev.off()

# Option 2 - Manually order the cells and treatments they way you want it
DefaultAssay(Mono.integrated_LPS) <- "SCT"
Idents(Mono.integrated_LPS) <- Mono.integrated_LPS$cell_disease_weight_treatments #Format: "CD4-T-cell1_HC_Lean_Poly_Basal"
Idents(Mono.integrated_LPS)<-factor(Idents(Mono.integrated_LPS), 
                                     levels = c("CD14_Monocyte1_HC_Lean_Basal_Basal","CD14_Monocyte1_HC_Lean_Basal_LPS",
                                                "CD14_Monocyte1_HC_Obese_Basal_Basal","CD14_Monocyte1_HC_Obese_Basal_LPS",
                                                "CD14_Monocyte1_HD_Lean_Basal_Basal","CD14_Monocyte1_HD_Lean_Basal_LPS",
                                                "CD14_Monocyte1_HD_Obese_Basal_Basal","CD14_Monocyte1_HD_Obese_Basal_LPS",
                                                "CD14_Monocyte2_HC_Lean_Basal_Basal","CD14_Monocyte2_HC_Lean_Basal_LPS",
                                                "CD14_Monocyte2_HC_Obese_Basal_Basal","CD14_Monocyte2_HC_Obese_Basal_LPS",
                                                "CD14_Monocyte2_HD_Lean_Basal_Basal","CD14_Monocyte2_HD_Lean_Basal_LPS",
                                                "CD14_Monocyte2_HD_Obese_Basal_Basal","CD14_Monocyte2_HD_Obese_Basal_LPS",
                                                "CD16_Monocyte_HC_Lean_Basal_Basal","CD16_Monocyte_HC_Lean_Basal_LPS",
                                                "CD16_Monocyte_HC_Obese_Basal_Basal","CD16_Monocyte_HC_Obese_Basal_LPS",
                                                "CD16_Monocyte_HD_Lean_Basal_Basal","CD16_Monocyte_HD_Lean_Basal_LPS",
                                                "CD16_Monocyte_HD_Obese_Basal_Basal","CD16_Monocyte_HD_Obese_Basal_LPS",
                                                "DC_HC_Lean_Basal_Basal","DC_HC_Lean_Basal_LPS",
                                                "DC_HC_Obese_Basal_Basal","DC_HC_Obese_Basal_LPS",
                                                "DC_HD_Lean_Basal_Basal","DC_HD_Lean_Basal_LPS",
                                                "DC_HD_Obese_Basal_Basal","DC_HD_Obese_Basal_LPS"))

Idents(Mono.integrated_LPS)<-factor(Idents(Mono.integrated_LPS), 
                                    levels = c(
                                               "DC_HD_Obese_Basal_LPS","DC_HD_Obese_Basal_Basal",
                                               "DC_HD_Lean_Basal_LPS","DC_HD_Lean_Basal_Basal",
                                               "DC_HC_Obese_Basal_LPS","DC_HC_Obese_Basal_Basal",
                                               "DC_HC_Lean_Basal_LPS","DC_HC_Lean_Basal_Basal",
                                               
                                               "CD14_Monocyte1_HD_Obese_Basal_LPS","CD14_Monocyte1_HD_Obese_Basal_Basal",
                                               "CD14_Monocyte1_HD_Lean_Basal_LPS","CD14_Monocyte1_HD_Lean_Basal_Basal",
                                               "CD14_Monocyte1_HC_Obese_Basal_LPS","CD14_Monocyte1_HC_Obese_Basal_Basal",
                                               "CD14_Monocyte1_HC_Lean_Basal_LPS","CD14_Monocyte1_HC_Lean_Basal_Basal",
                                               
                                               "CD14_Monocyte2_HD_Obese_Basal_LPS","CD14_Monocyte2_HD_Obese_Basal_Basal",
                                               "CD14_Monocyte2_HD_Lean_Basal_LPS","CD14_Monocyte2_HD_Lean_Basal_Basal",
                                               "CD14_Monocyte2_HC_Obese_Basal_LPS","CD14_Monocyte2_HC_Obese_Basal_Basal",
                                               "CD14_Monocyte2_HC_Lean_Basal_LPS","CD14_Monocyte2_HC_Lean_Basal_Basal",
                                               
                                               "CD16_Monocyte_HD_Obese_Basal_LPS","CD16_Monocyte_HD_Obese_Basal_Basal",
                                               "CD16_Monocyte_HD_Lean_Basal_LPS","CD16_Monocyte_HD_Lean_Basal_Basal",
                                               "CD16_Monocyte_HC_Obese_Basal_LPS","CD16_Monocyte_HC_Obese_Basal_Basal",
                                               "CD16_Monocyte_HC_Lean_Basal_LPS","CD16_Monocyte_HC_Lean_Basal_Basal"
                                               ))

table(Idents(Mono.integrated_LPS))

# This is Figure 2D
pdf("DotPlot_CD14Mono1_HCHD_ISG_LPS_v5.pdf", height = 4, width = 16, useDingbats=FALSE)
DotPlot(Mono.integrated_LPS, features = rev(markers.to.plot), 
        idents = c("CD14_Monocyte1_HC_Lean_Basal_Basal","CD14_Monocyte1_HC_Lean_Basal_LPS",
                   "CD14_Monocyte1_HC_Obese_Basal_Basal","CD14_Monocyte1_HC_Obese_Basal_LPS",
                   "CD14_Monocyte1_HD_Lean_Basal_Basal","CD14_Monocyte1_HD_Lean_Basal_LPS",
                   "CD14_Monocyte1_HD_Obese_Basal_Basal","CD14_Monocyte1_HD_Obese_Basal_LPS"),
        dot.scale = 8, 
        assay = 'SCT'
) + RotatedAxis() + theme(axis.text.x = element_text(angle = 90))
dev.off()

# This is Figure 2D
pdf("DotPlot_CD14Mono2_HCHD_ISG_LPS_v5.pdf", height = 4, width = 16, useDingbats=FALSE)
DotPlot(Mono.integrated_LPS, features = rev(markers.to.plot), 
        idents = c("CD14_Monocyte2_HC_Lean_Basal_Basal","CD14_Monocyte2_HC_Lean_Basal_LPS",
                   "CD14_Monocyte2_HC_Obese_Basal_Basal","CD14_Monocyte2_HC_Obese_Basal_LPS",
                   "CD14_Monocyte2_HD_Lean_Basal_Basal","CD14_Monocyte2_HD_Lean_Basal_LPS",
                   "CD14_Monocyte2_HD_Obese_Basal_Basal","CD14_Monocyte2_HD_Obese_Basal_LPS"),
        dot.scale = 8, 
        assay = 'SCT'
) + RotatedAxis() + theme(axis.text.x = element_text(angle = 90))
dev.off()

# This is Figure 2D
pdf("DotPlot_CD16Mono_HCHD_ISG_LPS_v5.pdf", height = 4, width = 16, useDingbats=FALSE)
DotPlot(Mono.integrated_LPS, features = rev(markers.to.plot), 
        idents = c("CD16_Monocyte_HC_Lean_Basal_Basal","CD16_Monocyte_HC_Lean_Basal_LPS",
                   "CD16_Monocyte_HC_Obese_Basal_Basal","CD16_Monocyte_HC_Obese_Basal_LPS",
                   "CD16_Monocyte_HD_Lean_Basal_Basal","CD16_Monocyte_HD_Lean_Basal_LPS",
                   "CD16_Monocyte_HD_Obese_Basal_Basal","CD16_Monocyte_HD_Obese_Basal_LPS"),
        dot.scale = 8, 
        assay = 'SCT'
) + RotatedAxis() + theme(axis.text.x = element_text(angle = 90))
dev.off()

# This is an Unused Figure
pdf("DotPlot_DC_HCHD_ISG_LPS_v5.pdf", height = 4, width = 16, useDingbats=FALSE)
DotPlot(Mono.integrated_LPS, features = rev(markers.to.plot), 
        idents = c("DC_HC_Lean_Basal_Basal","DC_HC_Lean_Basal_LPS",
                   "DC_HC_Obese_Basal_Basal","DC_HC_Obese_Basal_LPS",
                   "DC_HD_Lean_Basal_Basal","DC_HD_Lean_Basal_LPS",
                   "DC_HD_Obese_Basal_Basal","DC_HD_Obese_Basal_LPS"),
        dot.scale = 8, 
        assay = 'SCT'
) + RotatedAxis() + theme(axis.text.x = element_text(angle = 90))
dev.off()


########################################
##                                    ##
##       Differential Expression      ##
##                                    ##
########################################

library(Libra)

DefaultAssay(pbmc.integrated) <- "SCT"

#data("hagai_toy")
#test<-hagai_toy@meta.data
#DE = run_de(hagai_toy)
#head(DE)
#matrices = to_pseudobulk(hagai_toy, meta = meta)

# Adam's notes on Libra
# Can be run on a seurat object
# Requires three metadata categories:
# 1. label = This is the comparison label - I believe it can only take 2 variables
# 2. celltype = This is the group of cells - any number of variables, even 1, and yes, every comparison is made 
# 3. replicate = source of data (I think we'll use patient name)

# I tested doing all three diseases at the same tive vs doing just HC alone, and stats are the same


pbmc.integrated$patient_cell_disease_weight_treatments <- paste(pbmc.integrated$celltype,
                                                                pbmc.integrated$disease,
                                                                pbmc.integrated$weight,
                                                                pbmc.integrated$orig.ident,
                                                                pbmc.integrated$first, 
                                                                pbmc.integrated$second, sep = "_")

# These names are messed up, technically I wouldn't want the word disease in the variable
pbmc.integrated$cell_disease_first <- paste(pbmc.integrated$celltype,
                                            pbmc.integrated$first, sep = "_")
pbmc.integrated$cell_disease_second <- paste(pbmc.integrated$celltype,
                                             pbmc.integrated$second, sep = "_")



getwd()
setwd("I:/Adam/Unix/scRNA-seq/Marti_PBMC/Ranalysis/DE")

#######################
# Basal vs Treatments #
#######################
for (spec_cell in c("CD4-T-cell1", "CD4-T-cell2", 
                    "CD4-T-cell3", "CD4-T-cell4", 
                    "CD4-T-cell5", "CD4-T-cell6", 
                    "CD8-T-cell1", "CD8-T-cell2", 
                    "B-cell",
                    "CD14_Monocyte1", "CD14_Monocyte2", "CD16_Monocyte",
                    "NK-cell",
                    "DC")){
  Cell_Basal <- (paste(spec_cell, "_Basal", sep = ""))
  Cell_Poly <- (paste(spec_cell, "_Poly", sep = ""))
  Cell_LPS <- (paste(spec_cell, "_LPS", sep = ""))
  
  Idents(pbmc.integrated)<-pbmc.integrated$cell_disease_first # Categorize by POLYIC status
  All_Mono <- subset(pbmc.integrated, idents = c(Cell_Basal)) # No PolyIC, +/-LPS
  All_Mono$replicate<-All_Mono$orig.ident # Patient number are replicates
  All_Mono$label<-All_Mono$second # Define the comparison between LPS and Basal
  All_Mono$cell_type<-All_Mono$disease_weight # Test all diseases
  DE_LPS = run_de(All_Mono)
  
  Comp_File <- (paste(spec_cell, "_Basal_LPS_DE.txt", sep = ""))
  write.table(DE_LPS, file=Comp_File)

  HClean_Pseudo_File <- (paste(spec_cell, "_HC_lean_Basal_LPS_pseudobulk.txt", sep = ""))
  HDlean_Pseudo_File <- (paste(spec_cell, "_HD_lean_Basal_LPS_pseudobulk.txt", sep = ""))
  AHlean_Pseudo_File <- (paste(spec_cell, "_AH_lean_Basal_LPS_pseudobulk.txt", sep = ""))
  HCobese_Pseudo_File <- (paste(spec_cell, "_HC_obese_Basal_LPS_pseudobulk.txt", sep = ""))
  HDobese_Pseudo_File <- (paste(spec_cell, "_HD_obese_Basal_LPS_pseudobulk.txt", sep = ""))
  AHobese_Pseudo_File <- (paste(spec_cell, "_AH_obese_Basal_LPS_pseudobulk.txt", sep = ""))
  matrices = to_pseudobulk(All_Mono, meta = meta)
  write.table(matrices$HC_lean, file=HClean_Pseudo_File)
  write.table(matrices$HD_lean, file=HDlean_Pseudo_File)
  write.table(matrices$AH_lean, file=AHlean_Pseudo_File)
  write.table(matrices$HC_obese, file=HCobese_Pseudo_File)
  write.table(matrices$HD_obese, file=HDobese_Pseudo_File)
  write.table(matrices$AH_obese, file=AHobese_Pseudo_File)
############################### 
  
  
  
###############################
  Idents(pbmc.integrated)<-pbmc.integrated$cell_disease_second # Categorize by LPS status
  All_Mono <- subset(pbmc.integrated, idents = c(Cell_Basal)) # No LPS, +/-PolyIC
  All_Mono$replicate<-All_Mono$orig.ident # Patient number are replicates
  All_Mono$label<-All_Mono$first # Define the comparison between PolyIC and Basal
  All_Mono$cell_type<-All_Mono$disease_weight # Test all diseases
  DE_Poly = run_de(All_Mono)

  Comp_File <- (paste(spec_cell, "_Poly_Basal_DE.txt", sep = ""))
  write.table(DE_Poly, file=Comp_File)

  HClean_Pseudo_File <- (paste(spec_cell, "_HC_lean_Poly_Basal_pseudobulk.txt", sep = ""))
  HDlean_Pseudo_File <- (paste(spec_cell, "_HD_lean_Poly_Basal_pseudobulk.txt", sep = ""))
  AHlean_Pseudo_File <- (paste(spec_cell, "_AH_lean_Poly_Basal_pseudobulk.txt", sep = ""))
  HCobese_Pseudo_File <- (paste(spec_cell, "_HC_obese_Poly_Basal_pseudobulk.txt", sep = ""))
  HDobese_Pseudo_File <- (paste(spec_cell, "_HD_obese_Poly_Basal_pseudobulk.txt", sep = ""))
  AHobese_Pseudo_File <- (paste(spec_cell, "_AH_obese_Poly_Basal_pseudobulk.txt", sep = ""))
  matrices = to_pseudobulk(All_Mono, meta = meta)
  write.table(matrices$HC_lean, file=HClean_Pseudo_File)
  write.table(matrices$HD_lean, file=HDlean_Pseudo_File)
  write.table(matrices$AH_lean, file=AHlean_Pseudo_File)
  write.table(matrices$HC_obese, file=HCobese_Pseudo_File)
  write.table(matrices$HD_obese, file=HDobese_Pseudo_File)
  write.table(matrices$AH_obese, file=AHobese_Pseudo_File)
###############################  
  
  
###############################
  Idents(pbmc.integrated)<-pbmc.integrated$cell_disease_first # Categorize by PolyIC status
  All_Mono <- subset(pbmc.integrated, idents = c(Cell_Poly)) # Only PolyIC, +/-LPS
  All_Mono$replicate<-All_Mono$orig.ident # Patient number are replicates
  All_Mono$label<-All_Mono$second # Define the comparison between PolyLPS and POlyICBasal
  All_Mono$cell_type<-All_Mono$disease_weight # Test all diseases
  DE_LPS = run_de(All_Mono)
  
  Comp_File <- (paste(spec_cell, "_PolyLPS_Poly_DE.txt", sep = ""))
  write.table(DE_LPS, file=Comp_File)
  
  HClean_Pseudo_File <- (paste(spec_cell, "_HC_lean_PolyLPS_Poly_pseudobulk.txt", sep = ""))
  HDlean_Pseudo_File <- (paste(spec_cell, "_HD_lean_PolyLPS_Poly_pseudobulk.txt", sep = ""))
  AHlean_Pseudo_File <- (paste(spec_cell, "_AH_lean_PolyLPS_Poly_pseudobulk.txt", sep = ""))
  HCobese_Pseudo_File <- (paste(spec_cell, "_HC_obese_PolyLPS_Poly_pseudobulk.txt", sep = ""))
  HDobese_Pseudo_File <- (paste(spec_cell, "_HD_obese_PolyLPS_Poly_pseudobulk.txt", sep = ""))
  AHobese_Pseudo_File <- (paste(spec_cell, "_AH_obese_PolyLPS_Poly_pseudobulk.txt", sep = ""))
  matrices = to_pseudobulk(All_Mono, meta = meta)
  write.table(matrices$HC_lean, file=HClean_Pseudo_File)
  write.table(matrices$HD_lean, file=HDlean_Pseudo_File)
  write.table(matrices$AH_lean, file=AHlean_Pseudo_File)
  write.table(matrices$HC_obese, file=HCobese_Pseudo_File)
  write.table(matrices$HD_obese, file=HDobese_Pseudo_File)
  write.table(matrices$AH_obese, file=AHobese_Pseudo_File)
############################### 
  
###############################
  Idents(pbmc.integrated)<-pbmc.integrated$cell_disease_second # Categorize by LPS status
  All_Mono <- subset(pbmc.integrated, idents = c(Cell_LPS)) # Only LPS, +/-PolyIC
  All_Mono$replicate<-All_Mono$orig.ident # Patient number are replicates
  All_Mono$label<-All_Mono$first # Define the comparison between PolyLPS and BasalLPS
  All_Mono$cell_type<-All_Mono$disease_weight # Test all diseases
  DE_LPS = run_de(All_Mono)
  
  Comp_File <- (paste(spec_cell, "_PolyLPS_LPS_DE.txt", sep = ""))
  write.table(DE_LPS, file=Comp_File)

  AH_Pseudo_File <- (paste(spec_cell, "_AH_PolyLPS_LPS_pseudobulk.txt", sep = ""))
  HClean_Pseudo_File <- (paste(spec_cell, "_HC_lean_PolyLPS_LPS_pseudobulk.txt", sep = ""))
  HDlean_Pseudo_File <- (paste(spec_cell, "_HD_lean_PolyLPS_LPS_pseudobulk.txt", sep = ""))
  AHlean_Pseudo_File <- (paste(spec_cell, "_AH_lean_PolyLPS_LPS_pseudobulk.txt", sep = ""))
  HCobese_Pseudo_File <- (paste(spec_cell, "_HC_obese_PolyLPS_LPS_pseudobulk.txt", sep = ""))
  HDobese_Pseudo_File <- (paste(spec_cell, "_HD_obese_PolyLPS_LPS_pseudobulk.txt", sep = ""))
  AHobese_Pseudo_File <- (paste(spec_cell, "_AH_obese_PolyLPS_LPS_pseudobulk.txt", sep = ""))
  matrices = to_pseudobulk(All_Mono, meta = meta)
  write.table(matrices$HC_lean, file=HClean_Pseudo_File)
  write.table(matrices$HD_lean, file=HDlean_Pseudo_File)
  write.table(matrices$AH_lean, file=AHlean_Pseudo_File)
  write.table(matrices$HC_obese, file=HCobese_Pseudo_File)
  write.table(matrices$HD_obese, file=HDobese_Pseudo_File)
  write.table(matrices$AH_obese, file=AHobese_Pseudo_File)
###############################   
}




#######################
#      HC vs HD       #
#        Lean         #
#######################
Idents(pbmc.integrated)<-pbmc.integrated$weight # Categorize by Disease
All_Mono <- subset(pbmc.integrated, idents = c("Lean")) # get rid of Obese

for (spec_cell in c("CD4-T-cell1", "CD4-T-cell2", 
                    "CD4-T-cell3", "CD4-T-cell4", 
                    "CD4-T-cell5", "CD4-T-cell6", 
                    "CD8-T-cell1", "CD8-T-cell2", 
                    "B-cell",
                    "CD14_Monocyte1", "CD14_Monocyte2", "CD16_Monocyte",
                    "NK-cell",
                    "DC")){

  ###############################

  Idents(All_Mono)<-All_Mono$celltype # Categorize by celltype
  All_Mono_sub <- subset(All_Mono, idents = c(spec_cell)) # No PolyIC, +/-LPS
  
  All_Mono_sub$replicate<-All_Mono_sub$orig.ident # Patient number are replicates
  All_Mono_sub$label<-All_Mono_sub$disease # Define the comparison between HC and HD
  All_Mono_sub$cell_type<-All_Mono_sub$treatments # Test all treatment groups (Basal_Basal format)
  DE_LPS = run_de(All_Mono_sub)
  
  Comp_File <- (paste(spec_cell, "_HCvHD_Lean_DE.txt", sep = ""))
  write.table(DE_LPS, file=Comp_File)
  
  Basal_Pseudo_File <- (paste(spec_cell, "_HCvHD_Lean_Basal_pseudobulk.txt", sep = ""))
  LPS_Pseudo_File <- (paste(spec_cell, "_HCvHD_Lean_LPS_pseudobulk.txt", sep = ""))
  Poly_Pseudo_File <- (paste(spec_cell, "_HCvHD_Lean_Poly_pseudobulk.txt", sep = ""))
  PolyLPS_Pseudo_File <- (paste(spec_cell, "_HCvHD_Lean_PolyLPS_pseudobulk.txt", sep = ""))
  matrices = to_pseudobulk(All_Mono_sub, meta = meta)
  write.table(matrices$Basal_Basal, file=Basal_Pseudo_File)
  write.table(matrices$Basal_LPS, file=LPS_Pseudo_File)
  write.table(matrices$Poly_Basal, file=Poly_Pseudo_File)
  write.table(matrices$Poly_LPS, file=PolyLPS_Pseudo_File)
  ############################### 
}

#######################
#      HC vs HD       #
#        Obese        #
#######################

Idents(pbmc.integrated)<-pbmc.integrated$weight # Categorize by Disease
All_Mono <- subset(pbmc.integrated, idents = c("Obese")) # get rid of Lean

for (spec_cell in c("CD4-T-cell1", "CD4-T-cell2", 
                    "CD4-T-cell3", "CD4-T-cell4", 
                    "CD4-T-cell5", "CD4-T-cell6", 
                    "CD8-T-cell1", "CD8-T-cell2", 
                    "B-cell",
                    "CD14_Monocyte1", "CD14_Monocyte2", "CD16_Monocyte",
                    "NK-cell",
                    "DC")){

  ###############################
  
  Idents(All_Mono)<-All_Mono$celltype # Categorize by celltype
  All_Mono_sub <- subset(All_Mono, idents = c(spec_cell)) 
  
  All_Mono_sub$replicate<-All_Mono_sub$orig.ident # Patient number are replicates
  All_Mono_sub$label<-All_Mono_sub$disease # Define the comparison between HC and HD
  All_Mono_sub$cell_type<-All_Mono_sub$treatments # Test all treatment groups (Basal_Basal format)
  DE_LPS = run_de(All_Mono_sub)
  
  Comp_File <- (paste(spec_cell, "_HCvHD_Obese_DE.txt", sep = ""))
  write.table(DE_LPS, file=Comp_File)
  
  Basal_Pseudo_File <- (paste(spec_cell, "_HCvHD_Obese_Basal_pseudobulk.txt", sep = ""))
  LPS_Pseudo_File <- (paste(spec_cell, "_HCvHD_Obese_LPS_pseudobulk.txt", sep = ""))
  Poly_Pseudo_File <- (paste(spec_cell, "_HCvHD_Obese_Poly_pseudobulk.txt", sep = ""))
  PolyLPS_Pseudo_File <- (paste(spec_cell, "_HCvHD_Obese_PolyLPS_pseudobulk.txt", sep = ""))
  matrices = to_pseudobulk(All_Mono_sub, meta = meta)
  write.table(matrices$Basal_Basal, file=Basal_Pseudo_File)
  write.table(matrices$Basal_LPS, file=LPS_Pseudo_File)
  write.table(matrices$Poly_Basal, file=Poly_Pseudo_File)
  write.table(matrices$Poly_LPS, file=PolyLPS_Pseudo_File)
  ############################### 
}



#######################
#      HC vs HD       #
#        Both         #
#######################

Idents(pbmc.integrated)<-pbmc.integrated$weight # Categorize by Disease
All_Mono <- subset(pbmc.integrated, idents = c("Lean","Obese")) # get rid of Lean

for (spec_cell in c("CD4-T-cell1", "CD4-T-cell2", 
                    "CD4-T-cell3", "CD4-T-cell4", 
                    "CD4-T-cell5", "CD4-T-cell6", 
                    "CD8-T-cell1", "CD8-T-cell2", 
                    "B-cell",
                    "CD14_Monocyte1", "CD14_Monocyte2", "CD16_Monocyte",
                    "NK-cell",
                    "DC")){
  
  ###############################
  
  Idents(All_Mono)<-All_Mono$celltype # Categorize by celltype
  All_Mono_sub <- subset(All_Mono, idents = c(spec_cell)) 
  
  All_Mono_sub$replicate<-All_Mono_sub$orig.ident # Patient number are replicates
  All_Mono_sub$label<-All_Mono_sub$disease # Define the comparison between HC and HD
  All_Mono_sub$cell_type<-All_Mono_sub$treatments # Test all treatment groups (Basal_Basal format)
  DE_LPS = run_de(All_Mono_sub)
  
  Comp_File <- (paste(spec_cell, "_HCvHD_Both_DE.txt", sep = ""))
  write.table(DE_LPS, file=Comp_File)
  
  Basal_Pseudo_File <- (paste(spec_cell, "_HCvHD_Both_Basal_pseudobulk.txt", sep = ""))
  LPS_Pseudo_File <- (paste(spec_cell, "_HCvHD_Both_LPS_pseudobulk.txt", sep = ""))
  Poly_Pseudo_File <- (paste(spec_cell, "_HCvHD_Both_Poly_pseudobulk.txt", sep = ""))
  PolyLPS_Pseudo_File <- (paste(spec_cell, "_HCvHD_Both_PolyLPS_pseudobulk.txt", sep = ""))
  matrices = to_pseudobulk(All_Mono_sub, meta = meta)
  write.table(matrices$Basal_Basal, file=Basal_Pseudo_File)
  write.table(matrices$Basal_LPS, file=LPS_Pseudo_File)
  write.table(matrices$Poly_Basal, file=Poly_Pseudo_File)
  write.table(matrices$Poly_LPS, file=PolyLPS_Pseudo_File)
  ############################### 
}





########################################
##                                    ##
##           Myeloid Subset           ##
##                                    ##
########################################

# This R object contains only the Monocytes in a Seurat object
# Not only were these cells subset from the original, 
# but they were re-put together as a new object
# and thus reclustered.
# Contact me if you'd like to discuss ways to share this object
# The code needed to make this is on a separate R script for simplicity.
# I typically would run that on the HPC, because clsutering can take a while.
# But you probably can do it on your computer, depends.

load("XXX/PolyICFigures/Mono_20220628.Robj")

# Example
#load("XXX/PolyICFigures/Mono_20220628.Robj")

# This is Figure 3A
pdf("UMAP_Mono_all.pdf", height = 6, width = 7, useDingbats=FALSE)
DimPlot(Mono.integrated, label = TRUE, label.size = 5) #700x600
dev.off()

# For your records
VlnPlot(Mono.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Preparatory stuff
DefaultAssay(Mono.integrated) <- "SCT"

Mono.integrated$mono_clus <- Idents(Mono.integrated)
Mono.integrated$both <- paste(Mono.integrated$first, Mono.integrated$second, sep = "_")
Mono.integrated$uber <- paste(Mono.integrated$disease, 
                              Mono.integrated$first, Mono.integrated$second, sep = "_")
Mono.integrated$disease_weight <- paste(Mono.integrated$disease, 
                                        Mono.integrated$weight, sep = "_")

# Cell numbers
Mono.integrated$Mono_patient_cell_disease_weight_treatments <- paste(Mono.integrated$mono_clus,
                                                                     Mono.integrated$disease,
                                                                     Mono.integrated$weight,
                                                                     Mono.integrated$orig.ident,
                                                                     Mono.integrated$first, 
                                                                     Mono.integrated$second, sep = "_")

cell.num <- as.data.frame(table(Mono.integrated$Mono_patient_cell_disease_weight_treatments))
write.table(cell.num, file="CellTypes_HCandHD_Patient_Mono_scRNA.txt")

Mono.integrated$Mono_patient_cell_disease_weight_treatments <- paste(Mono.integrated$mono_clus,
                                                                     Mono.integrated$disease,
                                                                     Mono.integrated$weight,
                                                                     Mono.integrated$orig.ident,
                                                                     Mono.integrated$first, 
                                                                     Mono.integrated$second, sep = "_")

cell.num <- as.data.frame(table(Mono.integrated$Mono_patient_cell_disease_weight_treatments))
write.table(cell.num, file="CellTypes_HCandHD_Patient_Mono_scRNA.txt")


########################################
##                                    ##
##                UMAPS               ##
##                                    ##
########################################

# This is Supplemental Figure 2A

Idents(Mono.integrated) <- Mono.integrated$mono_clus

pdf("UMAP_allCond.pdf", height = 10, width = 12, useDingbats=FALSE)
DimPlot(Mono.integrated, reduction = 'umap', split.by = 'ultra', ncol = 4)
dev.off()

Idents(Mono.integrated) <- Mono.integrated$seurat_clusters

# This is Supplemental Figure 2B
DimPlot(Mono.integrated, reduction = 'umap', split.by = 'seurat_clusters', ncol = 4)

Idents(Mono.integrated) <- Mono.integrated$seurat_clusters



########################################
#     Cluster Markers for Monocytes    #
########################################
# Subset only the Monocytes and Basal 
# Because the other treatments would influence the baseline markers

# Only Basal
DefaultAssay(Mono.integrated) <- "integrated"
Idents(Mono.integrated) <- Mono.integrated$both 
table(Idents(Mono.integrated))
Mono.integrated_subset <- subset(Mono.integrated, idents = c("Basal_Basal"))
Idents(Mono.integrated_subset)  <- Mono.integrated_subset$mono_clus
table(Idents(Mono.integrated_subset))

# find markers for every cluster compared to all remaining cells, report only the positive ones
Mono.markers <- FindAllMarkers(Mono.integrated_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(Mono.markers, file="ClusterMarkers_Mono_Basal_scRNA.txt")

Idents(Mono.integrated) <- Mono.integrated$mono_clus
VlnPlot(Mono.integrated, features = c("GNLY", "CD14", "CLEC4E", "IL7R", "CD3D", "CLEC7A"), 
        pt.size = 0, ncol = 3)

# This is Figure 3B
# Dot plot
markers.to.plot <- c(
  "SYNE2","CD2","BCL11B",
  "TMSB4X","VMO1","SERPINA1",
  "RPL13","RPS18","RPS12",
  "CD163","TGFBI","SLC4A7",
  "CTSB","CTSD","CTSL",
  "CCL24","FBP1","CYP1B1",
  "CXCL8","SOD2","IL1B",
  "DSE","NAMPT","VCAN",
  "HLA-DRA","CD74","HLA-DRB1",
  "LAP3","PARP14","IFI44L",
  "CSTB","TXN","PTPRE",
  "HLA-DPA1","HLA-DQA1","LYZ",
  "DDX21","FN1","ANPEP",
  "IL1RN","LPL","MFSD12"
)

# This is Figure 3B
pdf("DotPlot_Mono_Only_all.pdf", height = 10, width = 15, useDingbats=FALSE)
DotPlot(Mono.integrated, features = rev(markers.to.plot), 
        #cols = c("blue","red", "lightblue", "pink"), 
        dot.scale = 8, 
        #split.by = "both",
        assay = 'SCT'
) + RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")
dev.off()

########################################
##                                    ##
##        Pathway Combinations        ##
##                                    ##
########################################

# https://www.biostars.org/p/52101/
library(GO.db)
library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                   filters = 'go_id', values = 'GO:0007507', mart = ensembl)

#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
gene.data <- getBM(attributes=c('hgnc_symbol'),
                   filters = 'go', values = 'GO:0071345', mart = ensembl)
x <- c(gene.data$hgnc_symbol)

GOBPOFFSPRING[["GO:0071345"]]
t <- c(GOBPOFFSPRING[["GO:0071345"]], "GO:0071345")
gene.data <- getBM(attributes=c('hgnc_symbol'),
                   filters = 'go', values = t, mart = ensembl)
x <- c(gene.data$hgnc_symbol)

DefaultAssay(Mono.integrated) <- "RNA"
Mono.integrated_temp <- PercentageFeatureSet(Mono.integrated, 
                                        features = x, 
                                        col.name = "interferon", 
                                        assay = "RNA")

VlnPlot(Mono.integrated_temp, 
        features = c("interferon"), 
        split.by = "ultra", 
        idents = c("0","1","2"),
        pt.size = 0, ncol = 1)


########################################
##                                    ##
##      Subset HD and HC by Treat     ##
##                                    ##
########################################

Mono.integrated$mono_disease_weight <- paste(Mono.integrated$seurat_clusters,
                                             Mono.integrated$disease, 
                                             Mono.integrated$weight, sep = "_")
Mono.integrated$mono_disease_weight_treatments <- paste(Mono.integrated$seurat_clusters,
                                             Mono.integrated$disease, 
                                             Mono.integrated$weight, 
                                             Mono.integrated$first, 
                                             Mono.integrated$second, 
                                             sep = "_")

Idents(Mono.integrated) <- Mono.integrated$disease
HCHD <- subset(Mono.integrated, idents = c("HC","HD"))
table(Idents(HCHD))

ISGgenes <- c(
  "RSAD2","GBP5","CNP","BST2",
  "OAS3","OAS2","OAS1",
  "APOBEC3H","APOBEC3G","APOBEC3F","APOBEC3D","APOBEC3C","APOBEC3A",
  "ADAR",
  "DDX58","IFIH1",
  "IFIT5","IFIT3","IFIT2","IFIT1",
  "IFI6",
  "SSBP3","PARP12","OASL","NT5C3A",
  "MS4A4A","MAP3K14","IFI44L","EIF2AK2",
  "DDX60","DDIT4","CGAS",
  "ZC3HAV1",
  "TRIM5","MX1",
  "ISG20","IFI16",
  "MOV10",
  "NCOA7",
  "IFITM3","IFITM2","IFITM1"
  #"IRF1",
  #"IRF7", "TRIM25"
)

# POLY IC ONLY
Idents(HCHD) <- HCHD$treatments
HCHD_Poly <- subset(HCHD, idents = c("Basal_Basal","Poly_Basal"))
table(Idents(HCHD_Poly))

DefaultAssay(HCHD_Poly) <- "SCT"
Idents(HCHD_Poly) <- HCHD_Poly$mono_disease_weight
table(Idents(HCHD_Poly))
Idents(HCHD_Poly) <- factor(Idents(HCHD_Poly), 
                            levels = c("0_HC_Lean","0_HC_Obese","0_HD_Lean","0_HD_Obese",
                                       "1_HC_Lean","1_HC_Obese","1_HD_Lean","1_HD_Obese",
                                       "2_HC_Lean","2_HC_Obese","2_HD_Lean","2_HD_Obese",
                                       "3_HC_Lean","3_HC_Obese","3_HD_Lean","3_HD_Obese",
                                       "4_HC_Lean","4_HC_Obese","4_HD_Lean","4_HD_Obese",
                                       "5_HC_Lean","5_HC_Obese","5_HD_Lean","5_HD_Obese",
                                       "6_HC_Lean","6_HC_Obese","6_HD_Lean","6_HD_Obese",
                                       "7_HC_Lean","7_HC_Obese","7_HD_Lean","7_HD_Obese",
                                       "8_HC_Lean","8_HC_Obese","8_HD_Lean","8_HD_Obese",
                                       "9_HC_Lean","9_HC_Obese","9_HD_Lean","9_HD_Obese",
                                       "10_HC_Lean","10_HC_Obese","10_HD_Lean","10_HD_Obese",
                                       "11_HC_Lean","11_HC_Obese","11_HD_Lean","11_HD_Obese",
                                       "12_HC_Lean","12_HC_Obese","12_HD_Lean","12_HD_Obese",
                                       "13_HC_Lean","13_HC_Obese","13_HD_Lean","13_HD_Obese"))

DefaultAssay(HCHD_Poly) <- "RNA"
HCHD_Poly_temp <- PercentageFeatureSet(HCHD_Poly, 
                                             features = ISGgenes, 
                                             col.name = "ISG", 
                                             assay = "RNA")
HCHD_Poly_temp$ISG_tot <- HCHD_Poly_temp$ISG*HCHD_Poly_temp$nCount_RNA

VlnPlot(HCHD_Poly_temp, 
        features = c("ISG_tot"), 
        split.by = "first", 
        pt.size = 0, ncol = 1)


# This is Figure 3C, top

pdf("VlnPlot_Mono_Only_Poly_ISG.pdf", height = 3, width = 15, useDingbats=FALSE)
plots <- VlnPlot(HCHD_Poly_temp, 
                 features = c("ISG_tot"), 
                 split.by = "first", 
                 pt.size = 0, combine = FALSE )
plots <- lapply(
  X = plots,
  FUN = function(p) p 
  + ggplot2::theme(axis.title.y = element_blank())
  + ggplot2::scale_fill_manual(values = c('black', 'red'))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()


# LPS ONLY
Idents(HCHD) <- HCHD$treatments
HCHD_LPS <- subset(HCHD, idents = c("Basal_Basal","Basal_LPS"))
table(Idents(HCHD_LPS))

DefaultAssay(HCHD_LPS) <- "SCT"
Idents(HCHD_LPS) <- HCHD_LPS$mono_disease_weight
table(Idents(HCHD_LPS))
Idents(HCHD_LPS) <- factor(Idents(HCHD_LPS), 
                            levels = c("0_HC_Lean","0_HC_Obese","0_HD_Lean","0_HD_Obese",
                                       "1_HC_Lean","1_HC_Obese","1_HD_Lean","1_HD_Obese",
                                       "2_HC_Lean","2_HC_Obese","2_HD_Lean","2_HD_Obese",
                                       "3_HC_Lean","3_HC_Obese","3_HD_Lean","3_HD_Obese",
                                       "4_HC_Lean","4_HC_Obese","4_HD_Lean","4_HD_Obese",
                                       "5_HC_Lean","5_HC_Obese","5_HD_Lean","5_HD_Obese",
                                       "6_HC_Lean","6_HC_Obese","6_HD_Lean","6_HD_Obese",
                                       "7_HC_Lean","7_HC_Obese","7_HD_Lean","7_HD_Obese",
                                       "8_HC_Lean","8_HC_Obese","8_HD_Lean","8_HD_Obese",
                                       "9_HC_Lean","9_HC_Obese","9_HD_Lean","9_HD_Obese",
                                       "10_HC_Lean","10_HC_Obese","10_HD_Lean","10_HD_Obese",
                                       "11_HC_Lean","11_HC_Obese","11_HD_Lean","11_HD_Obese",
                                       "12_HC_Lean","12_HC_Obese","12_HD_Lean","12_HD_Obese",
                                       "13_HC_Lean","13_HC_Obese","13_HD_Lean","13_HD_Obese"))

DefaultAssay(HCHD_LPS) <- "RNA"
HCHD_LPS_temp <- PercentageFeatureSet(HCHD_LPS, 
                                       features = ISGgenes, 
                                       col.name = "ISG", 
                                       assay = "RNA")
HCHD_LPS_temp$ISG_tot <- HCHD_LPS_temp$ISG*HCHD_LPS_temp$nCount_RNA

VlnPlot(HCHD_LPS_temp, 
        features = c("ISG_tot"), 
        split.by = "second", 
        pt.size = 0, ncol = 1)


# This is Figure 3C, bottom

pdf("VlnPlot_Mono_Only_LPS_ISG.pdf", height = 3, width = 15, useDingbats=FALSE)
plots <- VlnPlot(HCHD_LPS_temp, 
                 features = c("ISG_tot"), 
                 split.by = "second", 
                 pt.size = 0, combine = FALSE  )
plots <- lapply(
  X = plots,
  FUN = function(p) p 
  + ggplot2::theme(axis.title.y = element_blank())
  + ggplot2::scale_fill_manual(values = c('black', 'red'))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()



########################################
##                                    ##
##         IFNG Pathway Figure        ##
##                                    ##
########################################

# This code gives you Figure 4D
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0034341 - IFNg
gene.data <- getBM(attributes=c('hgnc_symbol'),
                   filters = 'go', values = 'GO:0034341', mart = ensembl)
x <- c(gene.data$hgnc_symbol)
GOBPOFFSPRING[["GO:0034341"]]

t <- c("GO:0060330", "GO:0060331", "GO:0060332", "GO:0060333", "GO:0060334",
       "GO:0060335", "GO:0060336",
       "GO:0071346", "GO:0034341")
gene.data <- getBM(attributes=c('hgnc_symbol'),
                   filters = 'go', values = 'GO:0071346', mart = ensembl)
x <- c(gene.data$hgnc_symbol)
gene.data <- getBM(attributes=c('hgnc_symbol'),
                   filters = 'go', values = t, mart = ensembl)
x <- c(gene.data$hgnc_symbol)
x <- c("FCAR",     "CCL4L2",   "CCL3L3",   "CCL18",    "CCL3",     "CCL4",     "CCL23",    "HLA-DPA1",
       "CCL16",    "AIF1",     "CCL15",    "CCL5",     "CCL14",    "MRC1",     "EDN1",     "RAB20",    "SIRPA",   
       "KIF5B",    "ACOD1",    "DAPK1",    "AQP4",     "PDE12",    "RPS6KB1",  "ACTR3",    "STXBP4",   "CCL17",   
       "CCL22",    "MYO1C",    "CX3CL1",   "STXBP1",   "CCL19",    "CCL21",    "STX8",     "STX4",     "XCL2",    
       "CD58",     "XCL1",     "IL12RB1",  "SYNCRIP",  "ASS1",     "HLA-DQB2", "TLR2",     "WAS",      "RPL13A",  
       "TDGF1",    "CCL20",    "STXBP3",   "CLDN1",    "ACTG1",    "WNT5A",    "GBP6",     "CCL25",    "GSN",     
       "VIM",      "TLR4",     "ADAMTS13", "GAPDH",    "ZYX",      "CD47",     "IL12B",    "RAB43",    "CCL26",   
       "CCL24",    "DAPK3",    "SLC26A6",  "RAB7B",    "TLR3",     "FLNB",     "TNF",      "CDC42EP4", "GBP3",    
       "GBP1",     "GBP2",     "GBP7",     "GBP4",     "CIITA",    "IRF8",     "CCL7",     "CCL11",    "CCL8",    
       "CCL13",    "CCL1",     "CCL2",     "LGALS9",   "NOS2",     "VPS26B",   "ACTR2",    "CASP1",    "CDC42EP2",
       "GBP5",     "EPRS1",    "STAT1",    "CDC42",    "FASLG",    "VAMP3")
#"CCL3L1" - This one gene seems to mess the whole thing up??? 

Idents(Mono.integrated) <- Mono.integrated$cell_disease
DefaultAssay(Mono.integrated) <- "RNA"
Mono.integrated_temp <- PercentageFeatureSet(Mono.integrated, 
                                             features = x, 
                                             col.name = "ifng", 
                                             assay = "RNA")
Mono.integrated_temp$ifng_tot <- Mono.integrated_temp$ifng*Mono.integrated_temp$nCount_RNA

Idents(Mono.integrated_temp)<-factor(Idents(Mono.integrated_temp), 
                                    levels = c("CD14_Monocyte1_HC","CD14_Monocyte2_HC","CD16_Monocyte_HC",
                                               "CD14_Monocyte1_HD","CD14_Monocyte2_HD","CD16_Monocyte_HD"))

VlnPlot(Mono.integrated_temp, 
        features = c("ifng"), 
        split.by = "ultra", 
        pt.size = 0, ncol = 1)


# Two ways to make the same figure
# Using monocytes only
#1200x600
plots <- VlnPlot(Mono.integrated_temp, features = c("ifng_tot"), split.by = "ultra",
                 pt.size = 0, combine = FALSE, ncol = 1)
plots <- lapply(
  X = plots,
  FUN = function(p) p 
  + ggplot2::theme(axis.title.y = element_blank())
  + ggplot2::scale_fill_manual(values = rep(c('black', 'red','grey','pink'),times=4))
)
CombinePlots(plots = plots, ncol = 1)


# Using the entire dataset, but focusing on monocytes
Idents(pbmc.integrated) <- pbmc.integrated$cell_disease
DefaultAssay(pbmc.integrated) <- "RNA"
pbmc.integrated_temp <- PercentageFeatureSet(pbmc.integrated, 
                                             features = x, 
                                             col.name = "ifng", 
                                             assay = "RNA")
pbmc.integrated_temp$ifng_tot <- pbmc.integrated_temp$ifng*pbmc.integrated_temp$nCount_RNA

Idents(pbmc.integrated_temp)<-factor(Idents(pbmc.integrated_temp), 
                                     levels = c("CD14_Monocyte1_HC","CD14_Monocyte2_HC","CD16_Monocyte_HC",
                                                "CD14_Monocyte1_HD","CD14_Monocyte2_HD","CD16_Monocyte_HD"))
#1200x600
pdf("ViolinPlot_IFNGall.pdf", height = 6, width = 12, useDingbats=FALSE)
plots <- VlnPlot(pbmc.integrated_temp, features = c("ifng_tot"), split.by = "ultra",
                idents = c("CD14_Monocyte1_HC","CD14_Monocyte2_HC","CD16_Monocyte_HC",
                           "CD14_Monocyte1_HD","CD14_Monocyte2_HD","CD16_Monocyte_HD"),
                 pt.size = 0, combine = FALSE, ncol = 1)
plots <- lapply(
  X = plots,
  FUN = function(p) p 
  + ggplot2::theme(axis.title.y = element_blank())
  + ggplot2::scale_fill_manual(values = rep(c('black', 'red','grey','pink'),times=4))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()


# Dot plot
y <- c("ACOD1","ACTG1","ACTR2","ACTR3","ADAMTS13","AIF1",
       "CASP1",
       "CCL1","CCL2","CCL3","CCL3L3","CCL4","CCL4L2","CCL5","CCL7","CCL8",
       "CCL13","CCL15","CCL17","CCL18","CCL19","CCL20","CCL22","CCL23","CCL24",          
       "CDC42","CDC42EP4","CD47","CD58","CIITA", 
       "DAPK1","DAPK3",             
       "EDN1","EPRS1", 
       "FASLG","FCAR","FLNB",   
       "GAPDH", "GBP1","GBP2","GBP3","GBP4","GBP5","GBP7","GSN",  
       "HLA-DPA1","HLA-DQB2",
       "IL12B","IL12RB1","IRF8", 
       "KIF5B","LGALS9",  
       "MRC1", "MYO1C","PDE12",        
       "RAB7B","RAB20","RAB43","RPL13A","RPS6KB1",
       "SIRPA","SLC26A6","STAT1","STX4","STX8","STXBP1","STXBP3","STXBP4","SYNCRIP", 
       "TLR2","TLR3","TLR4","TNF", 
       "VAMP3","VIM","VPS26B",          
       "WAS","WNT5A",
       "XCL1","XCL2","ZYX"
       )
#Removed for low expression - "CCL16","CCL14","AQP4","CX3CL1","CCL21","ASS1","TDGF1","CLDN1","GBP6","CCL25","CCL26","CCL11","NOS2","CDC42EP2",

Idents(pbmc.integrated)<-pbmc.integrated$celltype # Categorize by celltype
MonoSub <- subset(pbmc.integrated, idents = c("CD14_Monocyte1","CD14_Monocyte2","CD16_Monocyte","DC")) # only monocytes and DC

Idents(MonoSub) <- MonoSub$cell_disease_weight_treatments

Idents(MonoSub)<-factor(Idents(MonoSub), 
                        levels = c("DC_HD_Obese_Poly_LPS","DC_HD_Obese_Poly_Basal","DC_HD_Obese_Basal_LPS","DC_HD_Obese_Basal_Basal",
                                   "DC_HD_Lean_Poly_LPS","DC_HD_Lean_Poly_Basal","DC_HD_Lean_Basal_LPS","DC_HD_Lean_Basal_Basal",
                                   "DC_HC_Obese_Poly_LPS","DC_HC_Obese_Poly_Basal","DC_HC_Obese_Basal_LPS","DC_HC_Obese_Basal_Basal",
                                   "DC_HC_Lean_Poly_LPS","DC_HC_Lean_Poly_Basal","DC_HC_Lean_Basal_LPS","DC_HC_Lean_Basal_Basal",
                                                
                                   "CD14_Monocyte1_HD_Obese_Poly_LPS","CD14_Monocyte1_HD_Obese_Poly_Basal","CD14_Monocyte1_HD_Obese_Basal_LPS","CD14_Monocyte1_HD_Obese_Basal_Basal",
                                   "CD14_Monocyte1_HD_Lean_Poly_LPS","CD14_Monocyte1_HD_Lean_Poly_Basal","CD14_Monocyte1_HD_Lean_Basal_LPS","CD14_Monocyte1_HD_Lean_Basal_Basal",
                                   "CD14_Monocyte1_HC_Obese_Poly_LPS","CD14_Monocyte1_HC_Obese_Poly_Basal","CD14_Monocyte1_HC_Obese_Basal_LPS","CD14_Monocyte1_HC_Obese_Basal_Basal",
                                   "CD14_Monocyte1_HC_Lean_Poly_LPS","CD14_Monocyte1_HC_Lean_Poly_Basal","CD14_Monocyte1_HC_Lean_Basal_LPS","CD14_Monocyte1_HC_Lean_Basal_Basal",
                                                
                                   "CD14_Monocyte2_HD_Obese_Poly_LPS","CD14_Monocyte2_HD_Obese_Poly_Basal","CD14_Monocyte2_HD_Obese_Basal_LPS","CD14_Monocyte2_HD_Obese_Basal_Basal",
                                   "CD14_Monocyte2_HD_Lean_Poly_LPS","CD14_Monocyte2_HD_Lean_Poly_Basal","CD14_Monocyte2_HD_Lean_Basal_LPS","CD14_Monocyte2_HD_Lean_Basal_Basal",
                                   "CD14_Monocyte2_HC_Obese_Poly_LPS","CD14_Monocyte2_HC_Obese_Poly_Basal","CD14_Monocyte2_HC_Obese_Basal_LPS","CD14_Monocyte2_HC_Obese_Basal_Basal",
                                   "CD14_Monocyte2_HC_Lean_Poly_LPS","CD14_Monocyte2_HC_Lean_Poly_Basal","CD14_Monocyte2_HC_Lean_Basal_LPS","CD14_Monocyte2_HC_Lean_Basal_Basal",

                                   "CD16_Monocyte_HD_Obese_Poly_LPS","CD16_Monocyte_HD_Obese_Poly_Basal","CD16_Monocyte_HD_Obese_Basal_LPS","CD16_Monocyte_HD_Obese_Basal_Basal",
                                   "CD16_Monocyte_HD_Lean_Poly_LPS","CD16_Monocyte_HD_Lean_Poly_Basal","CD16_Monocyte_HD_Lean_Basal_LPS","CD16_Monocyte_HD_Lean_Basal_Basal",
                                   "CD16_Monocyte_HC_Obese_Poly_LPS","CD16_Monocyte_HC_Obese_Poly_Basal","CD16_Monocyte_HC_Obese_Basal_LPS","CD16_Monocyte_HC_Obese_Basal_Basal",
                                   "CD16_Monocyte_HC_Lean_Poly_LPS","CD16_Monocyte_HC_Lean_Poly_Basal","CD16_Monocyte_HC_Lean_Basal_LPS","CD16_Monocyte_HC_Lean_Basal_Basal"
                                   ))
table(Idents(MonoSub))

# This is Supplemental Figure 4B
pdf("DotPlot_IFNGPathway_Mono1.pdf", height = 6, width = 16, useDingbats=FALSE)
DotPlot(MonoSub, features = y,
        idents = c("CD14_Monocyte1_HC_Lean_Basal_Basal","CD14_Monocyte1_HC_Lean_Poly_Basal",
                   "CD14_Monocyte1_HC_Lean_Basal_LPS","CD14_Monocyte1_HC_Lean_Poly_LPS",
                   "CD14_Monocyte1_HC_Obese_Basal_Basal","CD14_Monocyte1_HC_Obese_Poly_Basal",
                   "CD14_Monocyte1_HC_Obese_Basal_LPS","CD14_Monocyte1_HC_Obese_Poly_LPS",
                   "CD14_Monocyte1_HD_Lean_Basal_Basal","CD14_Monocyte1_HD_Lean_Poly_Basal",
                   "CD14_Monocyte1_HD_Lean_Basal_LPS","CD14_Monocyte1_HD_Lean_Poly_LPS",
                   "CD14_Monocyte1_HD_Obese_Basal_Basal","CD14_Monocyte1_HD_Obese_Poly_Basal",
                   "CD14_Monocyte1_HD_Obese_Basal_LPS","CD14_Monocyte1_HD_Obese_Poly_LPS"
                   ),
          dot.scale = 4, 
          assay = 'SCT'
) + RotatedAxis() + theme(axis.text.x = element_text(angle = 90))
dev.off()

# This is Supplemental Figure 4B
pdf("DotPlot_IFNGPathway_Mono2.pdf", height = 6, width = 16, useDingbats=FALSE)
DotPlot(MonoSub, features = y,
        idents = c("CD14_Monocyte2_HC_Lean_Basal_Basal","CD14_Monocyte2_HC_Lean_Poly_Basal",
                   "CD14_Monocyte2_HC_Lean_Basal_LPS","CD14_Monocyte2_HC_Lean_Poly_LPS",
                   "CD14_Monocyte2_HC_Obese_Basal_Basal","CD14_Monocyte2_HC_Obese_Poly_Basal",
                   "CD14_Monocyte2_HC_Obese_Basal_LPS","CD14_Monocyte2_HC_Obese_Poly_LPS",
                   "CD14_Monocyte2_HD_Lean_Basal_Basal","CD14_Monocyte2_HD_Lean_Poly_Basal",
                   "CD14_Monocyte2_HD_Lean_Basal_LPS","CD14_Monocyte2_HD_Lean_Poly_LPS",
                   "CD14_Monocyte2_HD_Obese_Basal_Basal","CD14_Monocyte2_HD_Obese_Poly_Basal",
                   "CD14_Monocyte2_HD_Obese_Basal_LPS","CD14_Monocyte2_HD_Obese_Poly_LPS"
        ),
        dot.scale = 4, 
        assay = 'SCT'
) + RotatedAxis() + theme(axis.text.x = element_text(angle = 90))
dev.off()

# This is Figure was not used
pdf("DotPlot_IFNGPathway_CD16Mono.pdf", height = 6, width = 16, useDingbats=FALSE)
DotPlot(MonoSub, features = y,
        idents = c("CD16_Monocyte_HC_Lean_Basal_Basal","CD16_Monocyte_HC_Lean_Poly_Basal",
                   "CD16_Monocyte_HC_Lean_Basal_LPS","CD16_Monocyte_HC_Lean_Poly_LPS",
                   "CD16_Monocyte_HC_Obese_Basal_Basal","CD16_Monocyte_HC_Obese_Poly_Basal",
                   "CD16_Monocyte_HC_Obese_Basal_LPS","CD16_Monocyte_HC_Obese_Poly_LPS",
                   "CD16_Monocyte_HD_Lean_Basal_Basal","CD16_Monocyte_HD_Lean_Poly_Basal",
                   "CD16_Monocyte_HD_Lean_Basal_LPS","CD16_Monocyte_HD_Lean_Poly_LPS",
                   "CD16_Monocyte_HD_Obese_Basal_Basal","CD16_Monocyte_HD_Obese_Poly_Basal",
                   "CD16_Monocyte_HD_Obese_Basal_LPS","CD16_Monocyte_HD_Obese_Poly_LPS"
        ),
        #          cols = c("black","black", "black", "black"), 
        dot.scale = 4, 
        #          split.by = "treatments",
        assay = 'SCT'
) + RotatedAxis() + theme(axis.text.x = element_text(angle = 90))
dev.off()

# This is Figure was not used
pdf("DotPlot_IFNGPathway_DC.pdf", height = 6, width = 16, useDingbats=FALSE)
DotPlot(MonoSub, features = y,
        idents = c("DC_HC_Lean_Basal_Basal","DC_HC_Lean_Poly_Basal",
                   "DC_HC_Lean_Basal_LPS","DC_HC_Lean_Poly_LPS",
                   "DC_HC_Obese_Basal_Basal","DC_HC_Obese_Poly_Basal",
                   "DC_HC_Obese_Basal_LPS","DC_HC_Obese_Poly_LPS",
                   "DC_HD_Lean_Basal_Basal","DC_HD_Lean_Poly_Basal",
                   "DC_HD_Lean_Basal_LPS","DC_HD_Lean_Poly_LPS",
                   "DC_HD_Obese_Basal_Basal","DC_HD_Obese_Poly_Basal",
                   "DC_HD_Obese_Basal_LPS","DC_HD_Obese_Poly_LPS"
        ),
        #          cols = c("black","black", "black", "black"), 
        dot.scale = 4, 
        #          split.by = "treatments",
        assay = 'SCT'
) + RotatedAxis() + theme(axis.text.x = element_text(angle = 90))
dev.off()












