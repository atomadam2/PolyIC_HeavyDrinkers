# PolyIC_HeavyDrinkers
This is a repository of BASH and R scripts used for my PolyIC in Heavy Drinkers paper. 

IMPORTANT NOTES:

1. This manuscript is currently unpublished, but has been submitted.
2. The methods used to generate this dataset are important to understanding what was done for the analysis. All samples were multiplexed using the Biolegend Totalseq protocol for Cell Hasing using the 10x v3.1 single-cell RNA-seq kits. In short, this is a homebrew way to multiplex different samples into one lane of a single-cell run. First you buy Biolegend antibodies that label all cells with a different barcode. Then you combine different samples that will be multiplexed, then run the single-cell. Then, on the molecular biology side of things, these barcodes are much smaller than fractionated cDNA, so you separate them out during the bead cleanups, clea the barcodes, then everything gets sequenced together (you only need to sequenced these barcode libraries at about 5% of the total sequencing reads).
3. For these experiments, multiplexed samples consisted of all drug treatments for a single patient. So one lane would include Basal, PolyIC, LPS, and PolyIC+LPS for a single patient.
4. On top of that, the libraries were also separately illumina barcoded and multiplexed for sequencing. So one sequencing run had multiple illumina barcodes which needed to be separated, and each fastq consisted of cells that came from different experiments that had to be separated.
5. As a result, the fastq files contain cells from different experiments mixed together.
6. Because NCBI and SRA do not handle single-cell data well (I've learned this the hard way), we have decided to only upload the BAM files that come out of cellranger. This is useful because the BAM files are all nicely separated by experiment and labeled in an intuitive way.
7. The pros are you can go straight from these BAM files into R to make a seurat object of the data, and never have to run cellranger.
8. The other pro is that separating the cells into individual experiments with cellranger is very very weird.
9. The cons are I have never actually used BAM files to create a Seurat object, so I currently have no code to help with this. 

For this repository, I will include all R scripts used to generate the figures in the paper starting from the complete Robject containing all the data. 


