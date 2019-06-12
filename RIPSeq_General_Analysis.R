################################################################################
### R script to analyse RIPSeq and RNASeq files with the SARTools and edgeR packages
### Christophe Becavin
### December 12d, 2018
### designed to be executed with SARTools 1.4.0
################################################################################

################################################################################
###                parameters: to be modified by the user                    ###
################################################################################
rm(list=ls())                                        # remove all the objects from the R session
library(SARTools)

workDir <- "~/RNABindingProtein/ripseq-listeria"
setwd(workDir)      # working directory for the R session

projectName <- "RIPSeq_Zea"                         # name of the project
author <- "Christophe Becavin"                                # author of the statistical analysis/report

dirProject <- "../Results/"
if (!dir.exists(dirProject)){
  dir.create(dirProject)
}

targetFile <- "target_RIPSeq.txt"                           # path to the design/target file

annotFile = "NC_003210_Annot.txt"

rawDirRIPSeq <- "../Expression"                                      # path to the directory containing raw counts files
rawDirRNASeq <- "../Expression"                                      # path to the directory containing raw counts files
rawDirRIPSeqRigI <- "../Expression"                                      # path to the directory containing raw counts files


featuresToRemove <- c("alignment_not_unique",        # names of the features to be removed
                      "ambiguous", "no_feature",     # (specific HTSeq-count information and rRNA for example)
                      "not_aligned", "too_low_aQual")# NULL if no feature to remove

varInt <- "group"                                    # factor of interest
condRef <- "Fraction_BHI"                                      # reference biological condition
batch <- NULL                                        # blocking factor: NULL (default) or "batch" for example

alpha <- 0.05                                        # threshold of statistical significance
pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"

cpmCutoff <- 1                                       # counts-per-million cut-off to filter low counts
gene.selection <- "pairwise"                         # selection of the features in MDSPlot
normalizationMethod <- "upperquartile"                         # normalization method: "TMM" (default), "RLE" (DESeq) or "upperquartile"

colors <- c("dodgerblue","firebrick1",               # vector of colors of each biological condition on the plots
            "MediumVioletRed","SpringGreen")

################################################################################
###                             running script                               ###
################################################################################

# checking parameters
checkParameters.edgeR(projectName=projectName,author=author,targetFile=targetFile,
                      rawDir=rawDirRIPSeq,featuresToRemove=featuresToRemove,varInt=varInt,
                      condRef=condRef,batch=batch,alpha=alpha,pAdjustMethod=pAdjustMethod,
                      cpmCutoff=cpmCutoff,gene.selection=gene.selection,
                      normalizationMethod=normalizationMethod,colors=colors)

loadTargetFileModified <- function (targetFile, varInt, condRef, batch)
{
  target <- read.table(targetFile, header = TRUE, sep = "\t")
  if (!I(varInt %in% names(target)))
    stop(paste("The factor of interest", varInt, "is not in the target file"))
  if (!is.null(batch) && !I(batch %in% names(target)))
    stop(paste("The batch effect", batch, "is not in the target file"))
  target[, varInt] <- as.factor(target[, varInt])
  if (!I(condRef %in% as.character(target[, varInt])))
    stop(paste("The reference level", condRef, "is not a level of the factor of interest"))
  target[, varInt] <- relevel(target[, varInt], ref = condRef)
  target <- target[order(target[, varInt]), ]
  rownames(target) <- as.character(target[, 1])
  cat("Target file:\n")
  print(target)
  return(target)
}

# loading target file
target <- loadTargetFileModified(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# loading counts
counts <- loadCountData(target=target, rawDir=rawDirRIPSeq, featuresToRemove=featuresToRemove)

# description plots
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

run.edgdR.modified <- function (counts, target, varInt, condRef, batch = NULL, cpmCutoff = 1,
                                normalizationMethod = "TMM", pAdjustMethod = "BH", ...)
{
  minReplicates <- min(table(target[, varInt]))
  fcounts <- counts[rowSums(cpm(counts) >= cpmCutoff) >= minReplicates,
                    ]
  cat("Number of features discarded by the filtering:\n")
  cat(nrow(counts) - nrow(fcounts), "\n")
  design <- formula(paste("~", ifelse(!is.null(batch), paste(batch,
                                                             "+"), ""), varInt))
  dge <- DGEList(counts = fcounts, remove.zeros = TRUE)
  dge$design <- model.matrix(design, data = target)
  cat("\nDesign of the statistical model:\n")
  cat(paste(as.character(design), collapse = " "), "\n")
  dge <- calcNormFactors(dge, method = normalizationMethod)

  cat("\nNormalization factors:\n")
  print(dge$samples$norm.factors)

  return(list(dge = dge, results = results))
}

# Calc normalization factor
out.edgeR <- run.edgdR.modified(counts=counts, target=target, varInt=varInt, condRef=condRef,
                       batch=batch, cpmCutoff=cpmCutoff, normalizationMethod=normalizationMethod,
                       pAdjustMethod=pAdjustMethod)

# TMM normalization of counts
tmm <- out.edgeR$dge$samples$norm.factors
N <- colSums(out.edgeR$dge$counts)
f <- tmm * N/mean(tmm * N)
norm.counts <- scale(out.edgeR$dge$counts, center = FALSE, scale = f)
norm.counts <- removeNull(norm.counts)
norm.counts <- log2(norm.counts + 1)
norm.result = as.data.frame(norm.counts)

# calculate fold changes by substracting columns
norm.result$IP_LMO_BHI_C1_logFC = norm.result$IP_LMO_BHI_1 - norm.result$IP_PI_BHI_1
norm.result$IP_LMO_BHI_C2_logFC = norm.result$IP_LMO_BHI_2 - norm.result$IP_PI_BHI_2
norm.result$IP_LMO_BHI_Pooled_logFC = apply(norm.result[,c("IP_LMO_BHI_C1_logFC","IP_LMO_BHI_C2_logFC")],1,max)
norm.result$IP_LMO_Extract_FC = norm.result$IP_LMO_Extract_1 - norm.result$IP_PI_Extract_1

# Add annotation
annot <- read.delim("NC_003210_Annot.txt", row.names = 1)
annot$StartCodon = NULL
annot$SD.energy = NULL
norm.result <- merge(norm.result,annot,by="row.names",all.x=TRUE)
row.names(norm.result) <- norm.result$Row.names
write.table(norm.result,file=paste(dirProject,"RNABinding_All.xls",sep=""),sep = "\t",col.names = T, row.names = F, quote = F)

# move files
image_folder = paste(dirProject,"/figures",projectName,"",sep="")
if (!dir.exists(image_folder)){
  file.rename(from="figures/",to=image_folder)
}

# Heatmap of the paper were done with RNABinding_All.xls and in-house scripts from BACNET
# https://gitlab.pasteur.fr/bacnet/Bacnet-public

# Run RNASeq analysis
# long = AlienTrimmer with 45bp cutoff
source(file = "RunRNASeqDESeq2.R")

# Run RIPSeq of RIGI
source(file = "RunRIPSeqRigIEdgeR.r")
#  do circos plot of this RIPSEq RIGI
# "C:\Users\becav\OneDrive\Documents\RNABindingProtein\Circos\circos-0.69-6/bin/circos -conf C:\Users\becav\OneDrive\Documents\RNABindingProtein\ripseq-listeria\rigi_pulldown.conf"

# fisher test for RIP-Seq vs RNAseq
contingency_table = data.frame(c(12,280),c(25,2830))
fisher.test(contingency_table, alternative = "greater")

# Combine with Tiling_EGDe
norm.result <- read.delim("../Results/RNABinding_All.xls", row.names = 1)
norm.result$Row.names = NULL
keep_columns <- c("EGDe_030510.Tiling.","EGDe_Stat.Tiling.","IP_LMO_BHI_C1_logFC","IP_LMO_BHI_C2_logFC","IP_LMO_Extract_FC","IP_LMO_BHI_Pooled_logFC")
tiling_egde <- read.delim("Tiling_EGDe.txt", row.names = 1)
tiling_RipSeq <- merge(norm.result,tiling_egde,by="row.names")
row.names(tiling_RipSeq) <- tiling_RipSeq$Row.names
tiling_RipSeq$Locus_tag = row.names(tiling_RipSeq)
keep_columns <- c("Locus_tag","EGDe_030510.Tiling.","EGDe_Stat.Tiling.","IP_LMO_BHI_C1_logFC","IP_LMO_BHI_C2_logFC","IP_LMO_BHI_Pooled_logFC","IP_LMO_Extract_FC"
                  ,"begin","end","length","Strand","Info","Note")
tiling_RipSeq <- tiling_RipSeq[,keep_columns]
write.table(tiling_RipSeq,file="../Results/Tiling_vs_Enrich.xls",sep = "\t",col.names = T, row.names = F, quote = F)

