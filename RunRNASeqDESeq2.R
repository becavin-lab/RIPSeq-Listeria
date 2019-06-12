################################################################################
### R script to compare several conditions with the SARTools and DESeq2 packages
### Hugo Varet
### April 20th, 2015
### designed to be executed with SARTools 1.1.0
################################################################################

################################################################################
###                parameters: to be modified by the user                    ###
################################################################################

projectName <- "RNASeq"                         # name of the project
targetFile <- paste("target_","RNASeq",".txt",sep = "")                           # path to the design/target file

featuresToRemove <- c("alignment_not_unique",        # names of the features to be removed
                      "ambiguous", "no_feature",     # (specific HTSeq-count information and rRNA for example)
                      "not_aligned", "too_low_aQual",
                      "N_noFeature","N_multimapping","N_unmapped","N_ambiguous")

varInt <- "group"                                    # factor of interest
condRef <- "WT"                                      # reference biological condition
batch <- NULL                                # blocking factor: NULL (default) or "batch" for example

fitType <- "parametric"                              # mean-variance relationship: "parametric" (default) or "local"
cooksCutoff <- TRUE                                  # TRUE/FALSE to perform the outliers detection (default is TRUE)
independentFiltering <- TRUE                         # TRUE/FALSE to perform independent filtering (default is TRUE)
alpha <- 0.05                                        # threshold of statistical significance
pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"

typeTrans <- "VST"                                   # transformation for PCA/clustering: "VST" or "rlog"
locfunc <- "median"                                  # "median" (default) or "shorth" to estimate the size factors

colors <- c("dodgerblue","dodgerblue4","firebrick1","MediumVioletRed","tan3","SpringGreen")

################################################################################
###                             running script                               ###
################################################################################

# checking parameters
checkParameters.DESeq2(projectName=projectName,author=author,targetFile=targetFile,
                       rawDir=rawDirRNASeq,featuresToRemove=featuresToRemove,varInt=varInt,
                       condRef=condRef,batch=batch,fitType=fitType,cooksCutoff=cooksCutoff,
                       independentFiltering=independentFiltering,alpha=alpha,pAdjustMethod=pAdjustMethod,
                       typeTrans=typeTrans,locfunc=locfunc,colors=colors)

# loading target file
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# loading counts
counts <- loadCountData(target=target, rawDir=rawDirRNASeq, featuresToRemove=featuresToRemove)

# Load annotation description
print("Load Genes")


# description plots
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# analysis with DESeq2
out.DESeq2 <- run.DESeq2(counts=counts, target=target, varInt=varInt, batch=batch,
                         locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                         cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)

PCAPlot <- function (counts.trans, group, n = min(500, nrow(counts.trans)),
                     col = c("dodgerblue","dodgerblue4","firebrick1","MediumVioletRed","tan3","SpringGreen"),
                     outfile = TRUE)
{
  rv = apply(counts.trans, 1, var, na.rm = TRUE)
  pca = prcomp(t(counts.trans[order(rv, decreasing = TRUE),
                              ][1:n, ]))
  prp <- pca$sdev^2 * 100/sum(pca$sdev^2)
  prp <- round(prp[1:3], 2)
  if (outfile)
    png(filename = "figures/PCA.png", width = 1800 * 2, height = 1800,
        res = 300)
  par(mfrow = c(1, 2))
  abs = range(pca$x[, 1])
  abs = abs(abs[2] - abs[1])/25
  ord = range(pca$x[, 2])
  ord = abs(ord[2] - ord[1])/25
  plot(pca$x[, 1], pca$x[, 2], las = 1, cex = 2, pch = 16,
       col = col[as.integer(group)], xlab = paste0("PC1 (",
                                                   prp[1], "%)"), ylab = paste0("PC2 (", prp[2], "%)"),
       main = "Principal Component Analysis - Axes 1 and 2")
  abline(h = 0, v = 0, lty = 2, col = "lightgray")
  text(pca$x[, 1] - ifelse(pca$x[, 1] > 0, abs, -abs), pca$x[,
                                                             2] - ifelse(pca$x[, 2] > 0, ord, -ord), colnames(counts.trans),
       col = col[as.integer(group)])
  abs = range(pca$x[, 1])
  abs = abs(abs[2] - abs[1])/25
  ord = range(pca$x[, 3])
  ord = abs(ord[2] - ord[1])/25
  plot(pca$x[, 1], pca$x[, 3], las = 1, cex = 2, pch = 16,
       col = col[as.integer(group)], xlab = paste0("PC1 (",
                                                   prp[1], "%)"), ylab = paste0("PC3 (", prp[3], "%)"),
       main = "Principal Component Analysis - Axes 1 and 3")
  abline(h = 0, v = 0, lty = 2, col = "lightgray")
  text(pca$x[, 1] - ifelse(pca$x[, 1] > 0, abs, -abs), pca$x[,3] - ifelse(pca$x[, 3] > 0, ord, -ord), colnames(counts.trans),
       col = col[as.integer(group)])
  if (outfile)
    dev.off()
  return(invisible(pca$x))
}

exploreCounts <- function (object, group, typeTrans = "VST", gene.selection = "pairwise",
                           col = c("dodgerblue","dodgerblue4","firebrick1","MediumVioletRed","tan3","SpringGreen"))
{
  if (class(object) == "DESeqDataSet") {
    if (typeTrans == "VST")
      counts.trans <- assay(varianceStabilizingTransformation(object))
    else counts.trans <- assay(rlogTransformation(object))
    PCAPlot(counts.trans = counts.trans, group = group, col = col)
    clusterPlot(counts.trans = counts.trans, group = group)
  }
  else if (class(object) == "DGEList") {
    MDSPlot(dge = object, group = group, col = col, gene.selection = gene.selection)
    clusterPlot(counts.trans = cpm(object, prior.count = 2,
                                   log = TRUE), group = group)
  }
  else {
    stop("The object is not a DESeqDataSet nor a DGEList")
  }
}

exportResults.DESeq2 <- function (out.DESeq2, group, alpha = 0.05, export = TRUE) 
{
  genes <- read.delim(annotFile,header = TRUE)
  dds <- out.DESeq2$dds
  results <- out.DESeq2$results
  counts <- data.frame(Id = rownames(counts(dds)), counts(dds), 
                       round(counts(dds, normalized = TRUE)))
  colnames(counts) <- c("Id", colnames(counts(dds)), paste0("norm.", 
                                                            colnames(counts(dds))))
  bm <- data.frame(Id = rownames(results[[1]]), baseMean = round(results[[1]][, 
                                                                              "baseMean"], 2))
  base <- merge(counts, bm, by = "Id", all = TRUE)
  tmp <- base[, paste("norm", colnames(counts(dds)), sep = ".")]
  for (cond in levels(group)) {
    base[, cond] <- round(apply(as.data.frame(tmp[, group == 
                                                    cond]), 1, mean), 0)
  }
  complete <- list()
  for (name in names(results)) {
    complete.name <- base
    res.name <- data.frame(Id = rownames(results[[name]]), 
                           FoldChange = round(2^(results[[name]][, "log2FoldChange"]), 
                                              3), log2FoldChange = round(results[[name]][, 
                                                                                         "log2FoldChange"], 3), pvalue = results[[name]][, 
                                                                                                                                         "pvalue"], padj = results[[name]][, "padj"])
    complete.name <- merge(complete.name, res.name, by = "Id", 
                           all = TRUE)
    mcols.add <- data.frame(Id = rownames(counts(dds)), dispGeneEst = round(mcols(dds)$dispGeneEst, 
                                                                            4), dispFit = round(mcols(dds)$dispFit, 4), dispMAP = round(mcols(dds)$dispMAP, 
                                                                                                                                        4), dispersion = round(mcols(dds)$dispersion, 4), 
                            betaConv = mcols(dds)$betaConv, maxCooks = round(mcols(dds)$maxCooks, 
                                                                             4))
    complete.name <- merge(complete.name, mcols.add, by = "Id", 
                           all = TRUE)
    complete[[name]] <- complete.name
    complete.name <- merge(complete.name,genes,by.x="Id",by.y="Locustag",all.x=TRUE)
    #header_filter =c("Id","naive","mem","samhd1","FoldChange","log2FoldChange","pvalue","padj")
    #header_filter <- c(header_filter,colnames(genes))
    #header_filter = header_filter[!header_filter=="gene_id"]
    #print(colnames(complete.name))
    #print(header_filter)
    if (export) {
      up.name <- complete.name[which(complete.name$padj <= 
                                       alpha & complete.name$betaConv & complete.name$log2FoldChange >= 
                                       0), ]
      up.name <- up.name[order(up.name$padj), ]
      down.name <- complete.name[which(complete.name$padj <= 
                                         alpha & complete.name$betaConv & complete.name$log2FoldChange <= 
                                         0), ]
      down.name <- down.name[order(down.name$padj), ]
      name <- gsub("_", "", name)
      #write.table(complete.name[,header_filter], file = paste0("tables/", 
      #                                         name, ".complete.txt"), sep = "\t", row.names = FALSE, 
      #            dec = ".", quote = FALSE)
      write.table(complete.name, file = paste0("tables/", 
                                               name, ".complete.txt"), sep = "\t", row.names = FALSE, 
                  dec = ".", quote = FALSE)
      write.table(up.name, file = paste0("tables/", name, 
                                         ".up.txt"), row.names = FALSE, sep = "\t", dec = ".", 
                  quote = FALSE)
      write.table(down.name, file = paste0("tables/", name, 
                                           ".down.txt"), row.names = FALSE, sep = "\t", 
                  dec = ".", quote = FALSE)
    }
  }
  return(complete)
}

summarizeResults.DESeq2 <- function (out.DESeq2, group, independentFiltering = TRUE, cooksCutoff = TRUE, 
                                     alpha = 0.05, col = c("lightblue", "orange", "MediumVioletRed", 
                                                           "SpringGreen")) 
{
  if (!I("figures" %in% dir())) 
    dir.create("figures", showWarnings = FALSE)
  if (!I("tables" %in% dir())) 
    dir.create("tables", showWarnings = FALSE)
  dds <- out.DESeq2$dds
  results <- out.DESeq2$results
  diagSizeFactorsPlots(dds = dds, group = group, col = col)
  countsBoxplots(dds, group = group, col = col)
  dispersionsPlot(dds = dds)
  if (independentFiltering) {
    tabIndepFiltering <- tabIndepFiltering(results)
    cat("Number of features discarded by the independent filtering:\n")
    print(tabIndepFiltering, quote = FALSE)
  }
  else {
    tabIndepFiltering <- NULL
  }
  complete <- exportResults.DESeq2(out.DESeq2, group = group, 
                                   alpha = alpha)
  nDiffTotal <- nDiffTotal(complete = complete, alpha = alpha)
  cat("\nNumber of features down/up and total:\n")
  print(nDiffTotal, quote = FALSE)
  rawpHist(complete = complete)
  MAPlot(complete = complete, alpha = alpha)
  volcanoPlot(complete = complete, alpha = alpha)
  return(list(complete = complete, tabIndepFiltering = tabIndepFiltering, 
              nDiffTotal = nDiffTotal))
}


# PCA + clustering
exploreCounts(object=out.DESeq2$dds, group=target[,varInt], typeTrans=typeTrans, col=colors)

# summary of the analysis (boxplots, dispersions, diag size factors, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.DESeq2(out.DESeq2, group=target[,varInt], col=colors,
                                          independentFiltering=independentFiltering,
                                          cooksCutoff=cooksCutoff, alpha=alpha)

# save image of the R session
#save.image(file=paste0(projectName, ".RData"))

# generating HTML report
writeReport.DESeq2(target=target, counts=counts, out.DESeq2=out.DESeq2, summaryResults=summaryResults,
                   majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
                   targetFile=targetFile, rawDir=rawDirRNASeq, featuresToRemove=featuresToRemove, varInt=varInt,
                   condRef=condRef, batch=batch, fitType=fitType, cooksCutoff=cooksCutoff,
                   independentFiltering=independentFiltering, alpha=alpha, pAdjustMethod=pAdjustMethod,
                   typeTrans=typeTrans, locfunc=locfunc, colors=colors)

# move files
image_folder = paste(dirProject,"figures_",projectName,"/",sep="")
file.rename(from="figures/",to=image_folder)
if (!dir.exists(image_folder)){
  file.rename(from="figures/",to=image_folder)
}else{
  unlink(image_folder, recursive=TRUE)
  file.rename(from="figures/",to=image_folder)
}
table_folder = paste(dirProject,"tables_",projectName,"/",sep="")
if (!dir.exists(table_folder)){
  file.rename(from="tables/",to=table_folder)
}else{
  unlink(image_folder, recursive=TRUE)
  file.rename(from="tables/",to=table_folder)
}

file.rename(from="RNASeq_report.html",to=paste(dirProject,projectName,"_report.html",sep=""))
