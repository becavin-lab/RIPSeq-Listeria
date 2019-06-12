################################################################################
### R script to compare several conditions with the SARTools and edgeR packages
### Hugo Varet
### May 9th, 2016
### designed to be executed with SARTools 1.4.0
################################################################################

################################################################################
###                parameters: to be modified by the user                    ###
################################################################################
author <- "Christophe Becavin"                                # author of the statistical analysis/report

projectName <- "RIPSeq_RigI"                         # name of the project
targetFile <- paste("target_",projectName,".txt",sep = "")                           # path to the design/target file


featuresToRemove <- c("alignment_not_unique",        # names of the features to be removed
                      "ambiguous", "no_feature",     # (specific HTSeq-count information and rRNA for example)
                      "not_aligned", "too_low_aQual")# NULL if no feature to remove

varInt <- "group"                                    # factor of interest
condRef <- "CherryE"                                      # reference biological condition
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
                      rawDir=rawDirRIPSeqRigI,featuresToRemove=featuresToRemove,varInt=varInt,
                      condRef=condRef,batch=batch,alpha=alpha,pAdjustMethod=pAdjustMethod,
                      cpmCutoff=cpmCutoff,gene.selection=gene.selection,
                      normalizationMethod=normalizationMethod,colors=colors)



# loading target file
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# loading counts
counts <- loadCountData(target=target, rawDir=rawDirRIPSeqRigI, featuresToRemove=featuresToRemove)

# description plots
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)



# edgeR analysis
out.edgeR <- run.edgeR(counts=counts, target=target, varInt=varInt, condRef=condRef,
                       batch=batch, cpmCutoff=cpmCutoff, normalizationMethod=normalizationMethod,
                       pAdjustMethod=pAdjustMethod)

# MDS + clustering
exploreCounts(object=out.edgeR$dge, group=target[,varInt], gene.selection=gene.selection, col=colors)

# summary of the analysis (boxplots, dispersions, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.edgeR(out.edgeR, group=target[,varInt], counts=counts, alpha=alpha, col=colors)

# generating HTML report
writeReport.edgeR(target=target, counts=counts, out.edgeR=out.edgeR, summaryResults=summaryResults,
                  majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
                  targetFile=targetFile, rawDir=rawDirRIPSeqRigI, featuresToRemove=featuresToRemove, varInt=varInt,
                  condRef=condRef, batch=batch, alpha=alpha, pAdjustMethod=pAdjustMethod, colors=colors,
                  gene.selection=gene.selection, normalizationMethod=normalizationMethod, cpmCutoff=cpmCutoff)

# move files
image_folder = paste(dirProject,"/figures_",projectName,"/",sep="")
if (!dir.exists(image_folder)){
  file.rename(from="figures/",to=image_folder)
}else{
  unlink(image_folder, recursive=TRUE)
  file.rename(from="figures/",to=image_folder)
}
table_folder = paste(dirProject,"/tables_",projectName,"/",sep="")
if (!dir.exists(table_folder)){
  file.rename(from="tables/",to=table_folder)
}else{
  unlink(image_folder, recursive=TRUE)
  file.rename(from="tables/",to=table_folder)
}

file.rename(from=paste(projectName,"_report.html",sep=""),to=paste(dirProject,"/",projectName,"_report.html",sep=""))

# create bed file for Circos plot
diff_matrix = read.delim(file = paste(table_folder,"RIGIEvsCherryE.up.txt",sep=""))
annot = read.delim(file="NC_003210_Annot.txt")
diff_annot = merge(diff_matrix,annot,by.x="Id",by.y = "Locustag", all.x=TRUE,all.y=F)
write.table(diff_annot, file="../Results/RIGIEvsCherryE.txt", sep="\t", quote=FALSE, row.names = F, col.names = T)
diff_annot$rigi = "rigi"
diff_annot$start = 0
diff_annot$finish = 1500
diff_annot$chr = "chr"
bed_file = diff_annot[c("rigi","start","finish","chr","begin","end")]
write.table(bed_file, file="../Results/RIGIEvsCherryE.bed", sep="\t", quote=FALSE, row.names = F, col.names = F)

highlight = diff_annot[c("chr","begin","end")]
highlight$color = "fill_color=grey"
colnames(highlight) = c("chr","2360731","2402842","fill_color=blue")
write.table(highlight, file="../Results/RIGIEvsCherryE_highlight.conf", sep="\t", quote=FALSE, row.names = F, col.names = T)

            
            