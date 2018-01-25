#!/usr/bin/env Rscript

### RScript for DiffBind analysis of ATAC data. 
# This is designed to create all sample contrasts as permitted by the sample sheet.
# Currently uses DESEQ2 at the backend, with blocking model for batches. Could contrast with EdgeR.

# Package loading
library(DiffBind)
library(rtracklayer)

# Argument loading - Should direct to samples file
args = commandArgs(trailingOnly=TRUE)  # Command allows us to specify file using RScript on cmd line
SampleFile <- args[1]

# Create pathways and load required data
dir.create(file.path("./DiffBind_Output"), showWarnings = FALSE)
dir.create(file.path("./DiffBind_Results"))
samples <- read.csv(SampleFile)

##### Setup #####

# Import data using samples sheet (required!)
expData <- dba(sampleSheet = samples)

# Plot sample correlations, based on peak score
dev.off()
pdf("DiffBind_Output/CorrelationHeatmap_Occupancy-PeakCallerScore.pdf",8,8)
  dba.plotHeatmap(expData)
dev.off()

# Integrate count data for each peak, location based on sampleSheet.
# Main parameter to adjust is summits - i.e, the read pileup!
# If summits specified, peaks re-centered around consensus summit, and the value used to widen the range by, e.g, 250 bp's
expData <- dba.count(expData, summits=250)

# Plot sample correlations, based on read count
dev.off()
pdf("DiffBind_Output/CorrelationHeatmap_Affinity-ReadCount.pdf",8,8)
  dba.plotHeatmap(expData)
dev.off()

##### QC Graphs #####
pdf("DiffBind_Output/PCA_Plots.pdf",8,8)
  dba.plotPCA(expData, DBA_TISSUE)
  dba.plotPCA(expData, DBA_FACTOR)
  dba.plotPCA(expData, DBA_CONDITION)
  dba.plotPCA(expData, DBA_TREATMENT)
dev.off()

##### Contrast establishment - Factor flagging #####
# This works by checking how many different unique names are entered in each Sample sheet column. 
# 1 means no contrast to perform, 2 means a pairwise contrast, >2 means a complex analysis is required.
# FACTOR is reserved for batch conditions, to use for DESEQ2 blocking

Comps <- 0

if(length(unique(samples$Tissue)) == 2){
  FLAG_TISSUE <- TRUE
  Comps <- Comps + 1
}else if(length(unique(samples$Tissue)) > 2){
  FLAG_TISSUE <- "Complex"
  Tissue_min <- min(table(samples$Tissue))
  Comps <- Comps + 1
}else{FLAG_TISSUE <- FALSE}

if(length(unique(samples$Factor)) == 2){
  FLAG_FACTOR <- TRUE
  Comps <- Comps + 1
}else if(length(unique(samples$Factor)) > 2){
  FLAG_FACTOR <- "Complex"
  Factor_min <- min(table(samples$Factor))
  Comps <- Comps + 1
}else{FLAG_FACTOR <- FALSE}

if(length(unique(samples$Condition)) == 2){
  FLAG_CONDITION <- TRUE
  Comps <- Comps + 1
}else if(length(unique(samples$Condition)) > 2){
  FLAG_CONDITION <- "Complex"
  Condition_min <- min(table(samples$Condition))
  Comps <- Comps + 1
}else{FLAG_CONDITION <- FALSE}

if(length(unique(samples$Treatment)) == 2){
  FLAG_TREATMENT <- TRUE
  Comps <- Comps + 1
}else if(length(unique(samples$Treatment)) > 2){
  FLAG_TREATMENT <- "Complex"
  Treatment_min <- min(table(samples$Treatment))
  Comps <- Comps + 1
}else{FLAG_TREATMENT <- FALSE}

Flags <- data.frame(Flag = c("Tissue", "Factor", "Condition", "Treatment"), Status = c(FLAG_TISSUE, FLAG_FACTOR, FLAG_CONDITION, FLAG_TREATMENT))

##### Contrast specification ####

if(Comps == 0){
  print("No comparisons to be made")
  quit(save = "no")
}else{
  listComps <- vector("list", 4)
  names(listComps) <- c("Tissue", "Factor", "Condition", "Treatment")
}

if(FLAG_TISSUE == TRUE){
  listComps[[1]] <- dba.contrast(expData, categories=DBA_TISSUE, block=DBA_FACTOR)
}else if(FLAG_TISSUE == "Complex"){
  listComps[[1]] <- dba.contrast(expData, minMembers = Tissue_min, categories = DBA_TISSUE, block=DBA_FACTOR)
}else{listComps[[1]] <- FALSE}

if(FLAG_CONDITION == TRUE){
  listComps[[3]] <- dba.contrast(expData, categories=DBA_CONDITION, block=DBA_FACTOR)
}else if(FLAG_CONDITION == "Complex"){
  listComps[[3]] <- dba.contrast(expData, minMembers = Condition_min, categories = DBA_CONDITION, block=DBA_FACTOR)
}else{listComps[[3]] <- FALSE}

if(FLAG_TREATMENT == TRUE){
  listComps[[4]] <- dba.contrast(expData, categories=DBA_TREATMENT, block=DBA_FACTOR)
}else if(FLAG_TREATMENT == "Complex"){
  listComps[[4]] <- dba.contrast(expData, minMembers = Treatment_min, categories = DBA_TREATMENT, block=DBA_FACTOR)
}else{listComps[[4]] <- FALSE}

##### Performing Analysis #####
mapply(function(contrast, names){
  
  if(class(contrast) != "DBA"){
    return(FALSE)
    
  }else{
    dirOut <- paste("./DiffBind_Output/", names, sep="")
    dir.create(file.path(dirOut), showWarnings = FALSE)
    res <- dba.analyze(contrast, method=DBA_DESEQ2_BLOCK)
    comp_num <- length(res$contrasts)
    
    # plot HMs
    for(i in 1:comp_num){
      pdf(paste(dirOut, "/DiffCorrelationHM_Contrast_",res$contrasts[[i]]$name1,"_vs_",res$contrasts[[i]]$name2,".pdf", sep=""), 8,8)
        dba.plotHeatmap(res, contrast = i)
      dev.off()
    }
    
    # plot PCA
    for(i in 1:comp_num){
      pdf(paste(dirOut, "/Diff_PCA_", res$contrasts[[i]]$name1, "_vs_", res$contrasts[[i]]$name2,".pdf", sep=""), 8,8)
        dba.plotPCA(res, DBA_TISSUE, contrast=i)
        dba.plotPCA(res, DBA_FACTOR, contrast=i)
        dba.plotPCA(res, DBA_CONDITION, contrast=i)
        dba.plotPCA(res, DBA_TREATMENT, contrast=i)
      dev.off()
    }
    
    # plot MA
    pdf(paste(dirOut,"/MA_Plots.pdf", sep=""),8,8)
      for(i in 1:comp_num){
        dba.plotMA(res, contrast=i)
      }
    dev.off()
    
    # Retrieve diff peaks
    for(i in 1:comp_num){
      res.DB <- dba.report(res, contrast=i, method=DBA_DESEQ2_BLOCK)
      res.DB <- data.frame(PeakID = c(1:length(res.DB)), as.data.frame(res.DB))
      res.DB <- res.DB[,c(2:4,1,5:ncol(res.DB))]
      write.table(res.DB, paste(dirOut, "/DifferentiallyBound_Contrast_", res$contrasts[[i]]$name1, "_vs_", res$contrasts[[i]]$name2, ".bed", sep=""), quote=F, sep="\t", row.names=F)
      write.table(res.DB, paste("./DiffBind_Results/DifferentiallyBound_Contrast_", res$contrasts[[i]]$name1, "_vs_", res$contrasts[[i]]$name2, ".bed", sep=""), quote=F, sep="\t", row.names=F)
    }
    
    
    
    
  }
}, contrast = listComps, names = names(listComps))

write.table(Flags, "CompFlags.txt",sep = "\t",quote = F,row.names = F)
