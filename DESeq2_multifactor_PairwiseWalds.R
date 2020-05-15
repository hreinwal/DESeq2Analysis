### DESeq2 RNA Seq Analysis Script for multiple factors / treatments
# Author: Hannes Reinwald
# Contact: hannes.reinwald@ime.fraunhofer.de

# README ####
# This script is desgined to run DESeq2 normalization and statistical testing on RNAseq experiments
# in Ecotox testing with multiple concentrations. This script will automatically adopt to the k numbers
# of samples and N numbers of conditions (High,Mid,Low or Treat1, Treat2, Treat3, ... , Control) and will 
# generate Log2FC values with respect to control groups.

# Requiered Input: 
# - CountMatrix.csv
# - coldata.csv

# Main Analysis Steps are: 
# 1) RLE normalization (DESeq2) across all samples element of CountMatrix
# 2) Calc LFC (with respect to control) and LFcutOff (upper 90% quantile of abs(LFC values))
# 3) apeglm shrinking on LFC

# 4) Multiple t-testing with Benjamin-Hochberg correction (padj < 0.05) 
#    and independet hypothesis weighing (IHW) to identify DEGs
#    for LFC and apeglm(LFC) values for H0: LFC = 0; apeglm(LFC) = 0
# 4.1) Annotation of Genelists via org.Dr.eg.db [https://www.bioconductor.org/packages/release/data/annotation/html/org.Dr.eg.db.html]
# 4.2) DEGs with LFC / apeglm(LFC) > < abs(LFcutOff) are kept as potential molecular marker features

# 5) Plotting
# 5.1 MA- & Vulcano plots for N-1 conditions for LFC [LFcut & padj cut]
# 5.2 MA- & Vulcano plots for N-1 conditions for apeglm(LFC) [LFcut & padj cut]
# 5.3 DataQC
#     - Correlation bio. replicates (rlog(meanCounts))
#     - MeanSD plots (rlog(meanCounts))
#     - RLE_normalisation
#     - normCounts_rlogTransformation

# 6) Output: 
# 6.1 DESeq2 Results (xcel; still working on html output though ... )
#     - N-1 DESeq2 result files for LFC
#     - N-1 DESeq2 result files for apeglm(LFC)
# 6.2 DEGs (padj < 0.05 and LFcut for: 
#     - N-1 files for LFC
#     - N-1 files for apeglm(LFC)
############

### LOAD PACKAGES ####
library(DESeq2) #installed
library(IHW)    #installed
library(apeglm) #installed
#library(ashr)  #not needed yet .. failed to install though ... 
library(qvalue)
library(locfdr)
library(ggplot2)
library(ggpubr)
library(hexbin)
####################

### DATA IMPORT ####
# CountMatrix import
count.matrix <- read.csv2(file = "CountMatrix.csv", header = T, row.names = 1, fill = F)

# coldata import
coldata <- read.csv2("coldata.csv", header = TRUE, row.names = 1)
cols <- which(variable.names(coldata) != "X.1" & variable.names(coldata) != "X.2" & variable.names(coldata) != "X.3")
rows <- which(coldata$Condition != "") #this section is needed to avoid messed up data imports from xcel
coldata <- droplevels(coldata[rows,cols]) #clean levels of coldata file! 
rm(rows,cols)

coldata$SampleID <- row.names(coldata) #in case row names get lost down the line

substance <- levels(coldata$Substance) #Extract the name(s) of the tested substance(s)
condition <- levels(coldata$Condition)

# Extract Condition IDs
for (i in condition){
  x <- row.names(coldata[which(coldata$Condition %in% i),])
  assign(paste0("ID.",i),x)
  rm(i,x)
}

# Sort & check correct data import:
count.matrix <- count.matrix[,rownames(coldata)] #match() function could be propably used as well ...
stopifnot(all(rownames(coldata) == colnames(count.matrix))) #Must be TRUE
###################

### CREATE AND NAVIGATE TO OUTPUT FOLDERS ####
# Output folder
out <- "DESeq2_Pairwise"
dir.create(out)
setwd(paste0("./",out))
for (i in c("DEGs","Plots","Results")){
  dir.create(i)
}
for (i in c("DataQC","Heatmaps")){
  dir.create(paste0("Plots/",i))
}
rm(i,out)
############################################

### CREATE DESeq2 OBJECT #####
p <- 0.05 # padj cutoff

message(paste0(" 
 Starting DESeq2 Analysis - Pairwise Wald's t-test...
 padj < ",p,"
 Tested Substance: ",substance,"
               "))

# Import entire count matrix into a DESeq2 Object
dds <- DESeqDataSetFromMatrix(countData = count.matrix,
                              colData = coldata,
                              design = ~ Tank + Condition) #Multifactor Level
# log2 fold change and Wald test p value will be computed for the LAST variable in the design formula!

# Filter Matrix
#keep <- counts(dds)[apply(count.matrix[c(1:9)],1,function(z) !any(z==0)),] #Remove all Genes with a 0 value in their row
keep <- counts(dds)[rowSums(count.matrix) >= length(colnames(count.matrix)), ] 
# length(row.names(keep))  # use this line to check how many genes still remain in your dataset 

# create vector list of removed genes 
x <- row.names(count.matrix)
y <- row.names(keep)
rmv <- setdiff(x,y) # Names of removed genes
rm(x,y)

### NEW, filtered dds object ###
dds <- dds[row.names(keep),]
rm(keep)

# Set factor levels in right order:
# for more details check the DESeq2 Vignette on "Note on factor levels"
if (all(grepl("Exposure",condition[2:4]))==T & length(condition)==4){
  # use this line when dealing with: High, Mid and Low Exposure condition
  print("Sorting 3 Treatments (Low, Mid & High Exposure) and setting Control as reference level.")
  dds$Condition <- factor(dds$Condition, levels = c(condition[1],  # Control
                                                    condition[3],  # LExp
                                                    condition[4],  # MExp
                                                    condition[2])) # HExp
} else if (all(grepl("Exposure",condition[2:3]))==T & length(condition)==3) {
  # use this line when dealing with: High, Mid and Low Exposure condition
  print("Sorting 2 Treatments (Low & High Exposure) and setting Control as reference level.")
  dds$Condition <- factor(dds$Condition, levels = c(condition[1],  # Control
                                                    condition[3],  # LExp
                                                    condition[2])) # HExp
} else {
  # or just specifying the reference level (usually control)
  print( "Set Control to reference level for multiple treatments.")
  #dds$Condition <- factor(dds$Condition, ordered = FALSE )
  dds$Condition <- relevel(dds$Condition, ref = "Control")
}
# Set levels for coldata
coldata$Condition <- factor(coldata$Condition, levels = levels(dds$Condition))

# Get a list of different disp estimates for later QC plotting
fitType = c("parametric", "local", "mean")
estDisp.ls <- list()
for(i in fitType){
  x <- estimateSizeFactors(dds)
  x <- estimateDispersions(x, fitType = i)
  estDisp.ls[[i]] <- x
}
rm(x,i)

# FINALLY RUN DESeq function ----
dds <- DESeq(dds, test = "Wald") # Use only for pairwise comparison
# LRT might be better suited for multiple factor analysis
#dds <- DESeq(dds, test = "LRT", reduced = ~ Tank) #ANOVA-like approach

### DESeq2 counts and Objects ###
# Extract normalized read counts
df.counts <- counts(dds, normalized = T)

# Extract (n+1)log2 transformed mean read counts
ntd <- normTransform(dds)
df.ntd <- assay(ntd)

# Extract rlog transformed mean read counts
rld <- rlog(dds, blind = FALSE)
df.rld <- assay(rld)
###########################

### EXTRACT DESeq2 RESULTS ####
message(" Exporting DESeq normalized count matrix. This might take a while...
 Saving to tab delim txt ...")
# Extract DESEq2 normalized count matrix
write.table(as.data.frame(df.counts), paste0(substance,"_DESeqNormCounts.txt"))
# Extract rlog(norm counts) matrix
write.table(as.data.frame(df.rld), paste0(substance,"_rlogDESeqNormCounts.txt"))

# non shrunk
for (i in resultsNames(dds)[grepl("^Condition_",resultsNames(dds))]){
  x <- gsub("Condition_","", i)
  print(paste0("Extracting DESeq2 results for: ",x," [IHW, non-shrunk LFC]"))
  x <- gsub("_vs_Control","",x)
  x <- gsub("Exposure","",x)
  y <- results(dds, 
               name = i, 
               filterFun = ihw)
  assign(paste0("res.",x),y)
  rm(x,y,i)
}
# Same for reslfs (apeglm shrunk)
for (i in resultsNames(dds)[grepl("^Condition_",resultsNames(dds))]){
  x <- gsub("Condition_","", i)
  print(paste0("Extracting DESeq2 results for: ",x," [IHW, apeglm(LFC) shrunk]"))
  x <- gsub("_vs_Control","",x)
  x <- gsub("Exposure","",x)
  y <- lfcShrink(dds, 
                 coef = i,
                 type = "apeglm",
                 res = get(paste0("res.",x))
                 )
  assign(paste0("reslfs.",x),y)
  rm(x,y,i)
}

# Extract res. variable strings as resNames for loop functions
x <- resultsNames(dds)[grepl("^Condition_",resultsNames(dds))]
x <- gsub("Condition_","", x) # rmv Condition_ string
x <- gsub("_vs_Control","",x) # rmv _vs_Control string
resNames <- gsub("Exposure","",x) #rmv Exposure string
#############################

### Log2FC Cutoff SETTING ####
# via quantile calc:
# qunatile() can use 9 different methods for calculating the quantiles.
# Further details are provided in Hyndman and Fan (1996) who recommended type 8. 
# The default method is type 7, as used by S and by R < 2.0.0.
#
# Type 7:
# m = 1-p. p[k] = (k - 1) / (n - 1). In this case, p[k] = mode[F(x[k])]. This is used by S
#
# Type 8:
# m = (p+1)/3. p[k] = (k - 1/3) / (n + 1/3). Then p[k] =~ median[F(x[k])]. 
# The resulting quantile estimates are approximately median-unbiased regardless of the distribution of x.
##########################
# Only extracted from NON SHRUNK LFC values
for (i in resNames){
  x <- get(paste0("res.",i)) # call any of the res."" objects and assign to x 
  y <- quantile(abs(x$log2FoldChange), na.rm=T,
                type=8,
                probs=seq(0,1,0.05)) #Probability values [0,1] in 0.05 intervals
  z <- as.numeric(y["90%"]) # assign upper 90% quantile value as cut off to object z
  if (z < 1){
    message(paste0("
   |LFcut| for res.",i," [",substance,"] ",round(z,2)))
    assign(paste0("LFcut.",i),z) #assign the value z to new object LFcut."
  } else {
    message(paste0("
     |LFcut| for res.",i," [",substance,"] ",round(z,2),"
     LFcut.",i," > 1; Hence LFC cutoff set to: 1"))
    assign(paste0("LFcut.",i),1) #assign the value 1 to new object LFcut."
  }
  rm(x,y,z,i)
}
#############################

#### Append qvalue to res ####
message("
          Calculating qvalues
        ")
for (k in resNames){
  # non-shrunk
  res <- get(paste0("res.",k))
  res$qval <- qvalue(res$pvalue)$qvalue
  assign(paste0("res.",k),res)
}
for (k in resNames){
  # lfs shrunk
  res <- get(paste0("reslfs.",k))
  res$qval <- qvalue(res$pvalue)$qvalue
  assign(paste0("reslfs.",k),res)
}
rm(k,res)
#############################

### ANNOTATE GENES #####
library(AnnotationDbi)
library(org.Dr.eg.db)
message(" Annotate DESeq2 results ...
        ")
# organism annotation package ("org") for Danio rerio ("Dr") organized as 
# an AnnotationDbi database package ("db"), using Entrez Gene IDs ("eg") as primary key.
# to install: BiocManager::install("org.Dr.eg.db")
# https://www.bioconductor.org/packages/release/data/annotation/html/org.Dr.eg.db.html
# citation("org.Dr.eg.db")
# To get a list of all available key types, use:
# columns(org.Dr.eg.db)

# We can use the mapIds function to add individual columns to our results table
# We provide the row names of our results table as a key, and specify that keytype=ENSEMBL
db.keys.first <- c("ENTREZID",
                   "SYMBOL",
                   "GENENAME",
                   "ZFIN"
                   )
db.keys.unique <- c("IPI",
                    "ENZYME",
                    "ENSEMBLPROT"
                    #"UNIPROT",
                    #"ONTOLOGY",
                    #"GO",
                    #"GOALL"
                    )
# for res
for (k in resNames){
  x <- get(paste0("res.",k))
  for (i in db.keys.first){
    anno <- mapIds(org.Dr.eg.db,
                   keys=row.names(x),
                   keytype="ENSEMBL",
                   column= i,
                   multiVals= "first")
    x[,paste0(i)] <- anno
  }
  for (i in db.keys.unique){
    anno <- mapIds(org.Dr.eg.db,
                   keys=row.names(x),
                   keytype="ENSEMBL",
                   column= i,
                   multiVals= "asNA") #only unique values get ID
    x[,paste0(i)] <- anno
  }
  assign(paste0("res.",k),x)
}
# for reslfs
for (k in resNames){
  x <- get(paste0("reslfs.",k))
  for (i in db.keys.first){
    anno <- mapIds(org.Dr.eg.db,
                   keys=row.names(x),
                   keytype="ENSEMBL",
                   column= i,
                   multiVals= "first")
    x[,paste0(i)] <- anno
  }
  for (i in db.keys.unique){
    anno <- mapIds(org.Dr.eg.db,
                   keys=row.names(x),
                   keytype="ENSEMBL",
                   column= i,
                   multiVals= "asNA") #only unique values get ID
    x[,paste0(i)] <- anno
  }
  assign(paste0("reslfs.",k),x)
}
rm(anno,k,i,x)
detach("package:org.Dr.eg.db", unload = T)
#####################

### EXPORT RESULTS to Xcel #####
message(" Exporting DESeq2 result tables. This might take a while.
 Saving to Xcel...")

## DESEq2 res ##
for (k in resNames){
  res <- get(paste0("res.",k))
  #df <- as.data.frame(matrix(NA,length(rmv),ncol(res)),row.names = rmv)
  #names(df) <- colnames(res)
  #df <- rbind(res,df)
  #stopifnot(nrow(df) == nrow(count.matrix)) # check if dataset is complete!
  write.csv2(as.data.frame(res),
             file = paste0("Results/",substance,"_res_",k,".csv"))
  rm(k,res)
}
## DESEq2 reslfs ##
for (k in resNames){
  res <- get(paste0("reslfs.",k))
  #df <- as.data.frame(matrix(NA,length(rmv),ncol(res)),row.names = rmv)
  #names(df) <- colnames(res)
  #df <- rbind(res,df)
  #stopifnot(nrow(df) == nrow(count.matrix)) # check if dataset is complete!
  write.csv2(as.data.frame(res),
             file = paste0("Results/",substance,"_reslfs_",k,".csv"))
  rm(k,res)
}

message(paste0(" DESeq2 result tables were stored in:
 ",getwd(),"/Results"))

## DEGs ##
message(" Exporting DEGs in Xcel tables 
 Saving ...")

# Create output folders
home <- getwd()
setwd("DEGs/")
dir.create("lfs")
dir.create("lfsFCcut")
dir.create("unshrink")
dir.create("unshrinkFCcut")
setwd(home)

# padj cutoff 
for (k in resNames){
  # non-shrunk
  res <- get(paste0("res.",k)) 
  res1 <- subset(res[order(res$padj),], padj <= p)
  
  res2 <- get(paste0("res.",tail(resNames,1))) #this calls the last object of the resNames vector list; this is supposed to be the highest concentration tested!!!
  res2 <- subset(res2, padj <= p)
  elementHE <- is.element(rownames(res1),rownames(res2))
  
  res1 <- cbind(elementHE, as.data.frame(res1)) # TRUE if element part of DEGs in highest exposure (HE)
  write.csv2(res1,file = paste0("DEGs/unshrink/",substance,"_pcut_",k,".csv"))
  
  # apeglm shrunk
  res <- get(paste0("reslfs.",k)) 
  res1 <- subset(res[order(res$padj),], padj <= p)
  
  res2 <- get(paste0("reslfs.",tail(resNames,1))) #this calls the last object of the resNames vector list; this is supposed to be the highest concentration tested!!!
  res2 <- subset(res2, padj <= p)
  elementHE <- is.element(rownames(res1),rownames(res2))
  
  res1 <- cbind(elementHE, as.data.frame(res1)) # TRUE if element part of DEGs in highest exposure (HE)
  write.csv2(res1,file = paste0("DEGs/lfs/",substance,"_lfs_pcut",k,".csv"))
}

# padj & LFC cutoff
# assign
for (k in resNames){
  # non-shrunk
  res <- get(paste0("res.",k)) 
  LFcut <- get(paste0("LFcut.",k))
  res1 <- subset(res[order(res$padj),],log2FoldChange>LFcut | log2FoldChange<(-LFcut))
  assign(paste0("res.",k,".LFcut"),res1)
  
  # apeglm shrunk
  res <- get(paste0("reslfs.",k))
  LFcut <- get(paste0("LFcut.",k))
  res1 <- subset(res[order(res$padj),],log2FoldChange>LFcut | log2FoldChange<(-LFcut))
  assign(paste0("reslfs.",k,".LFcut"),res1)
}
# export
for(k in resNames){
  # non-shrunk
  res1 <- get(paste0("res.",k,".LFcut"))
  res1 <- subset(res1, padj <= p)
  res2 <- get(paste0("res.",tail(resNames,1),".LFcut")) #this calls the last object of the resNames vector list; this is supposed to be the highest concentration tested!!!
  res2 <- subset(res2, padj <= p)
  elementHE <- is.element(rownames(res1),rownames(res2))
  
  res1 <- cbind(elementHE, as.data.frame(res1)) # TRUE if element part of DEGs in highest exposure (HE)
  write.csv2(res1,file = paste0("DEGs/unshrinkFCcut/",substance,"_pcut_LFcut",k,".csv"))
  
  # apeglm shrunk
  res1 <- get(paste0("reslfs.",k,".LFcut"))
  res1 <- subset(res1, padj <= p)
  res2 <- get(paste0("reslfs.",tail(resNames,1),".LFcut")) #this calls the last object of the resNames vector list; this is supposed to be the highest concentration tested!!!
  res2 <- subset(res2, padj <= p)
  elementHE <- is.element(rownames(res1),rownames(res2))
  
  res1 <- cbind(elementHE, as.data.frame(res1)) # TRUE if element part of DEGs in highest exposure (HE)
  write.csv2(res1,file = paste0("DEGs/lfsFCcut/",substance,"_lfs_pcut_LFcut",k,".csv"))
}
rm(res,res1,res2,LFcut,k,elementHE)
message(paste0(" 
Finished DESeq2 Statistical testing.
Start with plotting the DESeq2 results.
 "))
##############################

############################################ PLOTTING ##############################################
message(" 
 Start Data Composition and QC plotting ...")
### QC - dispersion estimates for DESeq2's normalization model ----
pdf(file = paste0("Plots/DataQC/",substance,"_DispEstimates.pdf"), #print directly to pdf
    width = 7, height = 5.9,
    onefile = T,        #multiple figures in one file
    bg = "transparent" #Background color
)
for(i in fitType){
  print(
    plotDispEsts(estDisp.ls[[i]], main = paste(substance,'- DESeq2 fit type:',i))
    )
}
dev.off()
rm(estDisp.ls)

##################################################################

### QC - pvalue and LFC distribution #####
#pdf(file = paste0("Plots/DataQC/",substance,"_pvalue_LFC_distr.pdf"), #print directly to pdf
#    width = 4*length(resNames), height = 12, onefile = T, bg = "transparent")
## res -------
png(filename = paste0("Plots/DataQC/",substance,"_pvalue_LFC_distr.png"),
    width = 6*length(resNames), height = 18, units = "cm", bg = "white",
    pointsize = 7, res = 450)
par(mfrow=c(3,length(resNames)))
# loop for pvalue distr 
for (k in resNames){
  res <- get(paste0("res.",k))
  hist(res@listData[["pvalue"]],
       main= paste0(substance," - res.",k),
       col = "gray50",
       border = "gray50",
       ylab = "Frequency",
       xlab = "pvalues",
       breaks = 500)
}
# loop for pvalue vs padj (BH) & qval 
q <- 0.05 #specify Q value cut off
pval <- 0.05 # specify pval cut off 
for (k in resNames){
  res <- get(paste0("res.",k))
  plot(res@listData[["pvalue"]],res@listData[["padj"]],
       xlab="pvalues",
       ylab="conversion of pvalues",
       xlim= c(0,0.25),
       main= paste0(substance," - res.",k),
       col = "dodgerblue1",
       pch = 20
  )
  points(res@listData[["pvalue"]],res@listData[["qval"]], 
         col="springgreen3", 
         pch = 20)
  xpos <- 0.1
  abline(h=0.05, lwd=1, lty=2)
  abline(v=0.05, lwd=1, lty=2)
  legend(x=xpos,y=0.225,
         text.col = "black",
         bty = "n",
         legend = paste("Genes with pvalue <",pval,":",sum(res$pvalue <= pval, na.rm = T))
  )
  legend(x=xpos,y=0.135,
         text.col = "springgreen3",
         bty = "n",
         legend = paste("Genes with Qvalue <",q,":",sum(res$qval <= q, na.rm = T))
  )
  legend(x=xpos,y=0.045,
         text.col = "dodgerblue",
         bty = "n",
         legend = paste("Genes with padj (BH) <",p,":",sum(res$padj <= p, na.rm = T))
  )
}
rm(k,res)
  # loop for LFC distribution
for (k in resNames){
  res <- get(paste0("res.",k))
  LFcut <- get(paste0("LFcut.",k))
  hist(res$log2FoldChange,
       breaks = 1000,
       main=paste0(substance," - res.",k),
       col = "gray50",
       border = "gray50",
       xlab = "Log2 FC",
       #xlim = c(-7*LFcut,7*LFcut),
       #ylim = c(0,3000),
       xlim = c(-2,2)
  )
  abline(v=c(LFcut,-LFcut), 
         col="dodgerblue",
         lty = 2,
         lwd= 1 )
  legend(x="topright",
         text.col = "dodgerblue",
         bty = "n",
         legend = paste("LFcut =",round(LFcut, digits = 2)))
  rm(res,k)
}
dev.off()

## reslfs --------
png(filename = paste0("Plots/DataQC/",substance,"_pvalue_lfsLFC_distr.png"),
    width = 6*length(resNames), height = 18, units = "cm", bg = "white",
    pointsize = 7, res = 450)
par(mfrow=c(3,length(resNames)))
for (k in resNames){
  res <- get(paste0("reslfs.",k))
  hist(res@listData[["pvalue"]],
       main= paste0(substance," - reslfs.",k),
       col = "gray50",
       border = "gray50",
       ylab = "Frequency",
       xlab = "pvalues",
       breaks = 500)
}
# loop for pvalue vs padj (BH) & qval
q <- 0.05 #specify Q value cut off
pval <- 0.05 # specify pval cut off 
for (k in resNames){
  res <- get(paste0("reslfs.",k))
  plot(res@listData[["pvalue"]],res@listData[["padj"]],
       xlab="pvalues",
       ylab="conversion of pvalues",
       xlim= c(0,0.25),
       main= paste0(substance," - reslfs.",k),
       col = "dodgerblue1",
       pch = 20
  )
  points(res@listData[["pvalue"]],res@listData[["qval"]], 
         col="springgreen3", 
         pch = 20)
  xpos <- 0.1
  abline(h=0.05, lwd=1, lty=2)
  abline(v=0.05, lwd=1, lty=2)
  legend(x=xpos,y=0.225,
         text.col = "black",
         bty = "n",
         legend = paste("Genes with pvalue <",pval,":",sum(res$pvalue <= pval, na.rm = T))
  )
  legend(x=xpos,y=0.135,
         text.col = "springgreen3",
         bty = "n",
         legend = paste("Genes with Qvalue <",q,":",sum(res$qval <= q, na.rm = T))
  )
  legend(x=xpos,y=0.045,
         text.col = "dodgerblue",
         bty = "n",
         legend = paste("Genes with padj (BH) <",p,":",sum(res$padj <= p, na.rm = T))
  )
}
rm(k,res)
# loop for LFC distribution
for (k in resNames){
  res <- get(paste0("reslfs.",k))
  LFcut <- get(paste0("LFcut.",k))
  hist(res$log2FoldChange,
       breaks = 1000,
       main=paste0(substance," - reslfs.",k),
       col = "gray50",
       border = "gray50",
       xlab = "Log2 FC",
       #xlim = c(-7*LFcut,7*LFcut),
       #ylim = c(0,3000),
       xlim = c(-2,2)
  )
  abline(v=c(LFcut,-LFcut), 
         col="dodgerblue",
         lty = 2,
         lwd= 1 )
  legend(x="topright",
         text.col = "dodgerblue",
         bty = "n",
         legend = paste("LFcut =",round(LFcut, digits = 2)))
  rm(res,k)
}
dev.off()
#######################################

### QC - RLE COUNT NORMALISATION ####
# plotting raw counts and mean counts next to each other 
pdf(file = paste0("Plots/DataQC/",substance,"_RLE_normalisation.pdf"), #print directly to pdf
    width = 1.3*nrow(coldata),       
    height = 11,
    onefile = T,        #multiple figures in one file
    title = "",
    #paper = "a4",         #a4r = landscape DinA4 (r=rotated)
    bg = "transparent", #Background color
    fg ="black"         #Foreground color
)
par(mfrow=c(2,2))

barplot(colSums(counts(dds)),
        #xlab = "Samples",
        ylab = "Total gene counts",
        main = "Raw counts",
        las = 3, #rotating sample labels 90°
        ylim = c(0,1.2*max(colSums(counts(dds))))
)
barplot(colSums(counts(dds, normalized=T)),  #plot mean normalized counts
        #xlab = "Samples",
        ylab = "Total gene counts",
        main = "RLE normalized counts",
        las = 3, #rotating sample labels 90°
        ylim = c(0,1.2*max(colSums(counts(dds))))
)
boxplot(log10(counts(dds)+1),
        #xlab="Samples",
        ylab="Log10(Gene counts + 1)",
        las = 3, #rotating sample labels 90°
        main="Raw counts")
boxplot(log10(counts(dds, normalized=T)+1),
        #xlab="Samples",
        ylab="Log10(Gene counts + 1)",
        las = 3, #rotating sample labels 90°
        main="RLE normalized counts")
dev.off()
####################################

### QC - RLE NORMALIZED COUNT TRANSFORMATION ####
pdf(file = paste0("Plots/DataQC/",substance,"_NormCount_Transformation.pdf"), #print directly to pdf
    width = 5*length(resNames),       
    height = 4,
    onefile = T,        #multiple figures in one file
    title = "RLE normalized Count's transformation",
    #paper = "a4",         #a4r = landscape DinA4 (r=rotated)
    bg = "transparent", #Background color
    fg ="black"         #Foreground color
)
# plot...
par(mfrow=c(1,3))
boxplot(df.counts, notch = TRUE,
        las = 3, #rotating sample labels 90°
        main = "Normalized read counts - No zero counts",
        ylab = "read counts")

boxplot(df.ntd, notch = TRUE,
        las = 3, #rotating sample labels 90°
        main = "log2 Transformation",
        ylab = "log2 (norm. read counts + 1)")

boxplot(df.rld, notch = TRUE,
        las = 3, #rotating sample labels 90°
        main = "rlog Transformation",
        ylab = "rlog (norm. read counts)")

dev.off()
while (!is.null(dev.list()))  dev.off()
################################################

### QC - PearsonCor biol. replicates ####
skip <- T
if(skip == F){
  # This section of code assumes 3 biol. replicates! If this differs please change manually!
  pdf(file = paste0("Plots/DataQC/",substance,"_Correlation.pdf"), 
      width = 3.33333*length(ID.Control),       
      height = 3.6*length(condition),
      onefile = T,        #multiple figures in one file
      title = "",
      #paper = "a4",         #a4r = landscape DinA4 (r=rotated)
      bg = "transparent", #Background color
      fg ="black"         #Foreground color
  )
  par(mfrow=c(length(condition),length(ID.Control)))
  
  list <- c("df.counts","df.ntd","df.rld")
  for (k in list){
    df <- get(k)
    for (i in condition){
      ID <- get(paste0("ID.",i)) # Get object containing Sample IDs for a specific condition
      #1
      plot(df[,c(ID[1],ID[2])],
           cex =.1,
           main = paste("norm. reads -",i,"[",coldata[ID[1],"Tank"],"vs",coldata[ID[2],"Tank"],"]")
      )
      x <- df[,ID[1]]
      y <- df[,ID[2]]
      legend(x="bottomright", legend = paste("Pearson Cor =",round(cor(x,y),3)))
      #2
      plot(df[,c(ID[1],ID[3])],
           cex =.1,
           main = paste("norm. reads -",i,"[",coldata[ID[1],"Tank"],"vs",coldata[ID[3],"Tank"],"]")
      )
      x <- df[,ID[1]]
      y <- df[,ID[3]]
      legend(x="bottomright", legend = paste("Pearson Cor =",round(cor(x,y),3)))
      #3
      plot(df[,c(ID[2],ID[3])],
           cex =.1,
           main = paste("norm. reads -",i,"[",coldata[ID[2],"Tank"],"vs",coldata[ID[3],"Tank"],"]")
      )
      x <- df[,ID[2]]
      y <- df[,ID[3]]
      legend(x="bottomright", legend = paste("Pearson Cor =",round(cor(x,y),3)))
      rm(i,ID,x,y)
    }
  }
  rm(k,list,df)
  dev.off()
  while (!is.null(dev.list()))  dev.off()
}
########################################

### QC - PearsonCor biol. replicates - SHINY ####
## plot function ##
corplot.1 <- function(df1) ggplot(df1) +
  geom_point(aes(x, y, col = d), size = .8, alpha = .4) +
  labs(x=paste0(ID[p1],transf),y=paste0(ID[p2],transf),
       title = (paste0(substance,": ",i)),
       subtitle = (paste("[",coldata[ID[p1],"Tank"],"vs",coldata[ID[p2],"Tank"],"]"))) +
  annotate("text",
           label = paste0("Pearson = ",round(cor(df1$x,df1$y),2)),
           x = (max(df1$x)), 
           y = (min(df1$y)),
           hjust=1, vjust=0) + 
  annotate("text",
           label = paste0("R2 = ",round((cor(df1$x,df1$y))^2,2)),
           x = (min(df1$x)), 
           y = (max(df1$y)),
           hjust=0, vjust=1) + 
  scale_color_identity() +
  coord_equal(ratio=1) + 
  theme_bw()

# This section of code assumes 3 biol. replicates! If this differs please change manually!
#pdf(file = paste0("Plots/DataQC/",substance,"_Correlation_shiny.pdf"), 
#    width = 4.33333*length(ID.Control), height = 4.6*length(condition), onefile = T, bg = "transparent")

# Log10 mean counts -------------------------------------------------------------------
png(filename = paste0("Plots/DataQC/",substance,"_Correlation_log10.png"),
    width = 7.33333*length(ID.Control), height = 7.6*length(condition), units = "cm", bg = "white",
    pointsize = 1, res = 450)
df <- df.ntd #INPUT 
transf <- " - [log2(reads+1)]"
gg.list <- list()
for (i in condition){
  ID <- get(paste0("ID.",i)) # Get object containing Sample IDs for a specific condition
  p1 <- 1 #ID vector position
  p2 <- 2 #ID vector position
  df1 <- data.frame(x = df[,ID[p1]],
                    y = df[,ID[p2]],
                    d = densCols(df[,ID[p1]], df[,ID[p2]], colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
  )
  gg.list[[paste0(i,".",ID[p1],"vs",ID[p2])]] <- corplot.1(df1)
  
  p1 <- 1 #ID vector position
  p2 <- 3 #ID vector position
  df1 <- data.frame(x = df[,ID[p1]],
                    y = df[,ID[p2]],
                    d = densCols(df[,ID[p1]], df[,ID[p2]], colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
  )
  gg.list[[paste0(i,".",ID[p1],"vs",ID[p2])]] <- corplot.1(df1)
  
  p1 <- 2 #ID vector position
  p2 <- 3 #ID vector position
  df1 <- data.frame(x = df[,ID[p1]],
                    y = df[,ID[p2]],
                    d = densCols(df[,ID[p1]], df[,ID[p2]], colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
  )
  gg.list[[paste0(i,".",ID[p1],"vs",ID[p2])]] <- corplot.1(df1)
}
ggobj <- ls(gg.list)
print(ggpubr::ggarrange(plotlist = gg.list, ncol = 3, nrow = (length(ggobj)/3)))
dev.off()

# rlog mean counts --------------------------------------------------------------------
png(filename = paste0("Plots/DataQC/",substance,"_Correlation_rlog.png"),
    width = 7.33333*length(ID.Control), height = 7.6*length(condition), units = "cm", bg = "white",
    pointsize = 1, res = 450)
df <- df.rld #INPUT
transf <- " - [rlog(reads)]"
gg.list <- list()
for (i in condition){
  ID <- get(paste0("ID.",i)) # Get object containing Sample IDs for a specific condition
  p1 <- 1 #ID vector position
  p2 <- 2 #ID vector position
  df1 <- data.frame(x = df[,ID[p1]],
                    y = df[,ID[p2]],
                    d = densCols(df[,ID[p1]], df[,ID[p2]], colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
  )
  gg.list[[paste0(i,".",ID[p1],"vs",ID[p2])]] <- corplot.1(df1)
  
  p1 <- 1 #ID vector position
  p2 <- 3 #ID vector position
  df1 <- data.frame(x = df[,ID[p1]],
                    y = df[,ID[p2]],
                    d = densCols(df[,ID[p1]], df[,ID[p2]], colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
  )
  gg.list[[paste0(i,".",ID[p1],"vs",ID[p2])]] <- corplot.1(df1)
  
  p1 <- 2 #ID vector position
  p2 <- 3 #ID vector position
  df1 <- data.frame(x = df[,ID[p1]],
                    y = df[,ID[p2]],
                    d = densCols(df[,ID[p1]], df[,ID[p2]], colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
  )
  gg.list[[paste0(i,".",ID[p1],"vs",ID[p2])]] <- corplot.1(df1)
}
ggobj <- ls(gg.list)
print(ggpubr::ggarrange(plotlist = gg.list, ncol = 3, nrow = (length(ggobj)/3)))
dev.off()
################################################

### QC - meanSdPlot DATA VARIANCE ####
# BiocManager::install("vsn")
library(vsn)
# plotting SD of transformed data across samples against the mean count
sdp1 <- meanSdPlot(df.counts, ranks = T, plot = F)
sdp1 <- sdp1$gg + ggtitle("Mean counts - Ranked")+scale_y_continuous(trans='log10')+ylab("sd (log10 scale)") +
  theme_bw()

sdp2 <- meanSdPlot(df.ntd, ranks = T, plot = F)
sdp2 <- sdp2$gg + ggtitle("log2 (mean counts) - Ranked") +
  theme_bw() #+ ylim(0,3)

sdp3 <- meanSdPlot(df.rld, ranks = T, plot = F)
sdp3 <- sdp3$gg + ggtitle("rlog (mean counts) - Ranked") +
  theme_bw() #+ ylim(0,3)

sdp1b <- meanSdPlot(df.counts, ranks = F, plot = F) #original scale with ranks=F
sdp1b <- sdp1b$gg + ggtitle("Mean counts") +
  theme_bw()

sdp2b <- meanSdPlot(df.ntd, ranks = F, plot = F) #original scale with ranks=F
sdp2b <- sdp2b$gg + ggtitle("log2 (mean counts)") +
  theme_bw() #+ ylim(0,3)

sdp3b <- meanSdPlot(df.rld, ranks = F, plot = F) #original scale with ranks=F
sdp3b <- sdp3b$gg + ggtitle("rlog (mean counts)") + 
  theme_bw()#+ ylim(0,3)

pdf(file = paste0("Plots/DataQC/",substance,"_meanSD.pdf"), #print directly to pdf
    width = 17,       
    height = 9.6,
    onefile = T,           #multiple figures in one file
    title = "",
    #paper = "a4r",         #a4r = landscape DinA4 (r=rotated)
    bg = "transparent", #Background color
    fg ="black"         #Foreground color
)
print(ggarrange(sdp1, sdp2, sdp3, sdp1b, sdp2b, sdp3b,
                labels = c("A1", "B1", "C1","A2", "B2", "C2"),
                ncol = 3, nrow =2))
dev.off()
# Clean out env
rm(sdp1,sdp2,sdp3,sdp1b,sdp2b,sdp3b)
detach("package:vsn", unload = TRUE)


#####################################

### MA plots & Vulcano plots ####
MAfun <- function(res,title,topgenes, Symbol, shrink, LFcut){
  if(missing(title)){title = ''}
  if(missing(shrink)){shrink = F}
  if(missing(Symbol)){Symbol = F}
  if(missing(topgenes)){topgenes = 10}
  ggpubr::ggmaplot(
    res, 
    fdr = p, fc = 2^(LFcut), size = .6,
    top = topgenes,
    select.top.method = 'fc', # fc or padj
    palette = c("#B31B21", "#1465AC", "darkgray"),
    main = paste(substance,title,"[ LFcut:",round(LFcut,2),"]"),
    legend = "bottom",
    genenames = if(Symbol == F){NULL}else{as.vector(res$SYMBOL)},
    ylab = if(shrink == F){bquote(~Log[2]~ "fold change")}else{bquote("apeglm ("~Log[2]~"fc)")},
    xlab = bquote(~Log[2]~ "mean expression"),
    font.label = c("bold", 11), label.rectangle = F,
    font.legend = "bold",
    font.main = "bold",
    ggtheme = ggplot2::theme_light()
  ) + theme(aspect.ratio = 1, plot.title = element_text(size = 10, face = 'plain'), #line = element_line(size = .1),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            text = element_text(size = 1.5, face = 'plain')) +
    geom_text(aes(x = max(log2(res$baseMean))-.6, y = LFcut+.3, label = round(LFcut,2)), size = 3) + 
    geom_text(aes(x = max(log2(res$baseMean))-.65, y = -LFcut-.3, label = round(-LFcut,2)), size = 3)
  #+ geom_hline(yintercept = c(-LFcut,LFcut)) # to display lf cut line
}
VulcFun <- function(res,title,topgenes, Symbol, shrink, LFcut){
  if(missing(shrink)){shrink = F}
  if(missing(title)){title = ''}
  if(missing(topgenes)){topgenes = 10}
  if(missing(Symbol)){Symbol = F}
  
  if(Symbol == T){
    select <- res[order(res$padj),"SYMBOL"]}else{
      select <- rownames(res[order(res$padj),])}
  
  EnhancedVolcano::EnhancedVolcano(
    res, x = "log2FoldChange", y = "padj",
    title = paste(substance,title,"[ LFcut:",round(LFcut,2),"]"),
    #subtitle = if(shrink == F){"IHW; Non-shrunk"}else{"IHW; apeglm"},
    subtitle = NULL,
    lab = if(Symbol == F){rownames(res)}else{res$SYMBOL},
    selectLab = if(topgenes == 0){select[NULL]}else{select[1:topgenes]}, # select topgenes based on padj value
    legend = c('NS','Log2 FC','padj','padj & Log2 FC'),
    xlab =  if(shrink == F){bquote(~Log[2]~ "fold change")}else{bquote("apeglm ("~Log[2]~"fc)")},
    ylab = bquote(~-Log[10]~italic(padj)),
    FCcutoff = LFcut,  #default 1
    pCutoff = p,       #default 10e-6
    labSize = 3.0,
    pointSize = .85, #default 0.8
    col = c("grey30", "grey30", "royalblue", "red2"),
    shape = c(1, 0, 17, 19),   #default 19, http://sape.inf.usi.ch/quick-reference/ggplot2/shape for details
    colAlpha = 0.4, #default 0.5
    hline = c(0.01, 0.001),
    hlineCol = c('grey40','grey55'),
    hlineType = 'dotted',
    hlineWidth = 0.6,
    gridlines.major = T,
    gridlines.minor = T,
    legendVisible = F,
    drawConnectors = TRUE,
    #widthConnectors = 0.2,
    #colConnectors = 'grey15'
  ) + theme_light() +
    theme(aspect.ratio = 1, plot.title = element_text(size = 10, face = 'plain'), line = element_line(size = .6),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10))
}
# NON SHRUNK -----
gg.list <- list() # Create gg.list to store plotting objects
# ENSEMBL IDs
for (k in resNames){
  res <- get(paste0("res.",k))
  LFcut <- get(paste0("LFcut.",k))
  # MA plot
  gg <- MAfun(res, topgenes = 0, LFcut = LFcut,  title = paste("-",k))
  n <- paste0("MA.",k)
  gg.list[[n]] <- gg
  # Vulcano
  gg <- VulcFun(res, topgenes = 0, LFcut = LFcut, title = paste("-",k))
  n <- paste0("Vulc.",k)
  gg.list[[n]] <- gg
  #rm(gg,n,res,LFcut)
}
# SYMBOL
for (k in resNames){
  res <- get(paste0("res.",k))
  LFcut <- get(paste0("LFcut.",k))
  # MA plot
  gg <- MAfun(res, topgenes = 10,LFcut = LFcut, title = paste("-",k), Symbol = T)
  n <- paste0("symMA.",k)
  gg.list[[n]] <- gg
  # Vulcano
  gg <- VulcFun(res, topgenes = 10,LFcut = LFcut, title = paste("-",k), Symbol = T)
  n <- paste0("symVulc.",k)
  gg.list[[n]] <- gg
  #rm(gg,n,res,LFcut)
}
# plotting: 
message(" 
 Start Shiny MA & Vulcano plotting ...")

pdf(file = paste0("Plots/",substance,"_MAplots_Vulcano.pdf"), #print directly to pdf
    width = 4*length(resNames), height = 8, onefile = T, bg = "transparent")
x <- length(names(gg.list))/2
plotseq <- c(grep('MA.',names(gg.list)[1:x]),
             grep('Vulc.',names(gg.list)[1:x]))
print(ggpubr::ggarrange(plotlist = gg.list[plotseq], ncol = length(resNames), nrow =2))
print(ggpubr::ggarrange(plotlist = gg.list[plotseq+x], ncol = length(resNames), nrow =2))
dev.off()
while (!is.null(dev.list()))  dev.off()

# SHRUNK -----
# Create gg.list to store objects for plotting in there
gg.list <- list()
# ENSEMBL IDs
for (k in resNames){
  res <- get(paste0("reslfs.",k))
  LFcut <- get(paste0("LFcut.",k))
  # MA plot
  gg <- MAfun(res, topgenes = 0, LFcut = LFcut, title = paste("-",k), shrink = T)
  n <- paste0("MA.",k)
  gg.list[[n]] <- gg
  # Vulcano
  gg <- VulcFun(res, topgenes = 0, LFcut = LFcut, title = paste("-",k), shrink = T)
  n <- paste0("Vulc.",k)
  gg.list[[n]] <- gg
  rm(gg,n,res,LFcut)
}
# SYMBOL
for (k in resNames){
  res <- get(paste0("reslfs.",k))
  LFcut <- get(paste0("LFcut.",k))
  # MA plot
  gg <- MAfun(res, topgenes = 10, LFcut = LFcut, title = paste("-",k), Symbol = T, shrink = T)
  n <- paste0("symMA.",k)
  gg.list[[n]] <- gg
  # Vulcano
  gg <- VulcFun(res, topgenes = 10, LFcut = LFcut, title = paste("-",k), Symbol = T, shrink = T)
  n <- paste0("symVulc.",k)
  gg.list[[n]] <- gg
  #rm(gg,n,res,LFcut)
}
# plotting: 
message(" 
 Start Shiny MA & Vulcano plotting ...")

pdf(file = paste0("Plots/",substance,"_MAplots_Vulcano_lfs.pdf"), #print directly to pdf
    width = 4*length(resNames), height = 8, onefile = T, bg = "transparent")
x <- length(names(gg.list))/2
plotseq <- c(grep('MA.',names(gg.list)[1:x]),
             grep('Vulc.',names(gg.list)[1:x]))
print(ggpubr::ggarrange(plotlist = gg.list[plotseq], ncol = length(resNames), nrow =2))
print(ggpubr::ggarrange(plotlist = gg.list[plotseq+x], ncol = length(resNames), nrow =2))
dev.off()
while (!is.null(dev.list()))  dev.off()
################################

### Venn Diagram ####
# Display overlapping gene clusters for DEGs across all treatments
# 1) padj 
# 2) padj LFcut 
# 3) lfs padj 
# 4) lfs padj LFcut
message(" 
 Start Shiny Venn Diagram plotting ...")
pdf(file = paste0("Plots/",substance,"_Venn.pdf"),
    #paper = "a4",
    onefile = T,
    bg = "transparent", #Background color
    fg ="black",        #Foreground color
    width = 9.7,
    height = 8)

# make color palette 
#colors <- gplots::greenred(length(resNames))
#colors <- heat.colors(length(resNames))
#colors <- rainbow(length(resNames))
colors <- topo.colors(length(resNames))

# create empty list to store objects in there
venn.ls <- list()
for(k in resNames){
  res <- get(paste0("res.",k))
  x <- rownames(subset(res, padj < p))
  venn.ls[[k]] <- x
  rm(res,x)
}
gg <- plot(eulerr::euler(venn.ls),
           fills = list(fill = colors, alpha = 0.4),
           labels = list(col = "black", font = 4),
           legend = list(col = "black", font = 4),
           main = paste0(" [padj <",p,"]"),
           quantities = TRUE,
           shape = "ellipse",
           lty = 0)
print(gg)

#LFcut
venn.ls <- list()
for(k in resNames){
  res <- get(paste0("res.",k,".LFcut"))
  x <- rownames(subset(res, padj < p))
  n <- paste0(k,".LFcut")
  venn.ls[[n]] <- x
  rm(res,x,n)
}
gg <- plot(eulerr::euler(venn.ls),
           fills = list(fill = colors, alpha = 0.4),
           labels = list(col = "black", font = 4),
           legend = list(col = "black", font = 4),
           main = paste0(" [padj <",p,"]"),
           quantities = TRUE,
           shape = "ellipse",
           lty = 0)
print(gg)

# lfs
venn.ls <- list()
for(k in resNames){
  res <- get(paste0("reslfs.",k))
  x <- rownames(subset(res, padj < p))
  n <- paste0(k,".lfs")
  venn.ls[[n]] <- x
  rm(res,x,n)
}
gg <- plot(eulerr::euler(venn.ls),
           fills = list(fill = colors, alpha = 0.4),
           labels = list(col = "black", font = 4),
           legend = list(col = "black", font = 4),
           main = paste0(" [padj <",p,"]"),
           quantities = TRUE,
           shape = "ellipse",
           lty = 0)
print(gg)

#lfs & LFcut
venn.ls <- list()
for(k in resNames){
  res <- get(paste0("reslfs.",k,".LFcut"))
  x <- rownames(subset(res, padj < p))
  n <- paste0(k,".lfsLFcut")
  venn.ls[[n]] <- x
  rm(res,x,n)
}
gg <- plot(eulerr::euler(venn.ls),
           fills = list(fill = colors, alpha = 0.4),
           labels = list(col = "black", font = 4),
           legend = list(col = "black", font = 4),
           main = paste0(" [padj <",p,"]"),
           quantities = TRUE,
           shape = "ellipse",
           lty = 0)
print(gg)
dev.off()
rm(venn.ls,gg,colors)
####################

### HCLUST ####
message(" 
 Start Hierarchical Clustering ...")
## Pearson
pdf(paste0("Plots/",substance,"_hclust.pdf"),
    paper = "a4r",
    onefile = T,
    bg = "transparent", #Background color
    fg ="black",        #Foreground color
    width = 11.6,
    height = 8.5)
par(mfrow=c(2,3))
for (k in c("average","complete","mcquitty")){
  for(i in c("euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard")){
    dis <- vegan::vegdist(1-cor(df.rld, method = "pearson"),
                          method = i)
    plot(hclust(dis, method = k), # clustering: "single", "complete", "average", "mcquitty", "median" or "centroid"
         labels = paste(rld$Condition, rld$Tank, sep="-"),
         ylab = "Pearson correlation",
         main = paste(i))
  }
}
## Spearman
for (k in c("average","complete","mcquitty")){
  for(i in c("euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard")){
    dis <- vegan::vegdist(1-cor(df.rld, method = "spearman"),
                          method = i)
    plot(hclust(dis, method = k), # clustering: "single", "complete", "average", "mcquitty", "median" or "centroid"
         labels = paste(rld$Condition, rld$Tank, sep="-"),
         ylab = "Spearman correlation",
         main = paste(i))
  }
}
dev.off()
##############

### Sample Distance Mtx ####
message(" 
 Start Sample Distance Mtx Heatmap  ...")
pdf(file = paste0("Plots/",substance,"_SampleDistMtx.pdf"), #print directly to pdf
    width = ncol(df.rld), height = (2/3)*ncol(df.rld), onefile = T, bg = "transparent",
)
for(method in c('euclidean','maximum','manhattan')){
  # transpose input, calculate sample distance and create a distance matrix
  disMethod <- method
  sampleDist <- dist(t(df.rld), method = disMethod) #must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
  distMtx <- as.matrix(sampleDist)
  rownames(distMtx) <- paste(rld$Condition, rld$Tank, sep="-")
  # create annotation object for heatmap
  ann_col <- subset(coldata, select = 'Condition')
  ann_row <- subset(coldata, select = 'Tank')
  rownames(ann_row) <- paste(rld$Condition, rld$Tank, sep="-")
  # set colors
  Tanks <- levels(coldata$Tank)
  Conditions <- levels(coldata$Condition)
  col.Cond <- colorRampPalette(c('green',"gold1","indianred3"))(length(Conditions))
  col.Tank <- colorRampPalette(c("gray95","gray50"))(length(Tanks))
  ann_colors <- list(Condition = setNames(col.Cond, Conditions),
                     Tank = setNames(col.Tank, Tanks))
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(65)
  # plot heatmap
  print(
    pheatmap::pheatmap(
      distMtx,
      angle_col = "90",
      treeheight_col = 40, #default 50
      fontsize = 12, #default 10
      drop_levels = T,
      display_numbers = T,
      number_format = if(method == 'manhattan'){'%1.0f'}else{'%.1f'},
      fontsize_number = 0.7*12,
      main = paste0(substance,' - SampleDist [',disMethod,']'),
      clustering_distance_rows = sampleDist,
      clustering_distance_cols = sampleDist,
      cellwidth = 26,
      cellheight = 26,
      annotation_col = ann_col, # Condition!
      annotation_row = ann_row, # Tank!
      annotation_colors = ann_colors,
      col = colors
    )
  )
}
rm(sampleDist,distMtx,disMethod,ann_col,ann_row,ann_colors,colors)
dev.off()
while(!is.null(dev.list())) dev.off()
###########################

### PCA & t-SNE ####
message(" 
 Start t-SNE & PCA plotting ...")
# rlog
X <- df.rld 
gg.list <- list()
for (top in c(100,500,1000,nrow(df.rld))){
  ### PCA -----------
  col <- colorRampPalette(c('green',"gold1","indianred3"))(length(Conditions))
  pca <- plotPCA(rld, intgroup = c("Condition","Tank"), returnData=T, ntop=top) 
  percentVar <- round(100 * attr(pca, "percentVar"),1)
  pca <- ggplot(pca, aes(PC1, PC2, color=Condition, shape=Tank),size=3) +
    geom_point(size=3, alpha = .7) +
    ggtitle(paste0(substance," rlog(counts) ntop: ",top," ; ",sum(percentVar),"% Var")) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    scale_colour_manual(values = col) +
    #coord_fixed() +
    theme_bw() +
    theme(aspect.ratio = 1)
  n <- paste0("PCA_",top)
  gg.list[[n]] <- pca
  
  ### t-SNE -----------
  # Subsetting the top genes with max variance
  Xvar <- apply(df.rld,1,var)
  Xvar <- sort(Xvar, decreasing = T)
  select <- names(Xvar[1:top])
  X <- df.rld[select,]
  Xt <- t(X)
  # perform Rtsne calc.
  perplex <- (nrow(Xt)-1)/3
  set.seed(42) #to make results reproducable! 
  tsne <- Rtsne::Rtsne(Xt, 
                dims = 2, 
                initial_dims = ncol(X),
                #initial_dims = 50,
                perplexity = perplex, #(should not be bigger than 3 * perplexity < nrow(X) - 1, see details for interpretation)
                theta = 0.0, #Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE (default: 0.5)
                check_duplicates = F,
                pca = TRUE, 
                partial_pca = F, #(requires the irlba package). This is faster for large input matrices (default: FALSE)
                max_iter = 10000, #number of iterations
                is_distance = FALSE,
                pca_center = TRUE, 
                pca_scale = FALSE,
                normalize = F, #Default True; Set to F as RLE norm was performed prior! 
                verbose = T,
                num_threads = 2)
  # ggplot
  ggtsne <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2],
                       row.names = rownames(coldata),
                       Condition = coldata$Condition,
                       Tank = coldata$Tank,
                       ExtrDate = coldata$SamplingDate)
  gg <- ggplot(ggtsne)+geom_point(aes(x=x,y=y,color = Condition,shape = Tank),size=3, alpha = .7)+
    ggtitle(paste0(substance,": t-SNE on rlog[mean counts]; ntop: ",top)) +
    xlab('t-SNE 1') + ylab('t-SNE 2') +
    scale_colour_manual(values = col) +
    theme_bw() +
    theme(aspect.ratio = 1)
  n <- paste0("tSNE_",top)
  gg.list[[n]] <- gg
}
rm(Xvar,Xt,X,perplex,top,select,gg,ggtsne,tsne,pca,percentVar,n)
# Plotting
pdf(file = paste0("Plots/",substance,"_PCA_tSNE.pdf"),
    width = 22,       
    height = 10,
    onefile = T,           #multiple figures in one file
    #title = "",
    #paper = "a4r",         #a4r = landscape DinA4 (r=rotated)
    bg = "transparent", #Background color
    fg ="black"         #Foreground color
)
print(ggpubr::ggarrange(plotlist = gg.list[c(1,3,5,7,2,4,6,8)], ncol = 4, nrow =2))
dev.off()
rm(gg.list)
###################

### Shiny DEG correlation plots ####
# Plotting functions ----------------
# Correlation
corplot <- function(df){
  z <- length(levels(df$Type))
  ggplot(df, aes(x=LFC_HE, y=LFC_LE, color = Type)) +
    geom_point(size = 3, alpha = .4, shape = 16) +
    geom_point(data = df[which(df$Type %in% 'DEG Overlap'),], #2nd layer
               aes(x=LFC_HE, y=LFC_LE),
               shape = 16, size = 3, alpha = .5) +
    labs(title = paste0(substance," - DEG correlation"), 
         subtitle = paste(name,"vs", name.HE), 
         x = paste0('Log2FC (Control vs ',tail(resNames,1),")"), 
         y = paste0('Log2FC (Control vs ',k,")")) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = .4) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = .4) +
    theme_bw() + coord_equal(ratio=1) + 
    theme(aspect.ratio = 1, legend.position = 'bottom') +
    geom_rug(alpha = .5) +
    #geom_point(aes(x=df[df$Type %in% "DEG Overlap","LFC_HE"], y=df[df$Type %in% "DEG Overlap","LFC_LE"]), color="black") +
    #scale_colour_manual(values=c("dodgerblue2","darkgreen","red3")) #not really color blind friendly!
    #scale_shape_manual(values=c(15, 17, 16)) +
    #scale_colour_brewer(palette = "Dark2")+
    scale_colour_manual(values = if(z < 3){c("mediumblue","red2")}else{c("deepskyblue1","mediumblue","red2")})+
    #scale_colour_manual(values=c("gold1","darkorange1","red2"))+
    annotate("label",
             label = paste0("Pearson Cor = ",round(cor(df[df$Type %in% "DEG Overlap","LFC_HE"],df[df$Type %in% "DEG Overlap","LFC_LE"]),2)),
             color = "red2",
             x = (min(df$LFC_HE)), 
             y = (max(df$LFC_LE)),
             hjust=0, vjust=.3) + 
    annotate("label",
             label = paste0("Pearson Cor = ",round(cor(df$LFC_HE,df$LFC_LE),2)),
             x = (min(df$LFC_HE)), 
             y = (max(df$LFC_LE)-.45),
             hjust=0, vjust=.3) + 
    annotate("label",
             label = paste0("R2 = ",round(cor(df$LFC_HE,df$LFC_LE)^2,2)),
             x = (max(df$LFC_HE)), 
             y = (min(df$LFC_LE)),
             hjust=.75, vjust=0)+
    annotate("label",
             label = paste0("R2 = ",round(cor(df[df$Type %in% "DEG Overlap","LFC_HE"],df[df$Type %in% "DEG Overlap","LFC_LE"])^2,2)),
             color = "red2",
             x = (max(df$LFC_HE)), 
             y = (min(df$LFC_LE)+.45),
             hjust=.75, vjust=0)
}
# Venn Diagrams 
# venn.ls: List object containing gene IDs to plot Venn with
venplot <- function(venn.ls){plot(eulerr::euler(venn.ls,shape="ellipse",lty = 0),
                                  fills = list(fill = c("deepskyblue1","mediumblue"), alpha = .7),
                                  labels = list(col = "black", font = 4),
                                  #legend = list(col = "black", font = 4),
                                  quantities = TRUE, #to display values
                                  main = paste0(" [padj <",p,"; Stress: ",round(eulerr::euler(venn.ls,shape="ellipse")$stress,4),"; DiagError: ",round(eulerr::euler(venn.ls,shape="ellipse")$diagError,3),"]"))
}
# Plotting --------------------------
message(" 
 Start shiny DEG Correlation plotting ...")
pdf(file = paste0("Plots/",substance,"_DEG_cor.pdf"), onefile = T, bg = "transparent", #Background color
    width = 8*(length(resNames)-1), height = 14)

gg.list <- list() #empty list for storing plotting objects
# res
for(k in resNames[1:length(resNames)-1]){
  name <- paste0("res.",k)
  res <- get(name)
  name.HE <- paste0("res.",tail(resNames,1))
  res.HE <- get(name.HE) # Get highest exposure results
  
  LE <- setdiff(rownames(subset(res, padj <= p)),rownames(subset(res.HE, padj <= p))) #genes present in res and not in res.HE
  HE <- setdiff(rownames(subset(res.HE, padj <= p)),rownames(subset(res, padj <= p)))
  Ov <- intersect(rownames(subset(res.HE, padj <= p)),rownames(subset(res, padj <= p))) # Ov = overlap
  
  if(length(LE)!=0){
    df.LE <- data.frame(GeneID = LE, Type = "DEG in LE")
  }else{df.LE <- NULL}
  if(length(HE)!=0){
    df.HE <- data.frame(GeneID = HE, Type = "DEG in HE")
  }else{df.HE <- NULL}
  if(length(Ov)!=0){
    df.Ov <- data.frame(GeneID = Ov, Type = "DEG Overlap")
  }else{df.Ov <- NULL}
  
  df <- as.data.frame(rbind(df.LE,df.HE,df.Ov))
  rownames(df) <- df$GeneID
  df <- df[stringr::str_order(df$GeneID),]
  
  if (nrow(df) > 2){
    # Append LFC for HE and LE to df (there is propably an easier way to write this but can't figure it out... works for now)
    # HE
    res1 <- get(paste0("res.",tail(resNames,1))) 
    x <- which(rownames(res1) %in% df$GeneID == T) # Get position of the genes listed in df
    df1 <- res1[rownames(res1) %in% rownames(df),]
    df1 <- df1[stringr::str_order(rownames(df1)),]
    # Check if everything is correct 
    stopifnot(all(rownames(df1) == df$GeneID)) #Nice! 
    # Append HE LFC values
    df$LFC_HE <- df1$log2FoldChange
    
    # LE
    res1 <- get(paste0("res.",k))
    x <- which(rownames(res1) %in% df$GeneID == T) # Get position of the genes listed in df
    df1 <- res1[rownames(res1) %in% rownames(df),]
    df1 <- df1[stringr::str_order(rownames(df1)),]
    # Check if everything is correct 
    stopifnot(all(rownames(df1) == df$GeneID)) #Nice!
    # Append HE LFC values
    df$LFC_LE <- df1$log2FoldChange
    
    ## DEG Correlation plot
    gg <- corplot(df)
    gg.list[[paste0(k,".cor")]] <- gg
    
    ## Venn plot
    venn.ls <- list(LE = rownames(subset(res, padj <= p)),
                    HE = rownames(subset(res.HE, padj <= p)))
    gg <- venplot(venn.ls)
    gg.list[[paste0(k,".ven")]] <- gg
  }else{message(paste0(" Number of DEGs < 3 for: ",name," vs ",name.HE,"
 Skipped Correlation plotting"))}
}
if (nrow(df) > 2){
  # objects in gglist
  ggobj <- ls(gg.list)
  x <- ggobj[c(T, F)] #cor plots
  y <- ggobj[c(F, T)] #venn plots
  print(ggpubr::ggarrange(plotlist = gg.list[c(x,y)], ncol = length(resNames)-1, nrow =2))
}

# res.LFcut #
for(k in resNames[1:length(resNames)-1]){
  name <- paste0("res.",k,".LFcut")
  res <- get(name)
  name.HE <- paste0("res.",tail(resNames,1),".LFcut")
  res.HE <- get(name.HE) # Get highest exposure results
  
  LE <- setdiff(rownames(subset(res, padj <= p)),rownames(subset(res.HE, padj <= p))) #genes present in res and not in res.HE
  HE <- setdiff(rownames(subset(res.HE, padj <= p)),rownames(subset(res, padj <= p)))
  Ov <- intersect(rownames(subset(res.HE, padj <= p)),rownames(subset(res, padj <= p))) # Ov = overlap
  
  if(length(LE)!=0){
    df.LE <- data.frame(GeneID = LE, Type = "DEG in LE")
  }else{df.LE <- NULL}
  if(length(HE)!=0){
    df.HE <- data.frame(GeneID = HE, Type = "DEG in HE")
  }else{df.HE <- NULL}
  if(length(Ov)!=0){
    df.Ov <- data.frame(GeneID = Ov, Type = "DEG Overlap")
  }else{df.Ov <- NULL}
  
  df <- as.data.frame(rbind(df.LE,df.HE,df.Ov))
  rownames(df) <- df$GeneID
  df <- df[stringr::str_order(df$GeneID),]
  
  if (nrow(df) > 2){
    # Append LFC for HE and LE to df (there is propably an easier way to write this but can't figure it out... works for now)
    # HE
    res1 <- get(paste0("res.",tail(resNames,1))) 
    x <- which(rownames(res1) %in% df$GeneID == T) # Get position of the genes listed in df
    df1 <- res1[rownames(res1) %in% rownames(df),]
    df1 <- df1[stringr::str_order(rownames(df1)),]
    # Check if everything is correct 
    stopifnot(all(rownames(df1) == df$GeneID)) #Nice! 
    # Append HE LFC values
    df$LFC_HE <- df1$log2FoldChange
    
    # LE
    res1 <- get(paste0("res.",k))
    x <- which(rownames(res1) %in% df$GeneID == T) # Get position of the genes listed in df
    df1 <- res1[rownames(res1) %in% rownames(df),]
    df1 <- df1[stringr::str_order(rownames(df1)),]
    # Check if everything is correct 
    stopifnot(all(rownames(df1) == df$GeneID)) #Nice!
    # Append HE LFC values
    df$LFC_LE <- df1$log2FoldChange
    
    ## DEG Correlation plot
    gg <- corplot(df)
    gg.list[[paste0(k,".cor")]] <- gg
    
    ## Venn plot
    venn.ls <- list(LE = rownames(subset(res, padj <= p)),
                    HE = rownames(subset(res.HE, padj <= p)))
    gg <- venplot(venn.ls)
    gg.list[[paste0(k,".ven")]] <- gg
  }else{message(paste0(" Number of DEGs < 3 for: ",name," vs ",name.HE,"
 Skipped Correlation plotting"))}
}
if (nrow(df) > 2){
  # objects in gglist
  ggobj <- ls(gg.list)
  x <- ggobj[c(T, F)] #cor plots
  y <- ggobj[c(F, T)] #venn plots
  print(ggpubr::ggarrange(plotlist = gg.list[c(x,y)], ncol = length(resNames)-1, nrow =2))
}

# reslfs
for(k in resNames[1:length(resNames)-1]){
  name <- paste0("reslfs.",k)
  res <- get(name)
  name.HE <- paste0("reslfs.",tail(resNames,1))
  res.HE <- get(name.HE) # Get highest exposure results
  
  LE <- setdiff(rownames(subset(res, padj <= p)),rownames(subset(res.HE, padj <= p))) #genes present in res and not in res.HE
  HE <- setdiff(rownames(subset(res.HE, padj <= p)),rownames(subset(res, padj <= p)))
  Ov <- intersect(rownames(subset(res.HE, padj <= p)),rownames(subset(res, padj <= p))) # Ov = overlap
  
  if(length(LE)!=0){
    df.LE <- data.frame(GeneID = LE, Type = "DEG in LE")
  }else{df.LE <- NULL}
  if(length(HE)!=0){
    df.HE <- data.frame(GeneID = HE, Type = "DEG in HE")
  }else{df.HE <- NULL}
  if(length(Ov)!=0){
    df.Ov <- data.frame(GeneID = Ov, Type = "DEG Overlap")
  }else{df.Ov <- NULL}
  
  df <- as.data.frame(rbind(df.LE,df.HE,df.Ov))
  rownames(df) <- df$GeneID
  df <- df[stringr::str_order(df$GeneID),]
  
  if (nrow(df) > 2){
    # Append LFC for HE and LE to df (there is propably an easier way to write this but can't figure it out... works for now)
    # HE
    res1 <- get(paste0("res.",tail(resNames,1))) 
    x <- which(rownames(res1) %in% df$GeneID == T) # Get position of the genes listed in df
    df1 <- res1[rownames(res1) %in% rownames(df),]
    df1 <- df1[stringr::str_order(rownames(df1)),]
    # Check if everything is correct 
    stopifnot(all(rownames(df1) == df$GeneID)) #Nice! 
    # Append HE LFC values
    df$LFC_HE <- df1$log2FoldChange
    
    # LE
    res1 <- get(paste0("res.",k))
    x <- which(rownames(res1) %in% df$GeneID == T) # Get position of the genes listed in df
    df1 <- res1[rownames(res1) %in% rownames(df),]
    df1 <- df1[stringr::str_order(rownames(df1)),]
    # Check if everything is correct 
    stopifnot(all(rownames(df1) == df$GeneID)) #Nice!
    # Append HE LFC values
    df$LFC_LE <- df1$log2FoldChange
    
    ## DEG Correlation plot
    gg <- corplot(df)
    gg.list[[paste0(k,".cor")]] <- gg
    
    ## Venn plot
    venn.ls <- list(LE = rownames(subset(res, padj <= p)),
                    HE = rownames(subset(res.HE, padj <= p)))
    gg <- venplot(venn.ls)
    gg.list[[paste0(k,".ven")]] <- gg
  }else{message(paste0(" Number of DEGs < 3 for: ",name," vs ",name.HE,"
 Skipped Correlation plotting"))}
}
if (nrow(df) > 2){
  # objects in gglist
  ggobj <- ls(gg.list)
  x <- ggobj[c(T, F)] #cor plots
  y <- ggobj[c(F, T)] #venn plots
  print(ggpubr::ggarrange(plotlist = gg.list[c(x,y)], ncol = length(resNames)-1, nrow =2))
}

# reslfs.LFcut
for(k in resNames[1:length(resNames)-1]){
  name <- paste0("reslfs.",k,".LFcut")
  res <- get(name)
  name.HE <- paste0("reslfs.",tail(resNames,1),".LFcut")
  res.HE <- get(name.HE) # Get highest exposure results
  
  LE <- setdiff(rownames(subset(res, padj <= p)),rownames(subset(res.HE, padj <= p))) #genes present in res and not in res.HE
  HE <- setdiff(rownames(subset(res.HE, padj <= p)),rownames(subset(res, padj <= p)))
  Ov <- intersect(rownames(subset(res.HE, padj <= p)),rownames(subset(res, padj <= p))) # Ov = overlap
  
  if(length(LE)!=0){
    df.LE <- data.frame(GeneID = LE, Type = "DEG in LE")
  }else{df.LE <- NULL}
  if(length(HE)!=0){
    df.HE <- data.frame(GeneID = HE, Type = "DEG in HE")
  }else{df.HE <- NULL}
  if(length(Ov)!=0){
    df.Ov <- data.frame(GeneID = Ov, Type = "DEG Overlap")
  }else{df.Ov <- NULL}
  
  df <- as.data.frame(rbind(df.LE,df.HE,df.Ov))
  rownames(df) <- df$GeneID
  df <- df[stringr::str_order(df$GeneID),]
  
  if (nrow(df) > 2){
    # Append LFC for HE and LE to df (there is propably an easier way to write this but can't figure it out... works for now)
    # HE
    res1 <- get(paste0("res.",tail(resNames,1))) 
    x <- which(rownames(res1) %in% df$GeneID == T) # Get position of the genes listed in df
    df1 <- res1[rownames(res1) %in% rownames(df),]
    df1 <- df1[stringr::str_order(rownames(df1)),]
    # Check if everything is correct 
    stopifnot(all(rownames(df1) == df$GeneID)) #Nice! 
    # Append HE LFC values
    df$LFC_HE <- df1$log2FoldChange
    
    # LE
    res1 <- get(paste0("res.",k))
    x <- which(rownames(res1) %in% df$GeneID == T) # Get position of the genes listed in df
    df1 <- res1[rownames(res1) %in% rownames(df),]
    df1 <- df1[stringr::str_order(rownames(df1)),]
    # Check if everything is correct 
    stopifnot(all(rownames(df1) == df$GeneID)) #Nice!
    # Append HE LFC values
    df$LFC_LE <- df1$log2FoldChange
    
    ## DEG Correlation plot
    gg <- corplot(df)
    gg.list[[paste0(k,".cor")]] <- gg
    
    ## Venn plot
    venn.ls <- list(LE = rownames(subset(res, padj <= p)),
                    HE = rownames(subset(res.HE, padj <= p)))
    gg <- venplot(venn.ls)
    gg.list[[paste0(k,".ven")]] <- gg
  }else{message(paste0(" Number of DEGs < 3 for: ",name," vs ",name.HE,"
 Skipped Correlation plotting"))}
}
if (nrow(df) > 2){
  # objects in gglist
  ggobj <- ls(gg.list)
  x <- ggobj[c(T, F)] #cor plots
  y <- ggobj[c(F, T)] #venn plots
  print(ggpubr::ggarrange(plotlist = gg.list[c(x,y)], ncol = length(resNames)-1, nrow =2))
}
dev.off()
while (!is.null(dev.list()))  dev.off()
rm(x,y,ggobj,gg.list,venn.ls,df,df1,res,res1,res.HE,gg,df.LE,df.HE,df.Ov,LE,HE,Ov,name,name.HE)
###################################

### HEATMAPS of potential Biomarkers #####
# Plot heatmap with potential biomarkers genes
# Hence genes in the overlap between HE & LE
# if there are more than 2 treatment conditions. Take the overlap from 
# the highest treatment with the second highest for a heatmap;
# then the third highest and so on ...
## center rlog mean counts around the control's mean and scale for the gene's SD -----
centerForMean <- function(df){
  coldat <- coldata[colnames(df),]
  id <- rownames(coldat[coldat$Condition %in% 'Control',])
  ctrM <- apply(df[,id], 1, FUN = mean) #calcs the mean of control
  Sd <- apply(df, 1, FUN = sd) # calcs Sd of each row / gene
  (df - ctrM)/Sd # centers for the mean of control and scales for overall Sd of the row / gene
}
scaled.rld <- centerForMean(df.rld) #mean centered df for heatmap
## resort column order ---------------------------------------------------------------
id.order <- c()
for(i in levels(coldata$Substance)){
  for(k in levels(coldata$Condition)){
    id <- rownames(coldata[coldata$Condition %in% k & coldata$Substance %in% i,])
    id.order <- append(id.order,id)
  }
}
scaled.rld <- scaled.rld[,id.order]

## Get list of DE genes --------------------------------------------------------------
# non-shrunk
DEG <- list()
for (k in resNames){
  obj <- paste0("res.",k,".LFcut") # use this line to collect non-lfs LFcut DE genes
  x <- get(obj)
  x <- rownames(x[x$padj <= p,])
  DEG[[k]] <- x
}
# apeglm shrunk
DEGlfs <- list()
for (k in resNames){
  obj <- paste0("reslfs.",k,".LFcut")
  x <- get(obj)
  x <- rownames(x[x$padj <= p,])
  DEGlfs[[k]] <- x
}
# Get a list of potential biomarkers from overlaps between conditions
biomarker.ls <- list() # use this biomarker.ls as input for heatmap loop 
var <- length(resNames)
for (i in 1:(var-1)){
  # non shrunk
  name <- paste0(resNames[var],'.',resNames[i])
  biomarker.ls[[name]] <- intersect(DEG[[resNames[var]]],DEG[[resNames[i]]])
  # apeglm shrunk
  name <- paste0(resNames[var],'.',resNames[i],'.lfs')
  biomarker.ls[[name]] <- intersect(DEGlfs[[resNames[var]]],DEGlfs[[resNames[i]]])
  rm(name,obj,x)
}

## Create df with gene Symbols to correctly annotate the heatmap -------------------
# rowdata annotation file for heatmap with gene SYMBOLs
ref <- as.data.frame(subset(get(paste0('res.',resNames[1])), select = c('SYMBOL'))) #contains all ENSEMBL IDs
ref$GeneID <- as.factor(rownames(ref))
ref$SYMBOL <- factor(ref$SYMBOL)
# replace NAs in ref$SYMBOL with ENSEMBL IDs
#1) Add missing factor levels
Id <- rownames(ref[which(is.na(ref$SYMBOL)),]) # ENSEMBL IDs of missing SYMBOLs
levels <- levels(ref$SYMBOL)
for (i in c(1:length(Id))){
  levels[length(levels) + 1] <- Id[i]
}
ref$SYMBOL <- factor(ref$SYMBOL, levels = levels)
#2) replace NAs in ref$SYMBOL with ENSEMBL ID
ref[is.na(ref),'SYMBOL'] <- rownames(ref[which(is.na(ref$SYMBOL)),])
# length of levels should be = to nrow(ref); BUT can't be as multiple Gene IDs are assigned to a single Symbol
#length(levels) == nrow(ref)
# Get non unique SYMBOLS
#x <- ref[which(duplicated(ref$SYMBOL)),'SYMBOL']
#x[which(!is.na(x))]
#length(x[which(!is.na(x))])# many gene SYMBOLs not unique!!! 

## Define heatmap function --------------------------------------------------------
myheat2 <- function(df, colclust, rowclust, title, distM, Symbol){
  # Set anno colors
  Tanks <- levels(coldata$Tank)
  Conditions <- levels(coldata$Condition)
  col.Cond <- colorRampPalette(c('green',"gold1","indianred3"))(length(Conditions))
  col.Tank <- colorRampPalette(c("gray95","gray50"))(length(Tanks))
  ann_colors <- list(Condition = setNames(col.Cond, Conditions),
                     Tank = setNames(col.Tank, Tanks))
  colGaps <- c()
  for(k in Conditions[1:(length(Conditions)-1)]){
    x <- which(coldata[id.order,'Condition'] %in% k)
    x <- tail(x,1)
    colGaps <- append(colGaps,x)
  }
  # cluster fun
  callback <- function(hc, mat){
    sv = svd(t(mat))$v[,1]
    dend = reorder(as.dendrogram(hc), wts = sv)
    as.hclust(dend)
  }
  # set anno colors & breaks
  x <- 80 #length of color palette
  col <- colorRampPalette(c("mediumblue","white","red2"))(x) #defines color palette in heatmap
  s <- sd(df[,rownames(coldata[which(coldata$Condition %in% 'Control'),])]) # sets sd around 0 values from control to white
  m <- mean(df[,rownames(coldata[which(coldata$Condition %in% 'Control'),])])
  myBreaks <- c(seq(min(df), m-s, length.out=ceiling(x/2) + 1), 
                seq(m+s, max(df), length.out=floor(x/2)))
  # plot
  if(missing(colclust)){colclust = F}
  if(missing(rowclust)){rowclust = T}
  if(missing(distM)){distM = 'euclidean'} # euclidean, maximum, manhattan, ... 
  if(missing(title)){title = ''}
  if(missing(Symbol)){Symbol = F}
  if(Symbol == T){
    df1 <- df
    ref1 <- ref[row.names(df),]
    row.names(df1) <- ref1[row.names(ref1) %in% row.names(df),'SYMBOL']
  }else{df1 <- df}
  pheatmap::pheatmap(
    df1,
    angle_col = "90",
    gaps_col = if(colclust == T){NULL}else{colGaps},
    #cellwidth = 16,
    #cellheight = if(length(bioM) > 190){2}else{8}, #6
    clustering_callback = callback,
    treeheight_col = 30, #default 50
    border_color = NA, 
    cluster_rows= rowclust,
    cluster_cols= colclust,
    clustering_distance_rows = distM,
    clustering_distance_cols = distM, 
    clustering_method = "average",
    main = if(colclust == T){paste(title,'[mean cent]',distM)}else{paste(title,'[mean cent]')},
    show_rownames= T, 
    show_colnames = T,
    breaks = myBreaks,
    annotation_col = coldata[,c('Tank','Condition')],
    annotation_colors = ann_colors,
    color = col
  )
}
## Heatmap plotting loop ------------------------------------------------------------------------
#hist(as.data.frame(scaled.rld))
message(" 
 Start mean count centered heatmapping  ...")
for(k in names(biomarker.ls)){
  filename <- paste0(substance,'_Heatmap_',k)
  bioM <- biomarker.ls[[k]]
  if(length(bioM) < 3){
    message(paste0(
      " Number of DEGs for ",k," < 3. No heatmap was drawn."
    ))
  }else{
    title <- paste0(substance,': ',length(bioM),' DEG(',k,')')
    pdf(file = paste0("Plots/Heatmaps/",filename,".pdf"), #print directly to pdf
        width = 8, height = if(length(bioM) > 16){16}else{length(bioM)},
        onefile = T, bg = "transparent")
    myheat2(scaled.rld[bioM,], title=title)
    myheat2(scaled.rld[bioM,], title=title , Symbol = T)
    myheat2(scaled.rld[bioM,], title=title , colclust = T, distM = 'maximum')
    myheat2(scaled.rld[bioM,], title=title , colclust = T, distM = 'maximum', Symbol = T)
    dev.off()
  }
  while(!is.null(dev.list())) dev.off()
  rm(filename,bioM,title)
}
########################################

### ORA ...
### GSEA ...
### LDA ...
### PLS-DA ...


### Clear Env.####
rm(i,k,topgenes,xpos,dis)

message(paste0(" 
 Finished Wald's pairwise testing DESeq2 Analysis & Plotting.
 Saving R Environment in ",paste0(getwd(),"/.RData
 Saving ...
 ")))
save.image(paste0(getwd(),"/.RData"))
message(" .RData saved")

## Session Information ##
sink(paste0("SessionInfo_",substance,".txt"))
print(date())
print(devtools::session_info())
#print(xfun::session_info()) #use this instead when devtools breaks again ... 
sink()

# get back up to initial repo where analysis startet!
setwd("../")

message(" 
 Puh! What a ride! Finally END of SCRIPT :)
 ")
# Finally clear out memory (Helpful if you run this script in a loop to free some memory)
#rm(list = ls(all.names = TRUE)) #will clear all objects and hidden objects.

########### END OF SCRIPT ############