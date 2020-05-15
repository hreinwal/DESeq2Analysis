# DESeq2Analysis
R script for differential gene expression (DGE) analysis with DESeq2 and data visualization. 
This script requieres an unormalized _*CountMatrix.csv*_ and a _*coldata.csv*_ file as input.
* Raw gene count matrix can be downloaded from the public ArrayExpress repository (i.e. E-MTAB-9056)
* coldata file should contain at minimum the following columns:
  * Condition (Treatment conditions, i.e. HighExposure, LowExposure, Control)
  * Substance (Name of the tested Substance)
  * Tank (Tank number from which spaning group samples were collected)
  * Rownames = Column names from *CountMatrix.csv*
