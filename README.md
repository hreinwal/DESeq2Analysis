# DESeq2Analysis
R script for differential gene expression (DGE) analysis with DESeq2 and data visualization. 
This script requires a non-normalized _*CountMatrix.csv*_ and a _*coldata.csv*_ file as input.
* A raw gene count matrix can be downloaded from the public ArrayExpress repository (i.e. E-MTAB-9056)
* The coldata file should contain at least the following columns:
  * Condition (Treatment conditions, i.e. HighExposure, LowExposure, Control)
  * Substance (Name of the tested Substance)
  * Tank (Tank number from which spawning group samples were collected)
  * Row names = Column names from *CountMatrix.csv*
* Execute this script in the same folder where files are stored
