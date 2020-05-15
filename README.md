# DESeq2Analysis
R script for differential gene expression (DGE) analysis with DESeq2 and data visualization. 
This script requires a non-normalized **CountMatrix.csv** and a **coldata.csv** file as input.
* A raw gene count matrix can be downloaded from the public ArrayExpress repository (i.e. E-MTAB-9056)
* The coldata file should contain at least the following columns:
  * _Condition_ (Treatment conditions, i.e. HighExposure, LowExposure, Control)
  * _Substance_ (Name of the tested Substance)
  * _Tank_ (Tank number from which spawning group samples were collected)
  * _Row names_ = Column names from *CountMatrix.csv*
* Execute this script in the same folder where files are stored
