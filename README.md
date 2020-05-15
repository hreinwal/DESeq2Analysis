# DESeq2Analysis for full RNASeq Ecotoxicogenomic studies using zebrafish (*Danio rerio*)
R script for differential gene expression (DGE) analysis with [*DESeq2* (Love et al, 2014)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) for ecotoxicological testing on zebrafish embryos. 
This script requires a non-normalized **CountMatrix.csv** and a **coldata.csv** file as input.
* A raw gene count matrix can be downloaded from the public ArrayExpress repository (i.e. E-MTAB-9056)
* The coldata file should contain at least the following columns:
  * _Condition_ (Treatment conditions, i.e. HighExposure, LowExposure, Control)
  * _Substance_ (Name of the tested Substance)
  * _Tank_ (Tank number from which spawning group samples were collected)
  * _Row names_ = Column names from *CountMatrix.csv*
* Execute this script in the same folder where files are stored

This script will run DESeq2 with pairwise Wald's t-test with [*IHW* (Ignatiadis et al, 2016)](https://www.nature.com/articles/nmeth.3885) when correcting p-values for multiple testing after Benjamini-Hochberg. Effect size shrinking is applied through [*apeglm* (Zhu et al, 2019)](https://academic.oup.com/bioinformatics/article/35/12/2084/5159452) and if applied, the effect size cutoff (LFcut) is determined as the 90%-quantile of absolute log2-fold changes. The output will be annotated using R's *AnnotationDbi* with the *org.Dr.eg.db*. The script is designed to analyze one tested substance at a time.

Not the most beautiful code in world but it does the job ;)
