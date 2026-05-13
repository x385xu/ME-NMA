---
editor_options: 
  markdown: 
    wrap: 72
    
output: pdf_document
---

```{r setup, include=FALSE}
options(width = 80)
```

**Title**: Identifying Conditions Favouring Multiplicative Heterogeneity
Models in Network Meta-Analysis

**Author**: Xinlei Xu, Caitlin H. Daly, Audrey Béliveau

**Corresponding Author for Code**:

Name: Xinlei Xu Email:
[xinlei.xu\@uwaterloo.ca](mailto:xinlei.xu@uwaterloo.ca){.email}

**Configurations**:
R version 4.5.1 (2025-06-13)
Platform: aarch64-apple-darwin20
Running under: macOS Tahoe 26.3.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Toronto
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] tibble_3.3.0      ggpattern_1.3.1   gflnma_0.0.0.9000 ggplot2_4.0.3     dplyr_1.1.4      
[6] nmadb_1.2.0       netmeta_2.9-0     meta_8.2-1        metadat_1.4-0    

loaded via a namespace (and not attached):
 [1] gtable_0.3.6        xfun_0.54           CompQuadForm_1.4.4  lattice_0.22-7     
 [5] mathjaxr_1.8-0      tzdb_0.5.0          numDeriv_2016.8-1.1 vctrs_0.6.5        
 [9] tools_4.5.1         Rdpack_2.6.4        bitops_1.0-9        generics_0.1.4     
[13] pkgconfig_2.0.3     Matrix_1.7-3        RColorBrewer_1.1-3  S7_0.2.0           
[17] readxl_1.4.5        lifecycle_1.0.4     compiler_4.5.1      farver_2.1.2       
[21] stringr_1.5.2       textshaping_1.0.3   htmltools_0.5.8.1   RCurl_1.98-1.17    
[25] yaml_2.3.10         tidyr_1.3.1         pillar_1.11.1       nloptr_2.2.1       
[29] MASS_7.3-65         reformulas_0.4.1    boot_1.3-31         abind_1.4-8        
[33] nlme_3.1-168        tidyselect_1.2.1    digest_0.6.37       stringi_1.8.7      
[37] purrr_1.1.0         splines_4.5.1       magic_1.6-1         fastmap_1.2.0      
[41] grid_4.5.1          cli_3.6.5           metafor_4.8-0       magrittr_2.0.4     
[45] utf8_1.2.6          readr_2.1.5         withr_3.0.2         scales_1.4.0       
[49] rmarkdown_2.29      lme4_1.1-37         cellranger_1.1.0    ragg_1.5.0         
[53] hms_1.1.3           evaluate_1.0.5      knitr_1.50          rbibutils_2.3      
[57] rlang_1.1.6         Rcpp_1.1.0          glue_1.8.0          xml2_1.4.0         
[61] rstudioapi_0.17.1   minqa_1.2.8         R6_2.6.1            systemfonts_1.2.3  

**Code execution process**:

1.  Set folder ME_NMA as the current working directory.

2.  The following code files have no precedence for running.

    a)  Reproduce histograms of figure 1 and 6: run nmadb_analysis.R
    b)  Reproduce case studies 1 results: run casestudy1.R
    c)  Reproduce case studies 2~5 results: run casestudies2_5.R

**Source of data**: The original data used in this research can be found
at: <https://cran.r-project.org/src/contrib/Archive/nmadb/?C=D;O=A>
