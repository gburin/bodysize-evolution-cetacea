# The dynamic adaptive landscape of cetacean body size

Authors: [Gustavo Burin\*](mailto:gustavo.burin@nhm.ac.uk), [Travis Park\*](mailto:t.park@nhm.ac.uk), Tamora D. James, Graham J. Slater, [Natalie Cooper](mailto:natalie.cooper@nhm.ac.uk)

[*] These authors contributed equally to this work

This repository contains all the code and data used in the manuscript.

![](output/figs_manuscript/fig2.png)

## Data

* `/data` contains the cleaned data used in the analyses, including the phylogeny.

All raw and cleaned data are also available from the [Natural History Museum Data Portal (data.nhm.ac.uk)](https://doi.org/10.5519/vmbrpkuq). 

If you use the cleaned data please cite as follows: 
> Gustavo Burin; Travis Park; Tamora D James; Graham J Slater; Natalie Cooper (2022). Dataset: Data from the “Back to the water” project. Natural History Museum Data Portal (data.nhm.ac.uk). https://doi.org/10.5519/vmbrpkuq.

**Please cite the original sources of the data where possible.**

The phylogeny comes from Lloyd & Slater (2021). **Please cite this paper if you use these data!**

Lloyd, GT, Slater, GJ, 2021. A Total-Group Phylogenetic Metatree for Cetacea and the Importance of Fossil Data in Diversification Analyses. Syst. Biol. 70(5): 922- 939.

-------
## Analyses
The analysis code is divided into subfolders, one for each dataset used. The results are summarized into two `.Rmd` files that generate a PDF with all plots (both from the main text and the supplementary material). A script to generate the base elements for the main figures is also available. The results from bayou are provided into subfolders that follow the same structure as the scripts used to generate them, and the other results (e.g. phylogenetic imputation) are made available through `.RDS` files.

-------
## Other folders

* `output/figs_manuscript` contains the main figures
* `bayou_sensitivity_figs` contains the compiled [PDF](https://github.com/gburin/bodysize-evolution-cetacea/raw/master/bayou_sensitivity_figs/sensitivity_figures.pdf) and individual figures with the results for the _bayou_ sensitivity analyses

-------
## Session Info

```
R version 4.2.2 (2022-10-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Manjaro Linux

Matrix products: default
BLAS:   /usr/lib/libblas.so.3.11.0
LAPACK: /usr/lib/liblapack.so.3.11.0

locale:
 [1] LC_CTYPE=en_GB.utf8        LC_NUMERIC=C              
 [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.utf8     
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.utf8    
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] Rphylopars_0.3.9 bayou_2.2.0      coda_0.19-4      geiger_2.0.10   
[5] remotes_2.4.2    phytools_1.2-0   maps_3.4.1       ape_5.6-2       

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.10             subplex_1.8             mvtnorm_1.1-3          
 [4] lattice_0.20-45         phylolm_2.6.2           listenv_0.9.0          
 [7] tidyr_1.3.0             assertthat_0.2.1        digest_0.6.31          
[10] foreach_1.5.2           utf8_1.2.3              parallelly_1.34.0      
[13] R6_2.5.1                backports_1.4.1         ggplot2_3.4.1          
[16] pillar_1.8.1            rlang_1.0.6             phangorn_2.11.1        
[19] Matrix_1.5-3            combinat_0.0-8          splines_4.2.2          
[22] igraph_1.4.0            munsell_0.5.0           broom_1.0.3            
[25] compiler_4.2.2          numDeriv_2016.8-1.1     Deriv_4.1.3            
[28] pkgconfig_2.0.3         microbenchmark_1.4.9    mnormt_2.1.1           
[31] optimParallel_1.0-2     denstrip_1.5.4          globals_0.16.2         
[34] tidyselect_1.2.0        tibble_3.1.8            expm_0.999-7           
[37] quadprog_1.5-8          codetools_0.2-19        fitdistrplus_1.1-8     
[40] future_1.31.0           fansi_1.0.4             dplyr_1.1.0            
[43] MASS_7.3-58.2           grid_4.2.2              nlme_3.1-162           
[46] gtable_0.3.1            lifecycle_1.0.3         magrittr_2.0.3         
[49] scales_1.2.1            future.apply_1.10.0     cli_3.6.0              
[52] scatterplot3d_0.3-42    doBy_4.6.16             vctrs_0.5.2            
[55] generics_0.1.3          fastmatch_1.1-3         deSolve_1.34           
[58] iterators_1.0.14        tools_4.2.2             glue_1.6.2             
[61] purrr_1.0.1             plotrix_3.8-2           parallel_4.2.2         
[64] survival_3.5-3          colorspace_2.1-0        clusterGeneration_1.3.7
```
