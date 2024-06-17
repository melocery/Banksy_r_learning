## Banksy_r_learning
Repository for  
- learning and practicing Banksy
- reproducing the result of figure 3 in BANKSY paper

**Reference:**  
Singhal, V., Chou, N., Lee, J. et al. BANKSY unifies cell typing and tissue domain segmentation for scalable spatial omics data analysis. Nat Genet 56, 431–441 (2024). [https://doi.org/10.1038/s41588-024-01664-3](https://doi.org/10.1038/s41588-024-01664-3)  

## BANKSY
BANKSY is a method for clustering spatial omics data by augmenting the features of each cell with both an average of the features of its spatial neighbors along with neighborhood feature gradients. By incorporating neighborhood information for clustering, BANKSY is able to  

- improve cell-type assignment in noisy data  
- distinguish subtly different cell-types stratified by microenvironment  
- identify spatial domains sharing the same microenvironment

BANKSY is applicable to a wide array of spatial technologies (e.g. 10x Visium, Slide-seq, MERFISH, CosMX, CODEX) and scales well to large datasets.  

## code repository
[banksy-zenodo](https://github.com/jleechung/banksy-zenodo/tree/main): Scripts to reproduce BANKSY analyses.  
Focus on the **fig3-hypothalamus**.  
