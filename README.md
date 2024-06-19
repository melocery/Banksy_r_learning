## Banksy_r_learning
Repository for  
- learning and practicing Banksy
- reproducing the result of figure 3 in BANKSY paper

**Reference:**  
Singhal, V., Chou, N., Lee, J. et al. BANKSY unifies cell typing and tissue domain segmentation for scalable spatial omics data analysis. Nat Genet 56, 431–441 (2024). [https://doi.org/10.1038/s41588-024-01664-3](https://doi.org/10.1038/s41588-024-01664-3)  

## Scripts
- `process_hypothalamus_merfish.R`: reproduce the result of figure 3 using the provided banksy object.  
- `process_hypothalamus_merfish_allcells.R`: perform clustering using the original data.  

Need the (legacy) version of BANKSY (v 0.1.5) for recreating this analysis:  

```
remotes::install_github("prabhakarlab/Banksy@main")
# to install 0.1.5
# the actual installed version is 0.1.6 (Jun 17, 2024).
```

If the error message  

```
error in irlba::irlba(L, nv = n, nu = 0, maxit = iters) : 
  function 'as_cholmod_sparse' not provided by package 'Matrix'
```

appears, reinstall the packages `Matrix` and `irlba`:  

```
install.packages("Matrix")
install.packages("irlba")
```

## BANKSY
BANKSY is a method for clustering spatial omics data by augmenting the features of each cell with both an average of the features of its spatial neighbors along with neighborhood feature gradients. By incorporating neighborhood information for clustering, BANKSY is able to  

- improve cell-type assignment in noisy data  
- distinguish subtly different cell-types stratified by microenvironment  
- identify spatial domains sharing the same microenvironment

BANKSY is applicable to a wide array of spatial technologies (e.g. 10x Visium, Slide-seq, MERFISH, CosMX, CODEX) and scales well to large datasets.  

## Code Repository
[banksy-zenodo](https://github.com/jleechung/banksy-zenodo/tree/main): Scripts to reproduce BANKSY analyses.  
Focus on the **fig3-hypothalamus**.  

