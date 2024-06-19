## Download the following data into the `/data` directory

### Mouse Hypothalamus MERFISH data
Download the file Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv into the data directory.  
Rename the file merfish_all_cells.csv or change the corresponding codes in `process_hypothalamus_merfish.R` and `process_hypothalamus_merfish_allcells.R`  
Download link: [merfish dataset](https://datadryad.org/stash/dataset/doi:10.5061/dryad.8t8s248)  

### Semi processed data files
Download to help with exact reproducibility.  
`banksyObj_provided.rds` used in `process_hypothalamus_merfish.R`: [banksyObj_provided.rds](https://www.dropbox.com/scl/fi/eq05c2gip0g61vc0i6n1e/banksyObj_provided.rds?rlkey=z7w2tywtn4h8jiapeymv0l54e&dl=0)  
After downloading the datasets, the script `process_hypothalamus_merfish.R` can be run to regenerate the paper analyses, with the results stored in the /out directory.