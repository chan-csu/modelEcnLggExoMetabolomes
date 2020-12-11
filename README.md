ECN_20200729.mat/xml  
The draft model reconstructed using KBase and a genome for E. coli Nissle 1917 (NCBI Reference Sequence: NZ_CP007799)  

LGG_20200729.mat/xml  
Draft Cobra model reconstructed using KBase and a genome for Lactobacillus rhamnosus GG (NCBI Reference Sequence: NC_013198.1)  

.mat files are directly readable by Matlab  
.xml files are in Systems Biology Markup Language (SBML) level 3 version 1 with flux balance constraints (fbc) version 2  

metabolomicData.mat  
Matlab data file containing the fold change in metabolite abundance after culturing ECN or LGG relative to the blank media (MRS)  
variables therein:  
- `mappedMets` are the metabolites mapped to the models  
- `metsFoldChangeEc` and `metsFoldChangeLg` contains the fold change and the mapping between the mapped metabolites and the exchange reactions in the models  
- `metsFoldChangeHeader`: header for the columns in `metsFoldChangeEc/Lg`  

pFBAconstrainedByExoMetabolomes.m  
Matlab m-file for constructing the optimizing model for analysis from the metabolic models and metabolomic data.  
