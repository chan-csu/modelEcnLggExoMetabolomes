This repository contains supplementary files for the metabolic modeling analysis in the manuscript:  
Hove PR, Nealon NJ, Chan SHJ, Bover SM, Haberecht HB and Ryan EP. (Submitted) Metabolomics and proteomics of *L. rhamnosus* GG and *E. coli* Nissle  probiotic supernatants identify distinct pathways that mediate growth suppression of antimicrobial-resistant pathogens.

`ECN_20200729.mat/xml`  
Draft Cobra model reconstructed using KBase and a genome for *E. coli* Nissle 1917 (NCBI Reference Sequence: NZ_CP007799)  

`LGG_20200729.mat/xml`  
Draft Cobra model reconstructed using KBase and a genome for *Lactobacillus rhamnosus* GG (NCBI Reference Sequence: NC_013198.1)  

.mat files are directly readable by Matlab  
.xml files are in Systems Biology Markup Language (SBML) level 3 version 1 with flux balance constraints (fbc) version 2  

`metabolomicData.mat`
Matlab data file containing the fold change in metabolite abundance in ECN or LGG culture relative to the blank MRS media  
variables therein:  
- `mappedMets`: the metabolites mapped to the models  
- `metsFoldChangeEc`, `metsFoldChangeLg`: fold change of and mapping from the metabolites to the exchange reactions in the models  
- `metsFoldChangeHeader`: header for the columns in `metsFoldChangeEc/Lg`  

`pFBAconstrainedByExoMetabolomes.m`  
Matlab m-file for constructing the optimizing model for analysis from the metabolic models and metabolomic data.  
