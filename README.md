# placenta_cell_type_dna

Placental and Immune Cell DNA Methylation Reference Panel for Bulk Tissue Cell Composition Estimation in Epidemiological Studies

DOI: http://dx.doi.org/10.1080/15592294.2024.2437275

Authors: Kyle A. Campbell1, Justin A. Colacino2,3, John Dou1, Dana C. Dolinoy2,3, Sung Kyun Park1,2, Rita Loch-Caruso2, Vasantha Padmanabhan2,3,4,6, and Kelly M. Bakulski1

Affiliations:
1Epidemiology, School of Public Health, University of Michigan, Ann Arbor, MI 48109, USA.
2Environmental Health Sciences, School of Public Health, University of Michigan, Ann Arbor, MI 48109, USA.
3Nutritional Sciences, School of Public Health, University of Michigan, Ann Arbor, MI 48109, USA.
4Pediatrics, Michigan Medicine, University of Michigan, Ann Arbor, MI 48109, USA.
5Human Genetics, Michigan Medicine, University of Michigan, Ann Arbor, MI 48109, USA.
6Obstetrics and Gynecology, Michigan Medicine, University of Michigan, Ann Arbor, MI 48109, USA.

Corresponding author: Kelly M. Bakulski, PhD, 1415 Washington Heights, Ann Arbor, MI 48109, USA, bakulski@umich.edu, (734) 615-5899

Abstract
To distinguish DNA methylation (DNAm) from cell proportion changes in whole placental villous tissue research, we developed a robust cell type-specific DNAm reference to estimate cell composition. We collated new and existing cell type DNAm profiles quantified via Illumina EPIC or 450k microarrays. To estimate cell composition, we deconvoluted whole placental samples (n=36) with robust partial correlation based on the top 30 hyper- and hypomethylated sites identified per cell type. To test deconvolution performance, we evaluated RMSE in predicting principal components of DNAm variation in 204 external placental samples. We analyzed DNAm profiles (n=368,435 sites) from 12 cell types: cytotrophoblasts (n=18), endothelial cells (n=19), Hofbauer cells (n=26), stromal cells (n=21), syncytiotrophoblasts (n=4), six lymphocyte types (n=36), and nucleated red blood cells (n=11). Median cell composition was consistent with placental biology: 60.9% syncytiotrophoblast, 17.3% stromal, 8.8% endothelial, 3.7% cytotrophoblast, 3.7% Hofbauer, 1.7% nucleated red blood cells, and 1.2% neutrophils. Our expanded reference outperformed an existing reference in predicting DNAm variation (PC1, 15.4% variance explained, IQR=21.61) with cell composition estimates (RMSEP: 8.62 vs. 10.79, p-value<0.001). This cell type reference can robustly estimate cell composition from whole placental DNAm data to detect important cell types, reveal biological mechanisms, and improve causal inference.
Keywords: DNA methylation, placenta, blood cells, deconvolution, epigenetics

The reference matrix is available here as a .csv with read.csv("ref.csv") or R data object that can be loaded into R memory as "ref" with: > load("ref.rda") 
