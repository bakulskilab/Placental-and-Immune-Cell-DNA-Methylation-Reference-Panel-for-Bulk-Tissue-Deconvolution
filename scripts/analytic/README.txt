Scripts for "Placental and Immune Cell DNA Methylation Reference Panel for Bulk Tissue Cell Composition Estimation in Epidemiological Studies"

custom_helper_scripts Contains an .R script file that contains various helper functions employed throughout different code files
	10_custom_functions.R

01_qc Contains all the dataset quality control, preprocessing, and formatting for raw DNA methylation data for each study used in this analysis
	flow_sort_blood_epic_qc.rmd Quality control and formatting for FlowSorted.Blood.EPIC Bioconductor package
	merge_and_qc_datasets_campbell.rmd Quality control and formatting for the samples collected in this study
	merge_and_qc_datasets_nrbc_bakulski.rmd Quality control and formatting for the Bakulski et al., 2016 nucleated red blood cell samples collected in FlowSorted.CordBloodCombined.450k Bioconductor package
	merge_and_qc_datasets_nrbc_degoede.rmd Quality control and formatting for the deGoede et al., 2016 nucleated red blood cell samples collected in FlowSorted.CordBloodCombined.450k Bioconductor package
	merge_and_qc_datasets_yuan.rmd Quality control and formatting for the placental Yuan et al., 2021 study

02_dataset_management Contains the code used to merge quality controlled datasets for downstream analysis
	concatenate_qced_datasets_all_tissues.rmd Creates primary analytic datset by merging DNA methylation data and subsets to the sites common to all samples from all placental, immune cell, and nucleated red blood cell tissues
	concatenate_qced_datasets_placenta_only.rmd Merges DNA methylation data for the placental datasets only, used only for low_purity_testing

03_low_purity_testing Contains the code used to evaluate cell type fractions for inclusion/exclusion, as well as data visualizations of clustering and related processes
	deconvolution_sample_pruning_unified_qc.rmd

04_descriptive Contains the code necessary to create the analytic sample descriptive table
	sample_description_overlap.rmd

05_differential_methylation Contains the code used to test differential methylation between cell types and gene ontology enrichment
	compare_DNAm_cell_types_450k_beta.rmd Script to perform differential methylation across all cell type samples used in the reference panel
	differential_dnam_enrichment_450k_all_cell_types.rmd Script to perform biological process ontology enrichment analysis of differential DNA methylation results

06_deconvolution_Yuan_Campbell_whole_tissue Contains the code used to create and apply our deconvolution reference to the whole tissue samples collected in this study and Yuan et al., 2021
	deconvolution_450k_analytic.rmd

07_deconvolution_performance_compare Contains the the cross-validation code used to evaluate and test our deconvolution reference against the existing planet Bioconductor package
	xval_glm.rmd

	