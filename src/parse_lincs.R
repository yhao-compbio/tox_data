# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2021
## This script processes LINCS gene expression data after perturbation to obtain Tox21 compound-induced gene expression profile that matches input query time and dosage


## functions
library(cmapR);


## 0. Input arguments 
Args			<- commandArgs(T);
query_dose		<- Args[1];	# number that indicates treatment dosage of compounds for instances of interest 
query_time		<- Args[2];	# number that indicates measuring time (after treatment) for instances of interest
tox_sample_file		<- "data/tox21/tox21_assay_samples.tsv";			# name of Tox21 assay positive/negative compound file   
assay_meta_file		<- "downloads/tox21/tox21_assay_info.tsv"; 			# name of Tox21 assay info file  
library_meta_file	<- "downloads/tox21/tox21_10k_library_info.tsv";		# name of Tox21 compound library info file  
lincs_meta_file		<- "downloads/LINCS/siginfo_beta.txt";				# name of LINCS instance treatment meta file 
lincs_data_file		<- "downloads/LINCS/level5_beta_trt_cp_n720216x12328.gctx";	# name of LINCS level 5 compound-induced gene expression data file  
output_folder		<- "data/LINCS/tox21_assay";					# folder name of output processed file 

## 1. Process input Tox21 assay/compound data
# read in positive/negative compounds of Tox21 assays as data frame, select rows of data frame that contains positive compounds  
sample_df <- read.delim(file = tox_sample_file, header = T, sep = "\t");
outcome_id <- which(sample_df$assay_outcome == 1);
sample_df <- sample_df[outcome_id, c("assay_name", "compound_pubchem_cid")];
# read in Tox21 assay info as data frame, convert cell line names of Tox21 assays to lower case, create map data frame between assay name and cell line   
assay_df <- read.delim(file = assay_meta_file, header = T, sep = "\t");
assay_df$cell_line <- sapply(assay_df$cell_line, function(adcl){
	adcl_s <- strsplit(adcl, "-")[[1]];
	adcl1 <- paste(tolower(adcl_s), collapse = "");
	return(adcl1);	
});
assay_df <- assay_df[, c("protocol_name", "cell_line")];
# read in Tox21 compound info as data frame, convert compound names to lower case, creat map data frame between compound name PubChem ID
library_df <- read.delim(file = library_meta_file, header = T, sep = "\t");
library_df$SAMPLE_NAME <- tolower(library_df$SAMPLE_NAME);
library_df$PUBCHEM_CID <- sapply(library_df$PUBCHEM_CID, function(ldpc) paste("CID_", ldpc, sep = ""));
library_df <- unique(library_df[, c("SAMPLE_NAME", "PUBCHEM_CID")]);
# merge data frame of Tox21 assay positive compounds with assay info data frame and compond info data frame, combine cell line and compound names of each row    
sample_assay_df <- merge(sample_df, assay_df, by.x = "assay_name", by.y = "protocol_name");
sample_assay_df1 <- merge(sample_assay_df, library_df, by.x = "compound_pubchem_cid", by.y = "PUBCHEM_CID");
sample_assay_df1$compound_cell <- mapply(function(sadsn, sadcl) paste(sadsn, sadcl, sep = "@"), sample_assay_df1$SAMPLE_NAME, sample_assay_df1$cell_line);

## 2. Process input LINCS instance meta data to obtain the instances of Tox21 compounds that match the query time and dosage 
# read in LINCS instance meta data as data frame 
lincs_df <- read.delim(file = lincs_meta_file, header = T, sep = "\t");
# filter LINCS instances by perturbation type of compound treatment, select instance rows of interest  
trt_id <- which(lincs_df$pert_type %in% "trt_cp");
lincs_df <- lincs_df[trt_id, ];
# filter LINCS instances by quality control pass, select instance rows of interest    
qc_id <- which(lincs_df$qc_pass == 1);
lincs_df <- lincs_df[qc_id, ];
# filter LINCS instances by query compound dosage, select instance rows of interest  
query_idose <- paste(query_dose, "uM", sep = " ");
dose_id <- which(lincs_df$pert_idose %in% query_idose);
lincs_df <- lincs_df[dose_id, ];
# filter LINCS instances by query measuring time, select instance rows of interest  
query_itime <- paste(query_time, "h", sep = " ");
time_id <- which(lincs_df$pert_itime %in% query_itime);
lincs_df <- lincs_df[time_id, ];
# combine cell line and compound names of each instance row 
lincs_df$cmap_name <- tolower(lincs_df$cmap_name);
lincs_df$cell_iname <- tolower(lincs_df$cell_iname);
lincs_df$compound_cell <- mapply(function(ldcn, ldci) paste(ldcn, ldci, sep = "@"), lincs_df$cmap_name, lincs_df$cell_iname);
# merge LINCS instance data frame with Tox21 compound-assay data frame by combined cell line and compound names   
sample_assay_map <- merge(sample_assay_df1, lincs_df, by = "compound_cell");
# iterate by Tox21 assay of interest, select LINCS instances that perturb with compounds from the assay   
all_assays <- unique(sample_assay_map$assay_name);
all_assay_filter <- lapply(all_assays, function(aa){
	# select rows of interest from merged data frame to obtain LINCS instances that are relevant to the cell line of current Tox21 assay  
	aa_id <- which(sample_assay_map$assay_name %in% aa);
	sample_aa_map <- sample_assay_map[aa_id,];
	# iterate by compound from the selected data frame, select LINCS instance with the best data quality for each compound 
	sam_cids <- unique(sample_aa_map$compound_pubchem_cid);
	sc_filter <- sapply(sam_cids, function(sc){
		# select rows of interest from current LINCS instances data frame to obtain instances relevant to the current compound   
		sc_ids <- which(sample_aa_map$compound_pubchem_cid %in% sc);
		# select instance with the highest quality control score   
		sc_max <- max(sample_aa_map$cc_q75[sc_ids]);
		sc_max_id <- which(sample_aa_map$cc_q75[sc_ids] %in% sc_max)[[1]];
		max_id <- sc_ids[[sc_max_id]];
		return(max_id);
	});
	sample_aa_map_filter <- sample_aa_map[sc_filter,];
	return(sample_aa_map_filter);
});
# combine selected LINCS instance data frames from multiple Tox21 assays   
sample_assay_map1 <- do.call(rbind, all_assay_filter);
sample_assay_map2 <- sample_assay_map1[, c("assay_name", "compound_pubchem_cid", "sig_id")];
# obtain the unique instances from combined data frame, index the instances with numbers, add index of each instance as new column of the data frame  
all_sig_ids <- unique(sample_assay_map2$sig_id);
all_col_ids <- 1:length(all_sig_ids);
names(all_col_ids) <- all_sig_ids;
sample_assay_map2$col_id <- all_col_ids[sample_assay_map2$sig_id];
# write indexed data frame that contains LINCS instances of interest to output instance file  
output_map_file  <- paste(output_folder, "_siginfo_beta_", query_dose, "uM_", query_time, "h_map.tsv", sep = "");
write.table(sample_assay_map2, file = output_map_file, sep = "\t", col.names = T, row.names = F, quote = F);

## 3. Extract gene expression profile induced by selected LINCS instances of interest  
# use cmapR function to read in the meta data of all level 5 LINCS data columns    
lincs_col_meta <- read_gctx_meta(lincs_data_file, dim = "col");
# obtain the column indexes of selected LINCS instances of interest, use cmapR function to extract columns for whole level 5 data 
sample_assay_col_id <- which(lincs_col_meta$id %in% all_sig_ids);
sample_assay_lincs <- parse_gctx(lincs_data_file, cid = sample_assay_col_id); 
# write extracted gene expression data to output data file 
out_exp_mat <- sample_assay_lincs@mat;
out_exp_mat <- out_exp_mat[, all_sig_ids];
output_data_file <- paste(output_folder, "_level5_beta_trt_cp_", query_dose, "uM_", query_time, "h.tsv", sep = "");
write.table(out_exp_mat, file = output_data_file, row.names = T, col.names = F, sep = "\t", quote = F);
