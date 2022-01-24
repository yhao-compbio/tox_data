# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2021
## This script identifies positive/negative DSSTox compounds for drug-induced kidney injury (DIKI) phenotypes, and tests the association between each DIKI phenotype and HEK293 cell viability


## 0. Input arguments 
offsides_file		<- "data/offsides/offsides_all_adverse_event_samples.tsv";	# name of OFFSIDES adverse event compound file 
diki_file		<- "data/KidneyTox/diki_phenotype_terms.tsv";	# name of DIKI phenotype term name file 
dsstox_name_file	<- "/home/yunhao1/project/ComptoxAI/data/Dsstox_CAS_number_name.tsv";	# name of DSSTox compound info file  
dsstox_feature_file	<- "/home/yunhao1/project/ComptoxAI/data/dsstox_new_compounds_fingerprint_maccs.tsv";	# name of DSSTox compound fingerprint file
tox21_file		<- "downloads/tox21/tox21_10k_library_info.tsv";	# name of Tox21 compound library info file
tox21_sample_file	<- "data/tox21/tox21_assay_samples.tsv";	# name of input Tox21 assay to inactive/active compound file
query_assay_name	<- "tox21-rt-viability-hek293-p1";	# name of Tox21 HEK293 viability asssay
compound_path_file	<- "/home/yunhao1/project/DTox/data/compound_target_probability_tox21_interpret/gamma-epsilon_0.001_0.1/compound_target_fingerprint_maccs_probability_tox21-rt-viability-hek293-p1_whole_data.tsv_rt_14_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model.pt_gamma-epsilon_0.001_0.1_path_relevance_pv.tsv";	# name of input HEK293 cell vibiality DTox model interpretation result file
out_folder		<- "data/KidneyTox/diki_phenotype";	# folder name of output result files

## 1. Process input OFFSIDES adverse event data  
# read in positive/negative compounds of OFFSIDES adverse events as data frame, convert spaces in compound names to underscores
offsides_df <- read.delim(file = offsides_file, header = T, sep = "\t");
offsides_df$drug_name <- sapply(offsides_df$drug_name, function(oddn) paste(strsplit(oddn, "_")[[1]], collapse = " "));
# read in adverse event names of DIKI phenotypes as data frame 
diki_term_df <- read.delim(file = diki_file, header = T, sep = "\t");
# iterate by DIKI phenotype, obtain the compounds associated with each phenotype 
diki_drugs_list <- lapply(diki_term_df$term_name, function(dtdtn){
	dtdtn_id <- which(offsides_df$adverse_event_name %in% dtdtn);
	dtdtn_df <- offsides_df[dtdtn_id, ];
	return(dtdtn_df);
});
names(diki_drugs_list) <- diki_term_df$term_name;

## 2. Obtain the positive/negative DSSTox compounds of DIKI phenotypes 
# read in DSSTox compound info line by line as strings
dsstox_name_lines <- readLines(dsstox_name_file);
N_lines <- length(dsstox_name_lines);
# iterate by line, extract compound name and DSSTox ID from string of each line 
dsstox_name_df <- mapply(function(dnl){
	dnl_s <- strsplit(dnl, "\t")[[1]];
	return(dnl_s[2:3])
}, dsstox_name_lines[2:N_lines]);
dsstox_name_df <- data.frame(t(dsstox_name_df));
colnames(dsstox_name_df) <- c("dsstox_substance_id", "preferred_name");
dsstox_name_df$preferred_name <- tolower(dsstox_name_df$preferred_name);
# read in fingerprints of DSS compounds as data frame, filter DSSTox compound info data frame by compounds with computed fingerprints
dsstox_feature_df <- read.delim(file = dsstox_feature_file, header = T, sep = "\t");
dfd_id <- which(dsstox_name_df$dsstox_substance_id %in% rownames(dsstox_feature_df));
dsstox_name_df <- dsstox_name_df[dfd_id, ];
# iterate by DIKI phenotype, obtain the positive/negative DSSTox compounds of current phenotype 
dsstox_diki_list <- lapply(diki_drugs_list, function(ddl){
	# merge data frame of positive/negative compounds of current phenotype with DSSTox compound info to obtain the positive/negative DSSTox compounds of current phenotype
	ddl1 <- merge(dsstox_name_df, ddl, by.x = "preferred_name", by.y = "drug_name");
	ddl1_od <- order(ddl1$toxicity_label, decreasing = T);
	ddl1 <- ddl1[ddl1_od, c("adverse_event_name", "preferred_name", "dsstox_substance_id", "toxicity_label")];
	return(ddl1);
});
# aggregate positive/negative compound info data frames from all DIKI phenotypes, write to output file  
dsstox_diki_df <- do.call(rbind, dsstox_diki_list);
dsstox_diki_file <- paste(out_folder, "_dsstox_compounds.tsv", sep = "");
write.table(dsstox_diki_df, file = dsstox_diki_file, sep = "\t", row.names = F, col.names = T, quote = F);

## 3. Test the association between HEK293 viability and DIKI phenotypes by Fisher's exact test   
# read in Tox21 compound info as data frame, convert compound names to lower case, creat map data frame between compound name PubChem ID 
tox21_df <- read.delim(file = tox21_file, header = T, sep = "\t");
tox21_df$cid <- sapply(tox21_df$PUBCHEM_CID, function(tdpc) paste("CID", tdpc, sep = "_"));
tox21_df$compound_name <- tolower(tox21_df$SAMPLE_NAME);
tox21_df <- tox21_df[, c("cid", "compound_name")];
tox21_df <- unique(tox21_df);
# read in Tox21 assay active/inactive compounds as data frame, identify the active/inactive compounds of query HEK293 viability assay 
tox21_sample_df <- read.delim(file = tox21_sample_file, header = T, sep = "\t");
query_sample_id <- which(tox21_sample_df$assay_name %in% query_assay_name);
query_sample_df <- tox21_sample_df[query_sample_id, ];
# merge HEK293 viability assay active/inactive compound data frame with Tox21 compound info data frame to obtain HEK293 viability assay active/inactive compound info 
query_sample_df <- merge(query_sample_df, tox21_df, by.x = "compound_pubchem_cid", by.y = "cid");
# iterate by DIKI phenotype, test the association between HEK293 viability and each phenotype  
query_sample_association <- mapply(function(ddl){
	# merge HEK293 viability assay active/inactive compound info data frame with data frame of positive/negative compounds of current phenotype 
	query_sample_ddl <- merge(query_sample_df, ddl, by.x = "compound_name", by.y = "drug_name");
	# build 2 by 2 contingency table to test the association between HEK293 viability and current phenotype 
	qsd_assay_id <- which(query_sample_ddl$assay_outcome == 1);
	qsd_ae_id <- which(query_sample_ddl$toxicity_label == 1);
	N11 <- length(intersect(qsd_assay_id, qsd_ae_id));
	N12 <- length(qsd_assay_id) - N11;
	N21 <- length(qsd_ae_id) - N11;
	N22 <- nrow(query_sample_ddl) - length(qsd_assay_id) - N21;
	ddl_conti_mat <- matrix(c(N11, N12, N21, N22), 2, 2,  byrow = T);
	# compute test statistics: odds ratio and its 95% confidence interval 
	ddl_fisher <- fisher.test(ddl_conti_mat);
	return(c(ddl_fisher$estimate, ddl_fisher$conf.int));
}, diki_drugs_list);
# aggregate association test result data frames from all phenotypes, write to output result file 
diki_term_viab_assoc_df <- data.frame(diki_term_df, t(query_sample_association));
colnames(diki_term_viab_assoc_df) <- c(colnames(diki_term_df), "viability_odds_ratio", "viability_odds_ratio_95ci_lower", "viability_odds_ratio_95ci_upper");
diki_term_viab_assoc_file <- paste(out_folder, "_", query_assay_name, "_association.tsv", sep = "");
write.table(diki_term_viab_assoc_df, file = diki_term_viab_assoc_file, sep = "\t", row.names = F, col.names = T, quote = F);

## 4. Obtain the positive compound info of DIKI phenotypes 
# read in DTox HEK293 viability model interpretation results of compounds as data frame, filter Tox21 compound info data frame by compounds invovled in HEK293 viability
compound_path_df <- read.delim(file = compound_path_file, header = T, sep = "\t");
cpd_id <- which(tox21_df$cid %in% compound_path_df$cid);
tox21_df <- tox21_df[cpd_id, ]; 
# iterate by DIKI phenotype,  obtain the positive compound info of each phenotype
diki_tox21_df_list <- mapply(function(ddl, dtdtn, dtdts){
	# obtain the positive compounds of current phenotype, store compounds in a data frame along with phenotype info  
	ddl_pos_id <- which(ddl$toxicity_label == 1);
	dp_drugs <- ddl$drug_name[ddl_pos_id];
	dp_len <- length(ddl_pos_id);
	dp_df <- data.frame(rep(dtdtn, dp_len), rep(dtdts, dp_len), dp_drugs);
	colnames(dp_df) <- c("phenotype_name", "phenotype_name_short", "compound_name");
	# merge with Tox21 compound info data frame to obtain the PubChem IDs of positive compounds of current phenotype  
	dp_tox21_df <- merge(dp_df, tox21_df, by = "compound_name");
	return(dp_tox21_df);
}, diki_drugs_list, diki_term_df$term_name, diki_term_df$term_short, SIMPLIFY = F);
# aggregate positive compound info data frames from all phenotypes, write to output result file
diki_tox21_df <- do.call(rbind, diki_tox21_df_list);
diki_tox21_df <- diki_tox21_df[, c("phenotype_name", "phenotype_name_short", "cid", "compound_name")];
diki_tox21_file <- paste(out_folder, "_tox21_compounds.tsv", sep = ""); 
write.table(diki_tox21_df, file = diki_tox21_file, sep = "\t", row.names = F, col.names = T, quote = F);

