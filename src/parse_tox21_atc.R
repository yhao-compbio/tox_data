# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2021
## This script maps compounds from Tox21 HepG2 viability assay to their ATC drug categories


## 0. Input arguments 
tox21_file		<- "downloads/tox21/tox21_10k_library_info.tsv";	# name of Tox21 compound library info file  
compound_path_file	<- "/home/yunhao1/project/DTox/data/compound_target_probability_tox21_interpret/gamma-epsilon_0.001_0.1/compound_target_fingerprint_maccs_probability_tox21-rt-viability-hepg2-p1_whole_data.tsv_rt_25_ps_5_re_0_xs_20_al_0.5_ld_0.0001_model.pt_gamma-epsilon_0.001_0.1_path_relevance_pv.tsv";	# name of DTox HepG2 viability model interpretation result file
atc_file		<- "https://raw.githubusercontent.com/yhao-compbio/chemical/master/downloads/atc/ATC.csv";	# name of compound ATC code file  
atc_annotation_file	<- "https://raw.githubusercontent.com/yhao-compbio/chemical/master/downloads/atc/atc_annotation.tsv";	# name of ATC drug category annotation file
out_file		<- "data/tox21/tox21_10k_library_atc_code.tsv";	# name of output result file 

## 1. Process input Tox21 compound data  
# read in Tox21 compound info as data frame, convert compound names to lower case, creat map data frame between compound name PubChem ID 
tox21_df <- read.delim(file = tox21_file, header = T, sep = "\t");
tox21_df$cid <- sapply(tox21_df$PUBCHEM_CID, function(tdpc) paste("CID", tdpc, sep = "_"));
tox21_df$compound_name <- tolower(tox21_df$SAMPLE_NAME);
tox21_df <- tox21_df[, c("cid", "compound_name")];
tox21_df <- unique(tox21_df);
# read in DTox HepG2 viability model interpretation results of compounds as data frame, filter Tox21 compound info data frame by compounds invovled in HepG2 viability   
compound_path_df <- read.delim(file = compound_path_file, header = T, sep = "\t");
cpd_id <- which(tox21_df$cid %in% compound_path_df$cid);
tox21_df <- tox21_df[cpd_id, ];

## 2. Map compounds from Tox21 HepG2 viability assay to their ATC drug categories 
# read in ATC codes of compounds as data frame, obtain the ATC code of each compound 
atc_df <- read.csv(file = atc_file); 
atc_df$atc_id <- sapply(atc_df$Class.ID, function(adci){
	adci_s <- strsplit(adci, "/", fixed = T)[[1]];
	adci_atc <- adci_s[[length(adci_s)]];
	return(adci_atc);
});
atc_df$compound_name <- tolower(atc_df$Preferred.Label);
atc_df <- atc_df[, c("atc_id", "compound_name")];
atc_df <- unique(atc_df);
# merge ATC code data frame with HepG2 viability compound info data frame, obtain the ATC codes of compounds from HepG2 viability assay 
tox21_atc_df <- merge(tox21_df, atc_df, by = "compound_name");
# iterate by HepG2 viability compound, obtain detailed ATC code info of each code  
atc_s <- mapply(function(tadai){
	# first digit of ATC code 
	tadai_s <- strsplit(tadai, "")[[1]];
	tadai_p1 <- tadai_s[[1]];
	# first three digits of ATC code  
	tadai_p3 <- paste(tadai_s[1:3], collapse = "");
	# first five digits of ATC code 
	tadai_p5 <- paste(tadai_s[1:5], collapse = "");
	return(c(tadai_p1, tadai_p3, tadai_p5));
}, tox21_atc_df$atc_id);
atc_s <- t(atc_s);
colnames(atc_s) <- c("atc_id1", "atc_id3", "atc_id5");
tox21_atc_df1 <- data.frame(tox21_atc_df, atc_s);
# sort HepG2 viability compounds by first digit of ATC code  
tad_od <- order(tox21_atc_df1$atc_id1, tox21_atc_df1$atc_id3);
tox21_atc_df1 <- tox21_atc_df1[tad_od, ];
# read in annotations of ATC drug categories as data frame, merge it with HepG2 viability compound ATC code data frame to obtain the category annotation of each compound  
atc_annot_df <- read.delim(file = atc_annotation_file, header = T, sep = "\t");
tox21_atc_df2 <- merge(tox21_atc_df1, atc_annot_df, by = "atc_id1");
# write merged data frame to output result file 
tox21_atc_df2 <- tox21_atc_df2[, c("compound_name", "cid", "atc_id", "atc_id1", "atc_id1_annotation", "atc_id3", "atc_id5")];
write.table(tox21_atc_df2, file = out_file, sep = "\t", col.names = T, row.names = F, quote = F);
