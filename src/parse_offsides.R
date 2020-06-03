# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2019
## This script generates positive and negative samples for four types of organ toxicity

## 1. process drug-adverse event (AE) data 
# read in OFFSIDES data matrix
offsides <- read.csv(file = "downloads/offsides/OFFSIDES.csv.xz");
# obtain all adverse event names (for manual curation of AE~organ relationships)
all_ae_names <- unique(as.character(offsides$condition_concept_name));
writeLines(all_ae_names, "data/offsides/offsides_all_adverse_event_names.txt");
# convert data type of two columns
offsides$PRR <- as.numeric(as.character(offsides$PRR));
offsides$PRR_error <- as.numeric(as.character(offsides$PRR_error));

## 2. process AE~organ data
# read in AE~organ map 
ae_organ_map <- read.delim(file = "downloads/offsides/offsides_all_adverse_event_names_organ_map.tsv", header = T, sep = "\t");
non_void_id <- which(as.character(ae_organ_map$organ_type) != "");
ae_organ_map <- ae_organ_map[non_void_id,];
# merge offsides and organ map 
offsides_organ <- merge(offsides, ae_organ_map, by = "condition_concept_name");

## 3. calculate 95% confidence intervals (CIs) of PRR (lower bound - upper bound) 
# upper bound
offsides_organ$PRR_upper <- mapply(function(oop, oope){
	upper <- oop * exp(1.96 * oope);
	return(upper);
}, offsides_organ$PRR, offsides_organ$PRR_error);
# lower bound
offsides_organ$PRR_lower <- mapply(function(oop, oope){
	lower <- oop / exp(1.96 * oope);
	return(lower);
}, offsides_organ$PRR, offsides_organ$PRR_error);

## 4. obtain positive and negative samples 
# get all unique organ types
all_organ_types <- as.character(offsides_organ$organ_type);
unique_organ_types <- unique(all_organ_types);
# iterate by organ type
organ_type_samples <- lapply(unique_organ_types, function(uot){
	# get all PRR and CIs
	uot_id <- which(all_organ_types %in% uot);
	uot_df <- offsides_organ[uot_id,];
	# group by drugs
	all_uot_drugs <- tolower(as.character(uot_df$drug_concept_name));
	unique_uot_drugs <- unique(all_uot_drugs);
	# iterate by drug 
	unique_uot_drug_tox_label <- sapply(unique_uot_drugs, function(uud){
		# get all PRR and CIs of each drug
		uud_id <- which(all_uot_drugs %in% uud);
		uud_df <- uot_df[uud_id,];
		lower_len <- length(which(uud_df$PRR_lower > 2));
		# positive sample (lower bound of PPR > 2 for any associated adverse event)
		if(lower_len > 0)	return(1)
		# negative sample
		else	return(0)
	});
	# output format: data frame (column 1: organ type, column 2: drug name, column 3, toxicity label)
	uot_organ_type <- rep(uot, length(unique_uot_drugs));
	uot_output_df <- data.frame(organ_type = uot_organ_type, drug_name = unique_uot_drugs, toxicity_label = unique_uot_drug_tox_label);
	uot_output_od <- order(uot_output_df$toxicity_label, decreasing = T);
	uot_output_df <- uot_output_df[uot_output_od,];
	return(uot_output_df);
});
organ_type_sample_df <- do.call(rbind, organ_type_samples);

## 5. output 
# sample matrix 
write.table(organ_type_sample_df, file = "data/offsides/offsides_organ_toxicity_samples.tsv", sep = "\t", row.names = F, col.names = T, quote = F);
# all drugs involved
all_organ_type_drugs <- unique(as.character(organ_type_sample_df$drug_name));
writeLines(all_organ_type_drugs, "data/offsides/offsides_organ_toxicity_drugs.txt")


