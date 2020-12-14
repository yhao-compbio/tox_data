# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2019
## This script generates positive and negative samples for four types of organ toxicity

offsides_file	<- "downloads/offsides/OFFSIDES.csv.xz";
output_file	<- "data/offsides/offsides_all_adverse_event";

## 1. Process drug-adverse event (AE) data 
# read in OFFSIDES data matrix
offsides <- read.csv(file = offsides_file);
offsides <- unique(offsides);

## 2. Compute lower and upper bound of PRR 
# upper bound
offsides$PRR_upper <- mapply(function(op, ope){
	upper <- op * exp(1.96 * ope);
	return(upper);
}, as.numeric(offsides$PRR), as.numeric(offsides$PRR_error));
# lower bound
offsides$PRR_lower <- mapply(function(op, ope){
	lower <- op / exp(1.96 * ope);
	return(lower);
}, as.numeric(offsides$PRR), as.numeric(offsides$PRR_error));

## 3. Obtain positive and negative samples 
# get all unique organ types
all_aes <- offsides$condition_concept_name;
unique_aes <- unique(all_aes);
# iterate by adverse event 
ae_samples <- lapply(unique_aes, function(ua){
	# get all statistics of the adverse event  
	ua_id <- which(all_aes %in% ua);
	ua_df <- offsides[ua_id,];
	# identify positive and negative samples
	pos_id <- which(ua_df$PRR_lower > 1);
	pos_drugs <- unique(ua_df$drug_concept_name[pos_id]);
	ua_pos_len <- length(pos_drugs);
	neg_id <- intersect(which(ua_df$PRR_lower < 1), which(ua_df$PRR_upper > 1));
	neg_drugs <- unique(ua_df$drug_concept_name[neg_id]);
	ua_neg_len <- length(neg_drugs);
	sample_drugs <- c(pos_drugs, neg_drugs);
	ua_sample_len <- length(sample_drugs);
	if ((ua_pos_len < 100) || (ua_neg_len < 100))	return(logical(0))	
	else{
		# adverse event name
		ua_s1 <- strsplit(tolower(ua), "/", fixed = T)[[1]]
		ua1 <- paste(ua_s1, collapse = " ")
		ua_s2 <- strsplit(ua1, " ", fixed = T)[[1]]
		ua_name <- paste(ua_s2, collapse = "_")
		ua_name_vec <- rep(ua_name, ua_sample_len)
		# drug name 
		drug_vec <- sapply(sample_drugs, function(sds){
			# remove ',', '.', '(', ')', and "/"
			sds_s <- strsplit(tolower(sds), "", fixed = T)[[1]]
			c_id <- which(sds_s %in% c(",", ".", "(", ")", "/"))
			o_id <- setdiff(1:length(sds_s), c_id)
			sds_s <- sds_s[o_id]
			sds1 <- paste(sds_s, collapse = "")
			# substitue ' ' with '_'
			sds1_s <- strsplit(sds1, " ", fixed = T)[[1]]
			sds2 <- paste(sds1_s, collapse = "_")
			return(sds2)
		});
		# 
		ua_label_vec <- c(rep(1, length(pos_drugs)), rep(0, length(neg_drugs)));
		# output format: data frame (column 1: adverse event name, column 2: drug name, column 3, toxicity label) 
		ua_out_df <- data.frame(adverse_event_name = ua_name_vec, drug_name = drug_vec, toxicity_label = ua_label_vec)
		return(ua_out_df)
	}
});
# remove AEs that do not meet the minimum sample requirement 
ae_sample_len <- sapply(ae_samples, length);
ae_samples <- ae_samples[ae_sample_len > 0];
names(ae_samples) <- sapply(ae_samples, function(as) as$adverse_event_name[[1]]);
ae_sample_df <- do.call(rbind, ae_samples);
# make summary statistics table of all adverse events  
ae_sample_summary <- mapply(function(as){
	all_len <- nrow(as);
	pos_len <- sum(as$toxicity_label);
	neg_len <- all_len - pos_len;
	return(c(all_len, pos_len, neg_len));
}, ae_samples);
ae_sample_summary_df <- data.frame(names(ae_samples), t(ae_sample_summary));
colnames(ae_sample_summary_df) <- c("adverse_event", "n_all_samples", "n_positive_samples", "n_negative_samples");
assd_od <- order(ae_sample_summary_df$n_all_samples, decreasing = T);
ae_sample_summary_df <- ae_sample_summary_df[assd_od, ];

## 5. output 
# sample matrix 
write.table(ae_sample_df, file = paste(output_file, "_samples.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F);
# all drugs involved
all_ae_drugs <- unique(as.character(ae_sample_df$drug_name));
writeLines(all_ae_drugs, paste(output_file, "_drugs.txt", sep = ""));
# data summary table 
write.table(ae_sample_summary_df, file = paste(output_file, "_summary.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F);
