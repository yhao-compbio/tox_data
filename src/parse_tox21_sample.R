# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2021
## This script parses Tox21 screening dataset and generates active/inactive compounds for all assays from Tox21 


## 0. Input arguments
Args			<- commandArgs(T);
min_sample_size		<- as.integer(Args[1]);				# minimum sample size requirement for active/inactive compounds of each assay  
output_file		<- Args[2];					# output file name
tox21_data_files	<- "downloads/tox21/tox21_data_files.txt";	# tox21 data files

## 1. Obtain Tox21 assay information from data file names 
# read in data file names  
tox_files <- readLines(tox21_data_files);
# obtain assay information: whether measures viability, whether measures agnoists or antagonists
tox_hyper <- mapply(function(tf){
	# parse to get assay target name  
	tf_name <- strsplit(tf, "/", fixed = T)[[1]][[3]];
	tfn_s <- strsplit(tf_name, "-", fixed = T)[[1]];
	# return different results column for viability and non-viability assays 
	if("viability" %in% tfn_s) result_col <- "CHANNEL_OUTCOME"	
	else	result_col <- "ASSAY_OUTCOME"
	# return different outcome key words for agnoist and antagonist assays 
	if("agonist" %in% tfn_s)	result_type <- "^active agonist$"
	else if("antagonist" %in% tfn_s)	result_type <- "^active antagonist$"
	else result_type <- "^active "
	# return info 
	return(c(tf_name, result_col, result_type));
}, tox_files);

## 2. Extract active/inactive compounds from each assay data file  
# iterate by assay  
outcome_list <- mapply(function(tf, th1, th2, th3){
	# read data file line by line (read.delim can read in all the lines for some assays) 
	tox_lines <- readLines(tf);
	if(length(tox_lines) == 1)	return(logical(0))
	else{
		# parse line characters into columns and get column names  
		tox_data <- lapply(tox_lines, function(tl) strsplit(tl, "\t", fixed = T)[[1]])
		label_data <- tox_data[[1]]
		# for non-viability assays, only keep the lines of sample activity  
		if(th2 == "ASSAY_OUTCOME"){
			type_id <- which(label_data %in% "SAMPLE_DATA_TYPE")
			data_type <- sapply(tox_data, function(td) td[[type_id]])
			dt_id <- which(data_type %in% "activity")
			tox_data <- tox_data[dt_id]
		}
		else	tox_data <- tox_data[2:length(tox_lines)]
		# obtain assay screening results of all compounds based on the specified result column
		result_id <- which(label_data %in% th2)
		results <- sapply(tox_data, function(td) td[[result_id]])
		# obtain compound PubChem IDs 
		compound_id <- which(label_data %in% "PUBCHEM_CID")
		compounds <- sapply(tox_data, function(td) paste("CID", td[[compound_id]], sep = "_"))
		# obtain compounds with inactive screening results
		inactive_id <- which(results %in% "inactive")
		inactive_compounds <- unique(compounds[inactive_id])
		# obtain compounds with active screening results
		active_id <- grep(th3, results)
		active_compounds <- unique(compounds[active_id])
		# for non-viability assays, remove compound with ambiguous screening results  
		if(th2 == "ASSAY_OUTCOME"){
			inter_compounds <- intersect(inactive_compounds, active_compounds)
			inactive_compounds <- setdiff(inactive_compounds, inter_compounds)
			active_compounds <- setdiff(active_compounds, inter_compounds)
		}
		# for viability assays, count a compound as active as long as it appears active once, otherwise count it as inactive 
		else	inactive_compounds <- setdiff(inactive_compounds, active_compounds)
		# return active/inactive compounds for assays that pass minimum sample threshold  
		if((length(active_compounds) >= min_sample_size) && (length(inactive_compounds) >= min_sample_size)){
			outcome_vec <- c(rep(1, length(active_compounds)), rep(0, length(inactive_compounds)))
			outcome_df <- data.frame(rep(th1, length(outcome_vec)), c(active_compounds, inactive_compounds), outcome_vec)
			colnames(outcome_df) <- c("assay_name", "compound_pubchem_cid", "assay_outcome")
			return(outcome_df)
		}
		else	return(logical(0))
	}
}, tox_files, tox_hyper[1,], tox_hyper[2,], tox_hyper[3,], SIMPLIFY = F);
# output data frame of active/inactive compounds
tf_id <- which(sapply(outcome_list, length) > 0);
outcome_list <- outcome_list[tf_id];
output_df <- do.call(rbind, outcome_list);
write.table(output_df, file = output_file, sep = "\t", row.names = F, col.names = T, quote = F);

