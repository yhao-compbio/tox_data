## usr/bin/Rscript 
## created by Yun Hao @MooreLab 2019
## This script maps small molecules to their protein targets 

# Input arguments
Args		<- commandArgs(T);
compound_file	<- Args[1];	## File that contains names of compounds, one name per line.

# read in compound~target mapping from DrugBank
drugbank_df <- read.csv(file = "downloads/uniprot\ links.csv");
drugbank_target <- as.character(drugbank_df$UniProt.ID);
drugbank_name <- as.character(tolower(drugbank_df$Name));

# read in compound file
compounds <- readLines(compound_file);
compounds <- tolower(compounds);

# map compounds to protein targets
componud_target_list <- lapply(compounds, function(cpd){
	cpd_id <- which(drugbank_name %in% cpd);
	return(drugbank_target[cpd_id]);
});
names(componud_target_list) <- compounds

# create binary matrix of targets(columns) ~ compounds(rows) based on the mapping results.  
all_targets <- unique(unlist(componud_target_list));
N_targets <- length(all_targets);
compound_target_df <- mapply(function(ctl){
	tar_vec <- rep(NA, N_targets);
	names(tar_vec) <- all_targets;
	if(length(ctl) > 0){
		tar_vec[1:N_targets] <- 0
		tar_vec[ctl] <- 1
	}
	return(tar_vec);
}, componud_target_list);
compound_target_df <- t(compound_target_df);

# write mapping matrix to output file
out_file <- paste(compound_file, "_target.tsv", sep = "");
write.table(compound_target_df, file = out_file, row.names = F, col.names = T, quote = F, sep = "\t");
