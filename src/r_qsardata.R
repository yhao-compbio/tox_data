## usr/bin/Rscript 
## created by Yun Hao @MooreLab 2019
## This script extracts benchmark datasets from QSARdata R package

## Packages
# install.packages("QSARdata");
library(QSARdata);

## function 
remove.duplicate.molecules <- function(df){
	uni_molec <- as.character(df$Molecule);
	uni_molec_id <- lapply(uni_molec, function(um){
		um_id <- which(uni_molec %in% um);
		if(length(um_id) == 1)	return(um_id)
		else	return(logical(0))
	})
	uni_molec_id <- unlist(uni_molec_id);
	uni_df <- df[uni_molec_id,];
	return(uni_df);
};

combine.feature.response <- function(df1, df2){
	uni_df1 <- remove.duplicate.molecules(df1);
	uni_df2 <- remove.duplicate.molecules(df2);
	combine_df <- merge(uni_df1, uni_df2, by = "Molecule");
	return(combine_df);
};

## AquaticTox dataset
data(AquaticTox);
df_list <- list(AquaticTox_AtomPair, AquaticTox_Daylight_FP, AquaticTox_Dragon, AquaticTox_Lcalc, AquaticTox_moe2D, AquaticTox_moe2D_FP, AquaticTox_moe3D, AquaticTox_PipelinePilot_FP, AquaticTox_QuickProp);
names(df_list) <- c("atompair","daylight_fp","dragon","lcalc","moe2d","moe2d_fp","moe3d","pipelinepilot_fp","quickprop");
output_combined_dataset <- mapply(function(dl, ndl){
	data_df <- combine.feature.response(dl, AquaticTox_Outcome);
	writeLines(as.character(data_df$Molecule), "data/aquatic_dataset_molecules.txt");
	data_df <- data_df[-1];
	write.table(data_df, file = paste("data/aquatic_",ndl,"_dataset.tsv",sep=""), row.names = F, col.names = T, sep = "\t", quote = F);
	return(1);
}, df_list, names(df_list));

## Blood-Brain Barrier dataset
data(bbb2);
df_list <- list(bbb2_AtomPair, bbb2_Daylight_FP, bbb2_Dragon, bbb2_Lcalc, bbb2_moe2D, bbb2_moe2D_FP, bbb2_moe3D, bbb2_PipelinePilot_FP, bbb2_QuickProp);
names(df_list) <- c("atompair","daylight_fp","dragon","lcalc","moe2d","moe2d_fp","moe3d","pipelinepilot_fp","quickprop");
output_combined_dataset <- mapply(function(dl, ndl){
	data_df <- combine.feature.response(dl, bbb2_Outcome);
	writeLines(as.character(data_df$Molecule), "data/blood_brain_barrier_dataset_molecules.txt");
	data_df <- data_df[-1];
	write.table(data_df, file = paste("data/blood_brain_barrier_",ndl,"_dataset.tsv",sep=""), row.names = F, col.names = T, sep = "\t", quote = F);
	return(1);
}, df_list, names(df_list));

## Caco-2 Permeability dataset
data(caco);
df_list <- list(caco_AtomPair, caco_Dragon, caco_PipelinePilot_FP, caco_QuickProp);
names(df_list) <- c("atompair","dragon","pipelinepilot_fp","quickprop");
output_combined_dataset <- mapply(function(dl, ndl){
	data_df <- combine.feature.response(dl, caco_Outcome);
	writeLines(as.character(data_df$Molecule), "data/caco2_permeability_dataset_molecules.txt");
	data_df <- data_df[-1];
	write.table(data_df, file = paste("data/caco2_permeability_",ndl,"_dataset.tsv",sep=""), row.names = F, col.names = T, sep = "\t", quote = F);
	return(1);
}, df_list, names(df_list));

## Melting Point dataset
data(MeltingPoint);
data_df <- data.frame(MP_Descriptors, MP_Outcome);
write.table(data_df, file = paste("data/melting_point_moe2d3d_dataset.tsv",sep=""), row.names = F, col.names = T, sep = "\t", quote = F);

## Mutagen dataset
data(Mutagen);
data_df <- data.frame(Mutagen_Dragon, Mutagen_Outcome);
write.table(data_df, file = paste("data/mutagen_dragon_dataset.tsv",sep=""), row.names = F, col.names = T, sep = "\t", quote = F);

## Drug-Induced Phospholipidosis dataset
data(PLD);
df_list <- list(PLD_AtomPair, PLD_Dragon, PLD_PipelinePilot_FP, PLD_QuickProp, PLD_VolSurfPlus, PLD_LCALC);
names(df_list) <- c("atompair","dragon","pipelinepilot_fp","quickprop","volsurfplus","lcalc");
output_combined_dataset <- mapply(function(dl, ndl){
	data_df <- combine.feature.response(dl, PLD_Outcome);
	writeLines(as.character(data_df$Molecule), "data/drug_induced_pld_dataset_molecules.txt");
	data_df <- data_df[-1];
	write.table(data_df, file = paste("data/drug_induced_pld_",ndl,"_dataset.tsv",sep=""), row.names = F, col.names = T, sep = "\t", quote = F);
	return(1);
}, df_list, names(df_list));


