# Process toxicological resources and datasets

# Source datasets

+ [`downloads/ryan_reference_set/`](downloads/ryan_reference_set/) contains four curated datasets from [Defining a Reference Set to Support Methodological Research in Drug Safety](https://doi.org/10.1007/s40264-013-0097-8). Each dataset contains drugs that cause or do not cause a particular side effect.  

+ [`downloads/qsar_sets/`](downloads/qsar_sets/) contains QSAR benchmark datasets from [A Comparison of Methods for Modeling Quantitative Structure−Activity Relationships](https://doi.org/10.1021/jm0497141)

# Generated datasets 

+ [`data/aquatic_atompair_dataset.tsv`](data/aquatic_atompair_dataset.tsv), [`data/aquatic_daylight_fp_dataset.tsv`](data/aquatic_daylight_fp_dataset.tsv), [`data/aquatic_dragon_dataset.tsv`](data/aquatic_dragon_dataset.tsv), [`data/aquatic_lcalc_dataset.tsv`](data/aquatic_lcalc_dataset.tsv), [`data/aquatic_moe2d_dataset.tsv`](data/aquatic_moe2d_dataset.tsv), [`data/aquatic_moe2d_fp_dataset.tsv`](data/aquatic_moe2d_fp_dataset.tsv), [`data/aquatic_moe3d_dataset.tsv`](data/aquatic_moe3d_dataset.tsv), [`data/aquatic_pipelinepilot_fp_dataset.tsv`](data/aquatic_pipelinepilot_fp_dataset.tsv), [`data/aquatic_quickprop_dataset.tsv`](data/aquatic_quickprop_dataset.tsv) are nine datasets generated from "AqaticTox" in "QSARdata" R package.

+ [`data/blood_brain_barrier_atompair_dataset.tsv`](data/blood_brain_barrier_atompair_dataset.tsv), [`data/blood_brain_barrier_daylight_fp_dataset.tsv`](data/blood_brain_barrier_daylight_fp_dataset.tsv), [`data/blood_brain_barrier_dragon_dataset.tsv`](data/blood_brain_barrier_dragon_dataset.tsv), [`data/blood_brain_barrier_lcalc_dataset.tsv`](data/blood_brain_barrier_lcalc_dataset.tsv), [`data/blood_brain_barrier_moe2d_dataset.tsv`](data/blood_brain_barrier_moe2d_dataset.tsv), [`data/blood_brain_barrier_moe2d_fp_dataset.tsv`](data/blood_brain_barrier_moe2d_fp_dataset.tsv), [`data/blood_brain_barrier_moe3d_dataset.tsv`](data/blood_brain_barrier_moe3d_dataset.tsv), [`data/blood_brain_barrier_pipelinepilot_fp_dataset.tsv`](data/blood_brain_barrier_pipelinepilot_fp_dataset.tsv), [`data/blood_brain_barrier_quickprop_dataset.tsv`](data/blood_brain_barrier_quickprop_dataset.tsv) are nine datasets generated from "bbb2"in "QSARdata" R package.

+ [`data/caco2_permeability_atompair_dataset.tsv`](data/caco2_permeability_atompair_dataset.tsv), [`data/caco2_permeability_dragon_dataset.tsv`](data/caco2_permeability_dragon_dataset.tsv), [`data/caco2_permeability_pipelinepilot_fp_dataset.tsv`](data/caco2_permeability_pipelinepilot_fp_dataset.tsv), [`data/caco2_permeability_quickprop_dataset.tsv`](data/caco2_permeability_quickprop_dataset.tsv) are four datasets generated from "caco" in "QSARdata" R package. 

+ [`data/drug_induced_pld_atompair_dataset.tsv`](data/drug_induced_pld_atompair_dataset.tsv), [`data/drug_induced_pld_dragon_dataset.tsv`](data/drug_induced_pld_dragon_dataset.tsv), [`data/drug_induced_pld_lcalc_dataset.tsv`](data/drug_induced_pld_lcalc_dataset.tsv), [`data/drug_induced_pld_pipelinepilot_fp_dataset.tsv`](data/drug_induced_pld_pipelinepilot_fp_dataset.tsv), [`data/drug_induced_pld_quickprop_dataset.tsv`](data/drug_induced_pld_quickprop_dataset.tsv), [`data/drug_induced_pld_volsurfplus_dataset.tsv`](data/drug_induced_pld_volsurfplus_dataset.tsv) are six datasets generated from "PLD" in "QSARdata" R package.

+ [`data/melting_point_moe2d3d_dataset.tsv`](data/melting_point_moe2d3d_dataset.tsv) is a dataset generated from "MeltingPoint" in "QSARdata" R package.

+ [`data/mutagen_dragon_dataset.tsv`](data/mutagen_dragon_dataset.tsv) is a dataset generated from "Mutagen" in "QSARdata" R package.

# Source codes
 
+ [`src/r_qsardata.R`](src/r_qsardata.R) contains the codes that extract benchmark datasets from QSARdata R package. 
