# Collect and process toxicological resources and datasets

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

# References

+ Ryan, Patrick B., et al. "Defining a reference set to support methodological research in drug safety." Drug safety 36.1 (2013): 33-47.

+ DePriest, Scott A., et al. "3D-QSAR of angiotensin-converting enzyme and thermolysin inhibitors: a comparison of CoMFA models based on deduced and experimentally determined active site geometries." Journal of the American Chemical Society 115.13 (1993): 5372-5384.

+ Gohlke, Holger, and Gerhard Klebe. "DrugScore meets CoMFA: adaptation of fields for molecular comparison (AFMoC) or how to tailor knowledge-based pair-potentials to a particular protein." Journal of medicinal chemistry 45.19 (2002): 4153-4170.

+ Böhm, Markus, Jörg Stürzebecher, and Gerhard Klebe. "Three-dimensional quantitative structure− activity relationship analyses using comparative molecular field analysis and comparative molecular similarity indices analysis to elucidate selectivity differences of inhibitors binding to trypsin, thrombin, and factor Xa." Journal of medicinal chemistry 42.3 (1999): 458-477.

+ Sutherland, Jeffrey J., Lee A. O'Brien, and Donald F. Weaver. "A comparison of methods for modeling quantitative structure− activity relationships." Journal of Medicinal Chemistry 47.22 (2004): 5541-5554.

+ Burns et al. 'A mathematical model for prediction of drug molecule diffusion across the blood-brain barrier'. The Canadian Journal of Neurological Sciences (2004) vol. 31 (4) pp. 520-527

+ He and Jurs. 'Assessing the reliability of a QSAR model's predictions'. Journal of Molecular Graphics and Modelling (2005) vol. 23 (6) pp. 503-523

+ Karthikeyan et al. 'General melting point prediction based on a diverse compound data set and artificial neural networks'. Journal of chemical information and modeling (2005) vol. 45 (3) pp. 581-90

+ Kazius et al. 'Derivation and validation of toxicophores for mutagenicity prediction'. Journal of medicinal chemistry(Print) (2005) vol. 48 (1) pp. 312-320 
