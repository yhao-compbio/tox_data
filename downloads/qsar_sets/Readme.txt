Notes on electronic files

****************************************************************************************

The ACE, GPB, THER and THR sets were not compiled by us.  Please acknowledge the original
authors by citing their work in addition to this article.

ACE:

Depriest, S. A.; Mayer, D.; Naylor, C. B.; Marshall, G. R. 3D-QSAR of angiotensin-converting 
enzyme and thermolysin inhibitors - a comparison of CoMFA models based on deduced and 
experimentally determined active-site geometries. J. Am. Chem. Soc. 1993, 115, 5372-5384.

GPB and THER:

Gohlke, H.; Klebe, G. Drugscore meets CoMFA: Adaptation of fields for molecular comparison 
(AFMoC) or how to tailor knowledge-based pair-potentials to a particular protein. J. Med. Chem. 
2002, 45, 4153-4170.

THR:

Bohm, M.; Sturzebecher, J.; Klebe, G. Three-dimensional quantitative structure-activity 
relationship analyses using comparative molecular field analysis and comparative molecular 
similarity indices analysis to elucidate selectivity differences of inhibitors binding to trypsin, 
thrombin, and factor Xa. J. Med. Chem. 1999, 42, 458-477.

****************************************************************************************

Files in MDL SD and Sybyl mol2 formats are provided.  The files are tarred and gzipped 
contents of the directories listed below (see UNIX commands tar and gzip), stored within a 
Windows WinZip archive (only ZIP is supported by the ACS; use unzip for opening under 
Linux).  When .tar.gz files are opened in WinZip, the directory structure is not visible;  use 
"extract ..." in the action menu rather than drag-and-drop for preserving the directory structure.

3dqsar_mol2: Sybyl mol2 databases of aligned compounds, and tab-separated text files 
containing spreadsheet data.

3dqsar_sd: MDL SD files of aligned compounds (the format does not allow specification of 
charges; use the mol2 databases if you want our calculated charges).

corina_sd: MDL SD files of CORINA-generated structures.

25d_descriptors: contains tab-separated 2.5D descriptors for each set.

File contents:
- Min_dist indicates the Tanimoto coefficient calculated between each test set compound 
and the most similar training set compound; Avg_dist is the average value calculated 
over all training set compounds.  These may be viewed as measures for the degree of 
extrapolation from the training set.
- In the set column, 1, 2 and 3 indicate training, test and inactive compounds, respectively.
In the SD file of CORINA structures for the DHFR set, activities in uM are given for P. carinii 
(PC) and T. gondii (TG) in addition to rat liver (RL); only RL activities are used in this work.


