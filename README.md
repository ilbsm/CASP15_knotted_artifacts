# CASP15_knotted_artifacts
Programs used in an article:\
*"Knotted artifacts in predicted 3D RNA structures"\
Bartosz A. Gren, Maciej Antczak, Tomasz Zok, Joanna I. Sulkowska, Marta Szachniuk\
PLOS Computational Biology, 2024*

Required packages: topoly, tqdm\
`pip install -r requirements.txt` to install them.

# Project structure
`pdb_models/` -- .pdb files of analysed CASP15 predicted models.\
`pdb_targets/` -- .pdb files of CASP15 targets.\
`pipeline_knots/` -- programs used to study knots in RNA.\
`pipeline_lassos_interlaces/` -- programs used to study lassos and interlaces in RNA.\
`secondary_structure/` -- script and input files to generate pictures of secondary structures of models from figures 5 and 6.\
`requirements.txt` -- required packages.

## Content of `pipeline_knots/`
`1_pdb_to_xyz.py` -- program which extracts RNA backbones from .pdb files and saves them as .xyz files. It outputs to `xyz_models/` and `xyz_targets`.\
`2_check_topols.py` -- main program for determining knot type of .xyz extracted RNAs. It outputs to `res_models.txt` and `res_targets.txt`.\
`3_vmd_cyclethrough.py`-- program for manual checking of .pdb models listed in `knots_to_verify.txt`.\
`knots_to_verify.txt` -- list of RNAs with knots after probabilistic closure or TTCs.

## Content of `pipeline_lassos_interlaces/`

Pre-processing of 3D models:\
`clean-casp-headers.awk` -- remove CASP-specific headers.

Identification and classification of entanglements of structural elements of 3D RNA models:\
`R*/*_report.csv` -- entanglement identification and classification results obtained using RNAspider.\
`casp15-entanglements-summary.xlsx` -- a comprehensive analysis of entanglements identified.

Lasso depth analysis:\
`1_normalize_csv.py` -- a script normalizing CSV files.\
`2_analyze_depth.py` -- a script preparing lasso depth analysis summary.\
`lasso-depth-analysis.xlsx` -- the resultant summary.
