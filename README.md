# CASP15_knotted_artifacts
Programs used in an article:\
*"Knotted artifacts in predicted 3D RNA structures"\
Bartosz A. Gren, Maciej Antczak, Tomasz Zok, Joanna I. Sulkowska, Marta Szachniuk\
PLOS Computational Biology, 2024*

Required packages: topoly, tqdm (not crucial)

# Project structure
`pdb_models/` -- .pdb files of analysed CASP15 predicted models\
`pdb_targets/` -- .pdb files of CASP15 targets\
`pipeline_knots/` -- programs used to study knots in RNA\
`pipeline_lassos_interlaces/` -- programs used to study lassos and interlaces in RNA\
`requirements.txt` -- `pip install -r requirements.txt` to install required packages\

## Content of `pipeline_knots/`
`1_pdb_to_xyz.py` -- program extracting RNA backbones from .pdb files and saving them as .xyz files. Outputs to `xyz_models/` and `xyz_targets`\
`2_check_topols.py` -- main program for determining knot type of .xyz extracted RNAs. Outputs to `res_models.txt` and `res_targets.txt`\
`3_vmd_cyclethrough.py`-- program for manual checking of .pdb models listed in `knots_to_verify.txt`.\
`knots_to_verify.txt` -- list of RNAs with knots after probabilistic closure or TTCs.\

## Content of `pipeline_lassos_interlaces/`


