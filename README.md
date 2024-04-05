# CASP15_knotted_artifacts
Programs used in "Knotted artifacts in predicted 3D RNA structures"\
Bartosz A. Gren, Maciej Antcza, Tomasz Zok, Joanna I. Sulkowska, Marta Szachniuk

Required packages: topoly, tqdm, math, os, signal, sys

# Project structure
`pdb_models/` -- .pdb files of analysed CASP15 predicted models\
`pdb_targets/` -- .pdb files of CASP15 targets\
`pipeline_knots/` -- programs used to study knots in RNA\
`pipeline_lassos_interlaces/` -- programs used to study lassos and interlaces in RNA\

## pipeline_knots/ contents
`1_pdb_to_xyz.py` -- program extracting RNA backbones from .pdb files and saving them as .xyz files\
`2_check_topols.py` -- main program for determining knot type of .xyz extracted RNAs\
`3_vmd_cyclethrough.py`-- program for manual checking of .pdb models listed in `knots_to_verify.txt`\
`knots_to_verify.txt` -- list of RNAs with knots after probabilistic closure or TTCs\
`xyz_models/` -- .xyz files of extracted RNA backbones of analysed CASP15 predicted models\
`xyz_targets/` -- .xyz files of extracted RNA backbones of CASP15 targets\


