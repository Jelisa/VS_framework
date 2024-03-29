atomic_table_weights = {
    "h": 1,
    "c": 12,
    "n": 14,
    "o": 16,
    "f": 19,
    "na  ": 23,
    "mg": 12,
    "p": 31,
    "s": 32,
    "cl  ": 35,
    "k": 39,
    "ca  ": 20,
    "zn": 65,
    "br": 80,
    "i": 127,
    "b": 11
}

accepted_formats = "\.mae|\.pdb"  # This specifies the patterns for the accepted formats.

templates_keywords_pattern = "\$\{*(\w*_*w*)\}?"

complexes_pre_defined_name = ".*_complex_processed.pdb"
complexes_pre_created_name = "*_processed.pdb"

datalocal_templates = "DataLocal/Templates/OPLS2005/HeteroAtoms/"
datalocal_rotamerlib = "DataLocal/LigandRotamerLibs/"
datalocal_obc_param = "DataLocal/OBC/"

