import numpy as np
import scipy as sp
import uproot as up
import pandas as pd
import awkward as ak
import awkward_pandas as akpd
from IPython.display import display


tr_ph_filename = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/mcgpj/tr_ph v9/EnergySmearing/MCGPJ_kskl501_Merged.root"

with up.open(f"{tr_ph_filename}:tr_ph_merged") as ksTree: # type: ignore
    tree: pd.DataFrame = ksTree.arrays(['emeas', 'runnum', 'nt', 'tnhit', 'tdedx', 'tptot', 'is_coll', 'tth', 'tphi', 'nph', 'phen', 'phth', 'phphi'], library="pd")

tree = tree.query('0 < nt < 3')   
tree = tree.query('10 < tnhit[0]')   
# tree = tree.query('emeas > 500')   

display(tree.head())