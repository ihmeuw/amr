import pandas as pd
import numpy as np
import sys
import os
import argparse
import getpass
user = getpass.getuser()
simbinpath = FILEPATH
if not simbinpath in sys.path:
    sys.path.append(FILEPATH)
amrpath = FILEPATH
if not amrpath in sys.path:
    sys.path.append(FILEPATH)
from amr_prep.utils.amr_io import get_amr_data
from mcod_prep.utils.mcause_io import get_mcause_data
from simbin.random import sample_simplex, sample_bin_prob
from simbin.optim import adjust_prob_mat
import itertools

IN_DIR = FILEPATH

corrs = pd.read_csv(FILEPATH)
corrs_salmonella = pd.read_csv(FILEPATH)
corrs = corrs.append(corrs_salmonella)

studycombos = pd.read_csv(FILEPATH)

PATHOGEN = str(sys.argv[1])
LOCATION_ID = str(sys.argv[2])

resdf = pd.DataFrame()
stgprres = pd.DataFrame()
bugdf = corrs.loc[corrs['pathogen'] == PATHOGEN, :]
mtx = bugdf.pivot_table(index='abx_a', columns='abx_b', values='pearson')
bugdrugs = bugdf['abx_a'].sort_values().unique()

for drug in bugdrugs:
    bugfile = pd.read_csv(
        (
            FILEPATH
            + studycombos.loc[
                (studycombos['pathogen'] == PATHOGEN) & (studycombos['abx_class'] == drug),
                'resistance_run_id'
            ].astype('str')
            + "/draws_temp_0/" + str(LOCATION_ID) + ".csv"
        ).values[0]
    )
    bugfile['antibiotic'] = drug
    stgprres = stgprres.append(bugfile.loc[bugfile['year_id'] == 2018, :])

for draw in ["draw_" + str(i) for i in range(0, 1000)]:
    for drug in bugdrugs:
        mtx.loc[drug, drug] = stgprres.loc[stgprres['antibiotic'] == drug, draw].values[0]
    for drug1 in bugdrugs:
        for drug2 in [drug2 for drug2 in bugdrugs if drug2 not in drug1]:
            pearson = bugdf.loc[(bugdf['abx_a'] == drug1) & (
                bugdf['abx_b'] == drug2), 'pearson'].values
            mtx.loc[drug1, drug2] = pearson * np.sqrt(
                mtx.loc[drug1, drug1] * (1 - mtx.loc[drug1, drug1])
                * mtx.loc[drug2, drug2] * (1 - mtx.loc[drug2, drug2])
            ) + mtx.loc[drug1, drug1] * mtx.loc[drug2, drug2]
    probab_mat = mtx.values
    adjusted_prob_mat, prob = adjust_prob_mat(
        prob_mat=probab_mat, x0=sample_simplex(2**probab_mat.shape[0], size=1)[0])
    prob = prob.reshape(np.repeat([2], len(bugdrugs)))

    lst = list(itertools.product([0, 1], repeat=len(bugdrugs)))
    combos = pd.DataFrame.from_records(lst)
    combos.columns = bugdrugs
    for i in range(0, len(combos)):
        combos.loc[i, 'combinatoric_prop'] = prob[tuple(
            combos.loc[i, bugdrugs].astype('int').values)]
    combos.reset_index(inplace=True)
    combos['index'] = PATHOGEN + '-' + combos['index'].astype('str')
    combos = combos.rename(columns={'index': 'combinatoric_id'})
    melted = combos.melt(id_vars=['combinatoric_id', 'combinatoric_prop'],
                         var_name='abx_class', value_name='resistant')
    melted['draw'] = draw
    melted['location_id'] = LOCATION_ID
    melted['pathogen'] = PATHOGEN
    resdf = resdf.append(melted)
resdf2 = resdf.pivot_table(index=['combinatoric_id', 'pathogen', 'location_id',
                                  'abx_class', 'resistant'],
                           columns='draw', values='combinatoric_prop').reset_index()
outdir = FILEPATH +  str(PATHOGEN) + '/'
if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)
resdf2.to_csv(outdir + str(LOCATION_ID) + '.csv', index=False)
