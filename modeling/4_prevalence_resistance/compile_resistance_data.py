import pandas as pd
import numpy as np
import sys
import os
import argparse
import getpass
user = getpass.getuser()
amrpath = FILEPATH
if not amrpath in sys.path:
    sys.path.append(FILEPATH)
from amr_prep.utils.amr_io import get_amr_data
from mcod_prep.utils.mcause_io import get_mcause_data
from datetime import datetime

amr_dat = get_amr_data()

resdat = amr_dat.loc[(amr_dat['resistance'].notna()) & (amr_dat['resistance'].isin(['sensitive', 'resistant', 'susceptible'])), :]
resdat['resistance'] = resdat['resistance'].str.replace('sensitive', 'susceptible')

ressum = resdat.groupby(['location_id', 'year_id', 'source', 'pathogen', 'abx_class'], group_keys=False, as_index = False).apply(
    lambda x: pd.Series({'resistant':sum(x.loc[x['resistance'] == 'resistant', 'cases']),
                        'susceptible':sum(x.loc[x['resistance'] == 'susceptible', 'cases']),
                        'cases':sum(x['cases'])})).reset_index()

needrespath = ['acinetobacter_baumanii',
              'campylobacter',
              'chlamydia_spp',
              'citrobacter_spp',
              'clostridium_difficile',
              'cryptosporidiosis',
              'enterobacter_spp',
              'enterococcus_faecium',
              'enterococcus_faecalis',
              'enterococcus_spp',
              'escherichia_coli',
              'group_a_strep',
              'group_b_strep',
              'haemophilus_influenzae',
              'klebsiella_pneumoniae',
              'listeria',
              'moraxella_spp',
              'mycoplasma',
              'neisseria_meningitidis',
              'proteus_spp',
              'providencia_spp',
              'pseudomonas_aeruginosa',
              'pseudomonas_spp',
              'serratia_spp',
              'legionella_spp',
              'salmonella_enterica',
              'streptococcus_pneumoniae',
              'staphylococcus_aureus']

ressum2 = ressum.loc[ressum['pathogen'].isin(needrespath),:].reset_index()
ressum2.to_csv(FILEPATH, index = False)

fqgono = resdat.loc[(resdat['pathogen'] == 'neisseria_gonorrheae') & (resdat['abx_class'] == 'fluoroquinolone'),:]

fqgonosum = fqgono.groupby(['location_id', 'year_id', 'source', 'pathogen', 'abx_class'], group_keys=False, as_index = False).apply(
    lambda x: pd.Series({'resistant':sum(x.loc[x['resistance'] == 'resistant', 'cases']),
                        'susceptible':sum(x.loc[x['resistance'] == 'susceptible', 'cases']),
                        'cases':sum(x['cases'])})).reset_index()

ecdc = pd.read_csv(FILEPATH)
ecdc['pathogen'] = 'neisseria_gonorrheae'
ecdc = ecdc.rename(columns = {'nid':'source', 'n_resistant':'resistant', 'sample_size':'cases', 'antibiotic':'abx_class'})
ecdc['susceptible'] = ecdc['cases'] - ecdc['resistant']
ecdc2 = ecdc.loc[:,['location_id', 'year_id', 'source', 'pathogen', 'abx_class', 'resistant', 'susceptible', 'cases']]
fqgonosum = fqgonosum.append(ecdc2)

oxlit = pd.read_csv(FILEPATH)
oxlit = oxlit.rename(columns = {'nid':'source', 'n_resistant':'resistant', 'sample_size':'cases'})
oxlit['antimicrobial'] = 'fluoroquinolone'
oxlit['pathogen'] = 'neisseria_gonorrheae'
oxlit = oxlit.rename(columns = {'nid':'source', 'n_resistant':'resistant', 'sample_size':'cases', 'antimicrobial':'abx_class'})
oxlit['susceptible'] = oxlit['cases'] - oxlit['resistant']
oxlit2 = oxlit.loc[:,['location_id', 'year_id', 'source', 'pathogen', 'abx_class', 'resistant', 'susceptible', 'cases']]

oxlit2 = oxlit2.loc[oxlit2['location_id'].notna(),:]

fqgonosum = fqgonosum.append(oxlit2)

fqgonosum.to_csv(FILEPATH, index = False)