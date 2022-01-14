##################################################
## Project: Antimicrobial Resistance - Case Fatality Rates for Meningitis Pathogens
## Script purpose: Comnbine data sources for meningitis CFR model 
## Date: 9/4/2020
##################################################

import sys
import os

amrpath = FILEPATH
if not amrpath in sys.path:
    sys.path.append(FILEPATH)

import pandas as pd
import numpy as np
from db_queries import get_covariate_estimates, get_location_metadata
from cod_prep.claude.configurator import Configurator
from cod_prep.utils.misc import print_log_message
from cod_prep.downloaders import (get_current_location_hierarchy,
                                  get_map_version)
from mcod_prep.mcod_mapping import MCoDMapper
from amr_prep.utils.causes import map_to_estimated_pathogens
from mcod_prep.utils.mcause_io import get_mcause_data
from mcod_prep.utils.covariates import get_cov, merge_covariate
from amr_prep.utils.pathogen_formatting import PathogenFormatter
    
CONF = Configurator()
LSV_ID = CONF.get_id("location_set_version")
LS_ID = CONF.get_id("location_set")

hospdat = get_mcause_data("format_map", data_type_id = 3, is_active = True, sub_dirs = "infectious_syndrome", force_rerun=True, block_rerun=False)

meningdat = hospdat.loc[hospdat['infectious_syndrome'].str.contains("mening"),]

meningdat = meningdat.loc[meningdat['infectious_syndrome'].str.contains("unsp|ints_with", regex = True) == False,]

syndrome_to_pathogen = {
        'meningitis-streptococcal_': 'Streptococcus_pneumoniae',
        'meningitis-meningococ': 'Neisseria_meningitidis',
        'meningitis-staphylococcal': 'Staphylococcus_aureus',
        'meningitis-pneumococcal': 'Streptococcus_pneumoniae',
        'meningitis-listeriosis': 'Listeria_monocytogenes',
        'meningitis-hemophilus': 'Haemophilus_influenzae',
        'meningitis-viral': 'Viral'
    }
meningdat = meningdat.assign(pathogen = meningdat['infectious_syndrome'].map(syndrome_to_pathogen))

hospdatsum = meningdat.groupby(['pathogen', "location_id", "age_group_id", "sex_id", "year_id"], as_index = False).agg({
    'admissions':sum,
    'deaths':sum
})

SOURCE = pd.read_csv(FILEPATH)

cnsinf = PathogenFormatter('cns_infectious', 1)
cnsSOURCE = cnsinf.capture_syndrome(SOURCE)

cnsSOURCE2 = cnsSOURCE.loc[(cnsSOURCE['pathogen'].isin(['Staphylococcus aureus', 'Streptococcus pneumoniae', 'Klebsiella pneumoniae', 'Escherichia coli'])) &
                    (~cnsSOURCE['deaths'].isna()),]

SOURCEtopathogen = {
        'Staphylococcus aureus': 'Staphylococcus_aureus',
        'Streptococcus pneumoniae': 'Streptococcus_pneumoniae',
        'Klebsiella pneumoniae': 'Klebsiella_pneumoniae',
        'Escherichia coli': 'Escherichia_coli'
    }

cnsSOURCE2 = cnsSOURCE2.assign(pathogen = cnsSOURCE2['pathogen'].map(SOURCEtopathogen))

cnsSOURCEsum = cnsSOURCE2.groupby(['pathogen', "location_id", "age_group_id", "sex_id", "year_id"], as_index = False).agg({
    'admissions':sum,
    'deaths':sum
})

splitlit = pd.read_csv(FILEPATH)

hospdatsum['data_type'] = "hospital"
cnsSOURCEsum['data_type'] = "SOURCE"
allcfrdat = hospdatsum.append([cnsSOURCEsum, splitlit], ignore_index = True)

def agegrouper(row):
    if np.isin(row['age_group_id'], [2, 3, 4, 5]):
        return 'Neonatal_5'
    elif np.isin(row['age_group_id'], [range(6,18)]):
        return '5_65'
    elif row['age_group_id'] >= 18:
        return '65_plus' 
    else:
        return 'other'

def yeargrouper(row):
    if np.isin(row['year_id'], [range(1979,2000)]):
        return '<2000'
    elif np.isin(row['year_id'], [range(2000,2020)]):
        return '>2000'  
    else:
        return 'other'
    
allcfrdat['agecat'] = allcfrdat.apply(agegrouper, axis=1)
allcfrdat['yrcat'] = allcfrdat.apply(yeargrouper, axis=1)

allcfrdat.loc[allcfrdat['nid'].isna(), 'study_id'] = allcfrdat.loc[allcfrdat['nid'].isna(), 'data_type']
allcfrdat.loc[~allcfrdat['nid'].isna(), 'study_id'] = allcfrdat.loc[~allcfrdat['nid'].isna(), 'nid']

allcfrdat.loc[(allcfrdat['year_id'] == 1979) & (allcfrdat['location_id'] == 93), 'year_id'] = 1980
allcfrdat.loc[(allcfrdat['year_id'] == 1979) & (allcfrdat['location_id'] == 89), 'year_id'] = 1980

locs = get_location_metadata(location_set_id = LS_ID)
vstlnd_2_norway = np.int(locs.loc[locs['location_id'] == 60132, 'parent_id'])
allcfrdat.loc[allcfrdat['location_id'] == 60132, 'location_id'] = vstlnd_2_norway

allcfrdat = merge_covariate(allcfrdat, "haqi")

allcfrdatsum = allcfrdat.groupby(['pathogen', 'location_id', 'agecat', 'study_id', 'year_id', 'yrcat'], as_index = False).agg({
    'admissions':sum,
    'deaths':sum,
    'haqi':'mean'
})

allcfrdatsum.to_csv(FILEPATH, index = False)

