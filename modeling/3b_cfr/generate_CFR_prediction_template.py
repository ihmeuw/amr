"""
Generate Prediction Template for CFR Models from Config
"""
import sys
import os
import getpass
user = getpass.getuser()
amrpath = FILEPATH
if not amrpath in sys.path:
    sys.path.append(FILEPATH)
    
import pandas as pd
import subprocess
import numpy as np
from pathlib import Path
from amr_prep.utils.pathogen_formatting import PathogenFormatter
from amr_prep.utils.data_aggregation import prep_lit_data
from cod_prep.utils import report_duplicates, report_if_merge_fail, create_square_df
from db_queries import get_covariate_estimates
import itertools
from cod_prep.claude.claude_io import Configurator
from mcod_prep.utils.covariates import merge_covariate
from mcod_prep.utils.nids import add_nid_metadata
from cod_prep.downloaders import (
    add_age_metadata, get_current_location_hierarchy,
    add_location_metadata, get_cod_ages, pretty_print,
)
from mcod_prep.utils.mcause_io import get_mcause_data
from db_tools import ezfuncs

CONF = Configurator()

df = ezfuncs.query("SELECT * FROM shared.covariate", conn_def=CONF.get_database_setup('db'))
gbdcovs = list(df['covariate_name_short'])

inf_synd = SYNDROME
model_id = MODEL
outpath = FILEPATH

os.makedirs(outpath, exist_ok = True)

def year_categorizer(df, yearcol, yrcats):
    yrcats.insert(0,np.min(df[yearcol]))
    yrcats.append(np.max(df[yearcol]+1))
    yrcats = pd.Series(yrcats).drop_duplicates().tolist()
    df['yrcat'] =  pd.cut(df[yearcol], yrcats, right = False).astype(str)
    df['yrcat'] = df['yrcat'].str.replace("[", "")
    df['yrcat'] = df['yrcat'].str.replace(")", "")
    df['yrcat'] = df['yrcat'].str.replace(", ", "_")
    df['yrcat'] = df['yrcat'].str.replace(str(np.max(df[yearcol]+1)), str(np.max(df[yearcol])))
    return(df)

def create_predictions_template(pathogens, ages, years, sexes, covar_df):
    lh = get_current_location_hierarchy(
        location_set_version_id=CONF.get_id('location_set_version'))\
        .query("level == 3")
    locs = lh.location_id.unique().tolist()
    index = pd.MultiIndex.from_product(
        [locs, ages, pathogens, years, sexes],
        names=['location_id', 'age_group_id', 'pathogen', 'year_id', 'sex_id']
    )
    square_df = pd.DataFrame(index=index).reset_index()
    for cov in covar_df['covariate']:
        if (covar_df.loc[covar_df['covariate'] == cov, 'covariate_type'] == 'continuous').bool():
            if(cov in gbdcovs):
                square_df = merge_covariate(
                    square_df, cov, decomp_step=CONF.get_id("decomp_step"),
                    gbd_round_id=CONF.get_id("gbd_round"))
                cov_merge_cols = str(covar_df.loc[covar_df['covariate'] == cov, 'pred_merge_cols'].values[0]).split(", ")
                report_if_merge_fail(square_df, cov, cov_merge_cols)
    return square_df

globparams = pd.read_excel(FILEPATH, sheet_name = 0)
covs = pd.read_excel(FILEPATH, sheet_name = 1)
globparams = globparams.loc[(globparams['infectious_syndrome'] == inf_synd) & (globparams['model_id'] == model_id), :]
globparams = globparams.reset_index(drop = True)
covs = covs.loc[(covs['infectious_syndrome'] == inf_synd) & (covs['model_id'] == model_id), :]

if globparams['eval_pathogens'].isna().item():
    pred_path = str(globparams['model_pathogens'].values[0]).split(", ")
else:
    pred_path = str(globparams['eval_pathogens'].values[0]).split(", ")
years = globparams['years'].str.split(", ")[0]
pred_years = range(int(years[0]), int(years[1])+1)
if globparams['age_subset'].isna()[0]:
    pred_ages = str(globparams['age_categories'].values[0]).split(", ")
else:
    pred_ages = str(globparams['age_subset'].values[0]).split(", ")
if np.isnan(globparams['sex_val'].values[0]):
    pred_sexes = [3]
else:
    pred_sexes = str(globparams['sex_val'].values[0]).split(", ")
    
categories = set(covs.loc[covs['covariate_type'] == 'category', 'covariate'].tolist()) - set(['age_group_id', 'pathogen'])

pred_template = create_predictions_template(pred_path, pred_ages, pred_years, pred_sexes, covs)

pred_template['source'] = 'US_MedMined_BD'

if 'ICU' in categories:
    pred_template['ICU'] = 'mixed'
    
if 'hosp' in categories:
    pred_template_hosp = pred_template.copy()
    pred_template['hosp'] = 'community'
    pred_template_hosp['hosp'] = 'hospital'
    pred_template = pred_template.append(pred_template_hosp)

pred_template.to_csv(path_or_buf = outpath, index = False)