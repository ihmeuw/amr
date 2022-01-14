import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
from db_queries import get_outputs
from cod_prep.downloaders import get_current_location_hierarchy
from cod_prep.downloaders.ages import get_cod_ages, getcache_age_aggregate_to_detail_map
from cod_prep.claude.claude_io import Configurator
from db_queries import get_location_metadata
from mcod_prep.utils.covariates import merge_covariate

CONF = Configurator()

def agg_ages(df, agg_ages, outcol="agg_age_group_id"):
    age_agg_map = (
        getcache_age_aggregate_to_detail_map(force_rerun=False, block_rerun=True)
        .loc[lambda d: d["agg_age_group_id"].isin(agg_ages), :]
        .set_index("age_group_id", verify_integrity=True)['agg_age_group_id']
    )
    df[outcol] = df["age_group_id"].map(age_agg_map)
    return df

# get LRI incidence by age/country, to weight STREP VE data within detail age groups to aggregate age groups
incidence = get_outputs('cause', location_id = 'lvl3', cause_id = 322, year_id = 2019, compare_version_id = 7244, gbd_round_id = 6,
           decomp_step = 'step5', measure_id = 6, metric_id = 1, age_group_id = [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,30,31,32,235])
incidence = incidence.loc[incidence['sex_id'] == 3,['age_group_id', 'location_id', 'val']]
incidence['val_sum'] = incidence.groupby(['age_group_id', 'location_id'])['val'].transform('sum')
incidence = agg_ages(incidence, [42,179,239,25,26])
incidence['val_sum'] = incidence.groupby(['agg_age_group_id', 'location_id'])['val'].transform('sum')
incidence['incidence_weight'] = incidence['val']/incidence['val_sum']
iweights = incidence.loc[:,['age_group_id', 'agg_age_group_id', 'location_id', 'incidence_weight']]

# loop through files, apply incidence weights, and append
drawcols = ['draw_' + str(i) for i in range(1000)]
files = os.listdir("FILEPATH")
files = [s for s in files if "yld" in s]

vepafs = pd.DataFrame()
for file in files:
    print(files.index(file)/len(files))
    vepaf = pd.read_csv("FILEPATH")

    vepaf = vepaf.loc[~vepaf['age_group_id'].isin([2,3])]
    vepaf = vepaf.loc[vepaf['sex_id'] != 2]
    vepaf = vepaf.loc[vepaf['year_id'] == 2019]

    # apply incidence weights
    vepaf = agg_ages(vepaf, [42,179,239,25,26])
	vepaf = vepaf.merge(iweights, on = ['age_group_id', 'agg_age_group_id', 'location_id'])
    vepaf[drawcols] = vepaf[drawcols].multiply(vepaf['incidence_weight'], axis="index")

    # combine across draws
    vepaf = vepaf.groupby(['agg_age_group_id', 'location_id', 'year_id'])[drawcols].agg('sum').reset_index()
    vepaf['draw_mean'] = vepaf[drawcols].mean(axis = 1)
    vepaf['draw_se'] = vepaf[drawcols].sem(axis = 1)
    vepaf = vepaf.drop(drawcols, axis = 1)
    vepafs = vepafs.append(vepaf)

pathogens = pd.read_csv('FILEPATH')
altdorm = pathogens.loc[pathogens['infectious_syndrome'] == 'respiratory_infectious', 'pathogens'].str.replace('streptococcus_pneumoniae, ', '')
altdorm = altdorm.str.replace(', ', '-').values[0]
vepafs['log_ratio'] = np.log((1-vepafs['draw_mean'])/vepafs['draw_mean'])
vepafs['mod_se'] = vepafs['draw_se']/vepafs['draw_mean']
vepafs['hosp_continuous'] = 0
vepafs['pathogen_x'] = 'streptococcus_pneumoniae'
vepafs['pathogen_y'] = altdorm
vepafs['study_type'] = 'strep_ve'
vepafs['sex_id'] = 3

config = pd.read_excel("FILEPATH", sheet_name = 'run')
covs = [cov for cov in str(config.loc[0,'covariates']).split(',') if cov not in ['agg_age_group_id', 'hosp_continuous']]

for cov in covs:
    vepafs = merge_covariate(
                            vepafs, cov, decomp_step=CONF.get_id("decomp_step"),
                            gbd_round_id=CONF.get_id("gbd_round")
                        )
vepaf2 = vepafs.loc[:,['year_id', 'location_id', 'sex_id', 'agg_age_group_id', 'study_type', 'pathogen_x', 'pathogen_y', 'log_ratio', 'mod_se', 'hosp_continuous', 'haqi', 'PCV3_coverage_prop', 'Hib3_coverage_prop']]

vepaf2.to_csv("FILEPATH", index=False)
