import os
import sys
import getpass
import argparse
user = getpass.getuser()
amrpath = FILEPATH
if not amrpath in sys.path:
    sys.path.append(FILEPATH)
import pandas as pd
import numpy as np
from warnings import warn
from cod_prep.claude.claude_io import Configurator
from mcod_prep.utils.causes import get_infsyn_hierarchy, get_all_related_syndromes
from mcod_prep.utils.mcause_io import get_mcause_data
from amr_prep.utils.amr_io import get_amr_data
from amr_prep.utils.pathogen_formatting import PathogenFormatter
from cod_prep.utils import (
    print_log_message, report_duplicates, create_square_df,
    report_if_merge_fail
)
from db_queries import get_covariate_estimates, get_location_metadata
from db_tools import ezfuncs
from mcod_prep.utils.nids import add_nid_metadata, get_datasets
from mcod_prep.utils.covariates import get_cov, merge_covariate
from cod_prep.claude.relative_rate_split import relative_rate_split
from cod_prep.downloaders import (
    add_age_metadata, get_current_location_hierarchy,
    add_location_metadata, get_cod_ages, pretty_print,
    get_pop, add_population, getcache_age_aggregate_to_detail_map,
    get_country_level_location_id, get_ages,
    prep_age_aggregate_to_detail_map
)

CONF = Configurator()

parser = argparse.ArgumentParser(description="remap all or only new data")
parser.add_argument("syndrome", help="specific syndrome to compile, 'all' compiles every syndrome",)
arg = parser.parse_args()

if arg.syndrome == 'all':
	syndromes = ['cns_infectious', 'uti_plus', 'skin_infectious', 'blood_stream_infectious',
            'respiratory_infectious', 'cardiac_infectious', 'peritoneal_and_intra_abdomen_infectious',
            'bone_joint_infection', 'diarrhea']
else:
	syndromes = [arg.syndrome]

syndparams = pd.read_csv(FILEPATH)

for inf_synd in syndromes:
    print(inf_synd)
    syndpaths = syndparams.loc[syndparams['infectious_syndrome'] == inf_synd, 'pathogens'].str.split(", ").tolist()[0]
    excludeICU = not syndparams.loc[syndparams['infectious_syndrome'] == inf_synd, 'ICU_in_CFR'].values[0]

    formatter = PathogenFormatter(model_type = 'cfr', infectious_syndrome = inf_synd, 
                                  keep_pathogens = syndpaths,
                                  cache_kwargs = {'force_rerun': False, 'block_rerun': False, 'cache_results': False})
    formatdat = formatter.format_data()
    
    if inf_synd == 'cns_infectious':
        mening_lit = pd.read_csv(FILEPATH)
        mening_lit['pathogen'] = mening_lit['pathogen'].str.lower()
        mening_lit.loc[mening_lit['pathogen'] == 'other_unsp', 'pathogen'] = 'other'
        mening_lit.rename(columns={'nid': 'source', 'admissions':'cases'}, inplace=True)
        mening_lit = mening_lit.drop(['agg_admissions', 'data_type', 'agg_deaths'], axis = 1)

        mening_lit.loc[(mening_lit['year_id'] == 1979) & (mening_lit['location_id'] == 93), 'year_id'] = 1980
        mening_lit.loc[(mening_lit['year_id'] == 1979) & (mening_lit['location_id'] == 89), 'year_id'] = 1980

        LS_ID = CONF.get_id("location_set")
        locs = get_location_metadata(location_set_id = LS_ID, release_id = 10)
        vstlnd_2_norway = np.int(locs.loc[locs['location_id'] == 60132, 'parent_id'])
        mening_lit.loc[mening_lit['location_id'] == 60132, 'location_id'] = vstlnd_2_norway
        mening_lit = mening_lit.loc[mening_lit['cases'].notna(),:]
        mening_lit = mening_lit.loc[mening_lit['deaths'].notna(),:]
        mening_lit['ICU'] = 'mixed'

        formatdat = formatdat.append(mening_lit)
    
    format_w_covar = merge_covariate(
        formatdat, 'haqi', decomp_step=CONF.get_id("decomp_step"),
        gbd_round_id=CONF.get_id("gbd_round"))
    outpath = FILEPATH
    os.makedirs(outpath, exist_ok = True)
    format_w_covar.to_csv(path_or_buf = FILEPATH, index = False)                              