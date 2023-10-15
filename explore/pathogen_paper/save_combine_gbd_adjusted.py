# imports

import re
import os
import getpass
import glob
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import importlib

from cod_prep.downloaders import get_pop, get_current_location_hierarchy
from cod_prep.claude.configurator import Configurator
import db_queries.api.public as db

from cod_prep.utils import print_log_message
from cod_prep.downloaders import get_ages
from amr_prep.utils.amr_io import get_amr_results
from mcod_prep.utils.mcause_io import get_mcause_results
from db_queries import get_outputs as go

# directories

DIR = 'FILEPATH'
repo_dir = 'FILEPATH'
MAP_DIR = repo_dir + "/maps/"

suffix = '_draws'
paths_to_combine = [f'{DIR}gbd_causes_adjusted_regions{suffix}.csv', f'{DIR}gbd_causes_adjusted_global{suffix}.csv', f'{DIR}gbd_causes_adjusted_countries{suffix}.csv']
def combine_gbd_adjuted_datasets(df_paths):
	'''
	combine several csv resuls for gbd_fix
	'''
	dfs = []
	for i in df_paths:
		df = pd.read_csv(i)
		dfs.append(df)

	df_total = pd.concat(dfs)
	return df_total

def save_gbd_adjusted(df, output_dir, suffix):
	df.to_csv(f'{output_dir}gbd_causes_adjusted{suffix}.csv', index=False)

if __name__ == '__main__':
	total_results = combine_gbd_adjuted_datasets(paths_to_combine)
	save_gbd_adjusted(total_results, DIR, suffix)