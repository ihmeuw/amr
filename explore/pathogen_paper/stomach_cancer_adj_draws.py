
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
from get_draws.api import get_draws
from functools import reduce


DIR = 'FILEPATH'


lh=get_current_location_hierarchy()

locs = lh.loc[lh["level"].isin([3,2,1]), "location_id"].unique().tolist()+[1]
draws_dir=  'FILEPATH'
pops = get_pop(release_id=6, pop_run_id=192)

dalys_sc = get_draws("cause_id",414,
	               source="dalynator",
	        metric_id=1,
	        measure_id=2,
	        release_id=6,
	        location_id=locs,
	        sex_id=[1,2],
	        year_id=2019,
	       age_group_id=[1,22])

deaths_sc = get_draws("cause_id",414,
	               source="codcorrect",
	        metric_id=1,
	        measure_id=1,
	        release_id=6,
	        location_id=locs,
	        sex_id=[1,2],
	        year_id=2019,
	       age_group_id=[1,22])



x = 1000
draw_name = [f'draw_{i}' for i in list(range(0, x))]

deaths_sc_u5 = deaths_sc.copy()
deaths_sc_u5.loc[:,draw_name]=0
deaths_sc_u5.loc[:,'age_group_id']=1


df_sc = pd.concat([dalys_sc, deaths_sc, deaths_sc_u5])

IDs = ["measure_id",
       "location_id", 
       "metric_id",
       "sex_id",
    "age_group_id",
    "cause_id"
 ]

df_sc= df_sc.melt(
    id_vars=IDs,
    value_vars=draw_name, var_name='draws', value_name='val'
)

male_scalar = 0.8088235294117647
female_scalar = 0.8670520231213873

df_sc.loc[df_sc.sex_id==1,'val'] = df_sc.val*male_scalar
df_sc.loc[df_sc.sex_id==2,'val'] = df_sc.val*female_scalar

# combine the counts for both sexes 
df_sc = df_sc.groupby(['measure_id','location_id','metric_id','age_group_id','cause_id','draws']).sum().reset_index()

st_hpl = pd.read_csv(f'{draws_dir}hpylori_stomachcancer_draws.csv')
st_hpl['draw_name'] = [f'draw_{i}' for i in list(range(0, x))]

df_sc = df_sc[['measure_id','location_id', 'metric_id','age_group_id','cause_id','draws', 'sex_id','val']]
df_sc = df_sc.merge(st_hpl[['draw','draw_name']], how='left', left_on='draws', right_on= 'draw_name')


df_sc['val_adj'] = df_sc.val*df_sc.draw
df_sc = df_sc[['measure_id', 'location_id', 'metric_id', 'age_group_id', 'cause_id',
   'draws', 'sex_id', 'val_adj']]

df_sc = df_sc[['measure_id', 'location_id', 'metric_id', 'age_group_id', 'cause_id', 'sex_id', 'val_adj','draws']]


df_sc['pathogen'] = 'Stomach cancer (H.pylori)'

df_sc  = df_sc.pivot(index=['measure_id', 'location_id', 'metric_id', 'age_group_id', 'cause_id', 'sex_id','pathogen'], columns='draws').reset_index()

df_sc.columns = df_sc.columns.droplevel(1)

draws_name=[f'draw_{i}' for i in range(1000)]

cols = []
count = 0
for column in df_sc.columns:
    if column == 'val_adj':
        cols.append(f'draw_{count}')
        count+=1
        continue
    cols.append(column)


df_sc.columns = cols


df_sc.to_csv(f'{DIR}stomach_cancer_adj_draws.csv', index=False)