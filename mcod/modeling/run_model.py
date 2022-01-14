import os
import sys
import re
import math
import argparse
import pandas as pd
from mcod_prep.utils.covariates import merge_covariate
from mcod_prep.utils.mcause_io import get_mcause_data
from mcod_prep.utils.mcod_cluster_tools import submit_mcod
from mcod_prep.utils.nids import get_datasets
from mcod_prep.utils.causes import get_int_cause_hierarchy
from cod_prep.utils import report_if_merge_fail, print_log_message
from cod_prep.claude.configurator import Configurator
from cod_prep.downloaders import add_cause_metadata
from mcod_prep.utils.mcause_io import McauseResult, makedirs_safely
from cod_prep.downloaders import (
    create_age_bins, add_cause_metadata, add_location_metadata, pretty_print
)

CONF = Configurator('standard')


def fix_sepsis(df, int_cause):
    assert int_cause in ['explicit_sepsis', 'implicit_sepsis', 'sepsis']
    sepsis_type = int_cause.split('_')[0]
    if sepsis_type == 'sepsis':
        df[int_cause] = (df['sepsis'].isin(['explicit', 'implicit'])) * 1
    else:
        df[int_cause] = (df['sepsis'] == sepsis_type) * 1


def collapse_across_value_cols(df, end_product, int_cause, id_cols):
    value_cols=['successes', 'failures']
    denominator = 'cases'
    if end_product == 'attributable_burden':
        df['successes'] = (df[int_cause] != 1) * df['deaths']
        df['cases'] = (df[int_cause] != 1) * df['admissions']
    else:
        df['successes'] = (df[int_cause] == 1) * df['deaths']
        if end_product == 'incidence':
            df['cases'] = (df[int_cause] == 1) * df['admissions']
        elif end_product == 'mortality':
            denominator = 'deaths'
    if denominator == 'cases':
        df = df[df['cases'] > 0]
    df['failures'] = df[denominator] - df['successes']
    value_cols += [denominator]
    print_log_message(f"Collapsing {value_cols} across {id_cols}")
    df = df.groupby(id_cols, as_index=False)[value_cols].sum()
    df['obs_fraction'] = df['successes'] / df[denominator]
    return df


def format_for_model(df, int_cause, end_product, age_group_ids,
                     id_cols=['year_id', 'sex_id', 'age_group_id',
                              'location_id', 'level_1', 'level_2']):
    drop_rows = ((df['sex_id'] == 9) | ~(df['age_group_id'].isin(CONF.get_id('cod_ages'))) |
                 (df['cause_id'].isin([919, 744, 743])))
    df = df[~drop_rows]

    df = create_age_bins(df, age_group_ids)

    df = add_cause_metadata(
        df, ['yld_only', 'yll_only'], block_rerun=True, force_rerun=False,
        cause_set_version_id=CONF.get_id('cause_set_version')
    )
    if end_product == 'mortality':
        df = df[(df['deaths'] > 0) & (df['yld_only'] != 1)]
    else:
        df = df[(df['admissions'] > 0) & (df['yll_only'] != 1)]
    if int_cause in ['implicit_sepsis', 'explicit_sepsis', 'sepsis']:
        fix_sepsis(df, int_cause)
    assert int_cause in df.columns, "intermediate cause is missing!"
    assert set(df[int_cause].unique()) == {0, 1}, \
        f"expecting {int_cause} column to be 0 or 1"
    df = collapse_across_value_cols(df, end_product, int_cause, id_cols)
    df[id_cols] = df[id_cols].astype(int)
    return df


def merge_nested_cause_levels(df, int_cause):
    cause_meta_df = get_int_cause_hierarchy(
        int_cause, force_rerun=False, block_rerun=True, cache_dir='standard',
        cause_set_version_id=CONF.get_id('cause_set_version')
    )[['cause_id', 'level_1', 'level_2']]
    df = df.merge(cause_meta_df, how='left', on='cause_id')
    return df


def pull_covariates(int_cause):
    covariates_df = pd.read_csv("FILEPATH").query('int_cause == @int_cause')
    assert len(covariates_df) == 1
    return covariates_df['covariates'].str.split(', ').iloc[0]


def write_model_input_data(df, int_cause, end_product, description, output_dir, age_group_id=None):
    makedirs_safely(output_dir)
    print_log_message(f"Writing model input file to {output_dir}")
    if age_group_id:
        df = df.query(f'age_group_id == {age_group_id}')
    df.to_csv("FILEPATH", index=False)


def launch_modelworker(int_cause, end_product, description, output_dir, age_group_id=None):
    diag_dir = f'{output_dir}/DIRECTORY'
    in_sample_dir = f"{output_dir}/DIRECTORY"
    oo_sample_dir = f"{output_dir}/DIRECTORY"
    makedirs_safely(in_sample_dir)
    makedirs_safely(oo_sample_dir)
    if age_group_id:
        jobname = f'{int_cause}_{end_product}_{description}_{age_group_id}_modelworker'
    else:
        jobname = f'{int_cause}_{end_product}_{description}_modelworker'
    worker = "FILEPATH/run_modelworker.R"
    params = [int_cause, output_dir, diag_dir, CONF.get_id("gbd_round")]
    submit_mcod(jobname, 'r', worker, 2, '45G', params=params, verbose=True, logging=True,
                log_base_dir=output_dir, queue='long.q', runtime='07:00:00:00')


def get_SOURCENAME(int_cause):
    df = get_mcause_data(
        phase='format_map', sub_dirs=int_cause, assert_all_available=True,
        force_rerun=True, block_rerun=False, is_active=True, source='SOURCE_NAME'
    )
    group_cols = list(set(df.columns) - set(['code_id', 'deaths']))
    return df.query('age_group_id != 6 & cause_id != 743').groupby(
        group_cols, as_index=False)['deaths'].sum()


def main(description, int_cause, end_product, age_group_ids):
    data_kwargs = {
        'phase': 'redistribution', 'sub_dirs': int_cause, 'assert_all_available': False,
        'force_rerun': True, 'block_rerun': False, 'is_active': True,
        'data_type_id': [9, 3, 13], 'year_id': list(range(1980, 2050))
    }
    if end_product in ['incidence', 'attributable_burden']:
        data_kwargs.update({'data_type_id': 3})
    print_log_message("Pulling training data")
    df = get_mcause_data(**data_kwargs)
    if end_product == 'mortality':
        print_log_message("Excluding 2000-2010 COUNTRY Hospital data for linkage")
        COUNTRY_hosp_2000_2010 = [220786, 220787, 220788, 220789, 220790, 220791,
                              220792, 220793, 220794, 220795, 220796]
        df = df.query(f'nid not in {COUNTRY_hosp_2000_2010}')
        print_log_message("Appending SOURCE_NAME")
        df = pd.concat([df, get_SOURCENAME(int_cause)], sort=True)
    df = merge_nested_cause_levels(df, int_cause)
    df = format_for_model(df, int_cause, end_product, age_group_ids)
    covariates = pull_covariates(int_cause)
    for covariate in covariates:
        print_log_message('Merging on {}'.format(covariate))
        df = merge_covariate(df, covariate)
    assert df.notnull().values.all()
    if 'by_age' in description:
        for age_group_id in age_group_ids:
            output_dir_by_age = str(McauseResult(
                int_cause=int_cause, end_product=end_product, process='run_model',
                description=description, age_group_id=age_group_id, conf=CONF
            ).results_path)
            write_model_input_data(
                df, int_cause, end_product, description, output_dir_by_age, age_group_id)
            print_log_message("Launching model")
            launch_modelworker(int_cause, end_product, description, output_dir_by_age, age_group_id)
    else:
        output_dir = str(McauseResult(
            int_cause=int_cause, end_product=end_product, process='run_model',
            description=description, conf=CONF
        ).results_path)
        write_model_input_data(df, int_cause, end_product, description, output_dir)
        print_log_message("Launching model")
        launch_modelworker(int_cause, end_product, description, output_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Format data and launch model.')
    parser.add_argument('description', type=str)
    parser.add_argument('int_cause', type=str)
    parser.add_argument('end_product', type=str)
    parser.add_argument('--age_group_ids', action='append', type=int)
    args = parser.parse_args()
    main(**vars(args))