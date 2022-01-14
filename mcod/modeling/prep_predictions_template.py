import os
import pandas as pd
import argparse
import itertools
from db_queries import get_demographics
from cod_prep.utils import print_log_message
from mcod_prep.utils.causes import get_int_cause_hierarchy
from mcod_prep.mcod_mapping import MCoDMapper
from cod_prep.claude.configurator import Configurator
from cod_prep.downloaders import (
    get_current_cause_hierarchy, get_ages, add_cause_metadata, add_age_metadata,
    get_current_location_hierarchy
)
from cod_prep.downloaders.ages import getcache_age_aggregate_to_detail_map
from mcod_prep.utils.covariates import (
    drop_unmodeled_asc_aggregate_ages, enforce_redistribution_restrictions,
    get_int_cause_model_covariates, merge_covariate
)
from mcod_prep.utils.mcause_io import (
    McauseResult, makedirs_safely
    )
CONF = Configurator()

def drop_unmodeled_asc(df, cause_meta_df, age_meta_df, yld=False):
    if yld:
        prefix = 'yld'
    else:
        prefix = 'yll'
    add_cause_cols = [f'{prefix}_age_start', f'{prefix}_age_end',
                      'male', 'female']
    df = add_cause_metadata(df, add_cause_cols, cause_meta_df=cause_meta_df)
    df = add_age_metadata(df, ['simple_age'], age_meta_df=age_meta_df)
    df = drop_unmodeled_sex_causes(df)
    df = drop_unmodeled_age_causes(df, prefix)
    return df.drop(add_cause_cols + ['simple_age'], axis=1)


def drop_unmodeled_age_causes(df, prefix):
    too_old = df['simple_age'] > df[f'{prefix}_age_end']
    too_young = df['simple_age'] < df[f'{prefix}_age_start']
    return df[~(too_old | too_young)]


def drop_unmodeled_sex_causes(df):
    return df.query("~((sex_id == 1 & male == 0) | (sex_id == 2 & female == 0))")


def save_template(year_id, cause_id, int_cause, model_ages, end_product, description,
                  subnationals=True):
    dem_dict = get_demographics(gbd_team="cod", gbd_round_id=CONF.get_id('gbd_round'))
    assert isinstance(model_ages, list)
    dem_dict.update({'cause_id': [cause_id], 'year_id': [year_id], 'age_group_id': model_ages})
    if not subnationals:
        level_3_locations = get_current_location_hierarchy(
            location_set_version_id=CONF.get_id("location_set_version"),
            force_rerun=False, block_rerun=True, cache_results=False
        ).query("level == 3").location_id.unique().tolist()
        dem_dict.update({'location_id': level_3_locations})
    rows = itertools.product(*list(dem_dict.values()))
    template = pd.DataFrame.from_records(rows, columns=list(dem_dict.keys()))
    kwargs = {'force_rerun': False, 'block_rerun': True,
              'cause_set_version_id': CONF.get_id('cause_set_version')}
    cause_meta_df = get_current_cause_hierarchy(**kwargs)
    age_meta_df = get_ages(force_rerun=False, block_rerun=True)
    detail_ages = CONF.get_id('cod_ages')
    yld = end_product in ['incidence', 'attributable_burden']
    if end_product == 'mortality' and int_cause in MCoDMapper.infectious_syndromes: 
        if len(set(model_ages) - set(detail_ages)) != 0:
            age_detail_map = getcache_age_aggregate_to_detail_map()
            template.rename(columns={"age_group_id": "agg_age_group_id"}, inplace=True)
            template = template.merge(age_detail_map[['agg_age_group_id', 'age_group_id']], how='left', on='agg_age_group_id')
            template = drop_unmodeled_asc(template, cause_meta_df, age_meta_df, yld)
            template.rename(columns={"age_group_id":"detailed_age_group_id", "agg_age_group_id": "age_group_id"}, inplace=True)
        else: 
            template = drop_unmodeled_asc(template, cause_meta_df, age_meta_df, yld)
    else: 
        if len(set(model_ages) - set(detail_ages)) == 0:
            template = drop_unmodeled_asc(template, cause_meta_df, age_meta_df, yld)
        else:
            assert not yld, NotImplementedError
            template = drop_unmodeled_asc_aggregate_ages(
            template, cause_meta_df, age_meta_df, detail_ages)
    if end_product == 'rdp':
        template = enforce_redistribution_restrictions(
            int_cause, cause_meta_df, age_meta_df, template)
    else:
        custom_ch = get_int_cause_hierarchy(
            int_cause, **kwargs).set_index('level_2')['level_1'].to_dict()
        template['level_1'] = template['cause_id'].map(custom_ch)
        template = template.rename(columns={'cause_id': 'level_2'})
    covariates = get_int_cause_model_covariates(int_cause)
    for covariate in covariates:
        template = merge_covariate(
            template, covariate, **{'force_rerun': False, 'block_rerun': True}
        )
    result = McauseResult(
        int_cause=int_cause,
        end_product=end_product,
        process='predictions_template',
        year_id=year_id,
        cause_id=cause_id,
        description=description,
        age_group_id=model_ages[0],
        conf=CONF
    )
    result.write_results(template)
    makedirs_safely(str(result.parent_model_dir / str(year_id)))


def main(year_id, cause_id, description, int_cause, age_group_ids, end_product,
         subnationals=True):
    if 'by_age' in description:
        for age_group_id in age_group_ids:
            save_template(
                year_id, cause_id, int_cause, [age_group_id], end_product,
                description, subnationals=subnationals
            )
    else:
        save_template(
            year_id, cause_id, int_cause, age_group_ids, end_product,
            description, subnationals=subnationals
        )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Save csv that has all age/sex/location/year/causes '
        'for which to predict.'
    )
    parser.add_argument('int_cause', type=str)
    parser.add_argument('end_product', type=str)
    parser.add_argument('description', type=str)
    parser.add_argument('--age_group_ids', action='append', type=int)
    parser.add_argument('--subnationals', action='store_true', help='Predict for subnationals')
    args = parser.parse_args()
    argparse_dict = vars(args)
    task_id = os.environ.get('SGE_TASK_ID')
    if task_id:
        print(f'Running as array job, task_id: {task_id}')
        tmp_dir = McauseResult(
            int_cause=argparse_dict['int_cause'],
            end_product=argparse_dict['end_product'],
            process='run_model',
            description=argparse_dict['description'],
            age_group_id=argparse_dict['age_group_ids'][0],
            conf=CONF
        ).parent_model_dir
        task_row = pd.read_csv("FILEPATH").iloc[int(task_id) - 1]
        cause_id = int(task_row['cause_id'])
        year_id = int(task_row['year_id'])
    else:
        NotImplementedError
    argparse_dict.update({'year_id': year_id, 'cause_id': cause_id})
    print(f"Using the following arguments: {argparse_dict}")
    main(**argparse_dict)