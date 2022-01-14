from __future__ import division
import os
from builtins import str
from builtins import range
import sys
import pandas as pd
from cod_prep.downloaders import (
    get_current_location_hierarchy, add_location_metadata, add_population,
    create_age_bins, get_pop, get_age_weights, get_ages,
    getcache_age_aggregate_to_detail_map
)
from cod_prep.claude.configurator import Configurator
from cod_prep.utils import report_duplicates, report_if_merge_fail, print_log_message
from mcod_prep.utils.mcause_io import McauseResult

CONF = Configurator('standard')
NUM_DRAWS = 1000
DEM_COLS = ['cause_id', 'location_id', 'sex_id', 'year_id', 'age_group_id']
DRAW_COLS = ["draw_" + str(x) for x in range(0, NUM_DRAWS)]
COD_AGES = CONF.get_id('cod_ages')


def get_lvl2_lvl0_aggregates(df, pop_df=None, draw_cols=DRAW_COLS, id_cols=DEM_COLS,
                             num_draws=NUM_DRAWS):
    df = add_population(df, pop_df=pop_df)
    report_if_merge_fail(df, 'population', ['age_group_id', 'sex_id', 'location_id', 'year_id'])

    print_log_message("Creating global aggregates")
    global_df = df.copy()
    global_df['location_id'] = 1
    global_df = scale_location_aggregates(
        global_df, pop_df=pop_df, draw_cols=draw_cols, id_cols=id_cols,
        num_draws=num_draws)

    print_log_message("Creating region aggregates")
    region_df = df.copy()
    region_df['location_id'] = region_df['region_id']
    region_df = scale_location_aggregates(
        region_df, pop_df=pop_df, draw_cols=draw_cols, id_cols=id_cols,
        num_draws=num_draws)

    print_log_message("Creating super region aggregates")
    super_region_df = df.copy()
    super_region_df['location_id'] = super_region_df['super_region_id']
    super_region_df = scale_location_aggregates(
        super_region_df, pop_df=pop_df, draw_cols=draw_cols, id_cols=id_cols,
        num_draws=num_draws
    )

    scaled_locs_df = pd.concat([global_df, region_df, super_region_df],
                               ignore_index=True, sort=True)

    return scaled_locs_df


def scale_location_aggregates(df, pop_df=None, draw_cols=DRAW_COLS, id_cols=DEM_COLS,
                              num_draws=NUM_DRAWS):
    df = df.copy()
    value_cols = draw_cols + ['population']
    assert 'population' in df.columns

    df = df.groupby(
        id_cols, as_index=False
    )[value_cols].sum()

    for draw in range(0, num_draws):
        df['rate_' + str(draw)] = df[draw_cols[draw]] / df['population']

    df = df.drop('population', axis=1)
    df = add_population(df, pop_df=pop_df)
    report_if_merge_fail(df, 'population', ['age_group_id', 'sex_id', 'location_id', 'year_id'])

    for draw in range(0, num_draws):
        df[draw_cols[draw]] = df['rate_' + str(draw)] * df['population']

    rate_draws = ["rate_" + str(x) for x in range(0, num_draws)]
    df = df.drop(rate_draws + ['population'], axis=1)
    return df


def get_country_aggregate(df, lhh, draw_cols=DRAW_COLS, id_cols=DEM_COLS):
    country_agg_df = add_location_metadata(df, ['level', 'iso3'], location_meta_df=lhh)
    report_if_merge_fail(country_agg_df, 'iso3', 'location_id')
    country_agg_df = country_agg_df.query("level > 3")

    group_cols = ['iso3'] + [x for x in id_cols if x != 'location_id']
    country_agg_df = country_agg_df.groupby(group_cols, as_index=False)[draw_cols].sum()

    ihme_loc_id_dict = lhh.set_index('ihme_loc_id')['location_id'].to_dict()
    country_agg_df['location_id'] = country_agg_df['iso3'].map(ihme_loc_id_dict)
    country_agg_df.drop('iso3', axis=1, inplace=True)

    return country_agg_df


def aggregate_locs(df, lhh, draw_cols=DRAW_COLS, id_cols=DEM_COLS, num_draws=NUM_DRAWS):
    df = df.copy()
    assert len(draw_cols) == num_draws

    print_log_message("Reading in population")
    pop_df = get_pop(pop_run_id=CONF.get_id('pop_run'), force_rerun=False, block_rerun=True)
    pop_df = pop_df.query(f'age_group_id in {COD_AGES}')

    data_ages = list(df.age_group_id.unique())
    if len(set(data_ages) - set(COD_AGES)) > 0:
        print_log_message(f"Aggregating population age groups to match: {data_ages} ids")
        pop_df = create_age_bins(pop_df, data_ages, dropna=True)
        pop_df = pop_df.groupby(
            ['age_group_id', 'sex_id', 'location_id', 'year_id'], as_index=False
        )['population'].sum()

    df = add_location_metadata(df, ['region_id', 'super_region_id'], location_meta_df=lhh)
    lvl2_lvl0_df = get_lvl2_lvl0_aggregates(
        df, pop_df=pop_df, draw_cols=draw_cols, id_cols=id_cols, num_draws=num_draws)
    country_df = get_country_aggregate(df, lhh, draw_cols=draw_cols, id_cols=id_cols)

    loc_agg_df = pd.concat([lvl2_lvl0_df, country_df], ignore_index=True, sort=True)

    return loc_agg_df


def create_age_standardized_rates(df, draw_cols, dem_cols, location_set=40):
    cache_kwargs = {'force_rerun': False, 'block_rerun': True}
    df = df.copy()
    df = df.loc[df['age_group_id'].isin(COD_AGES)]
    age_weight_df = get_age_weights(force_rerun=False, block_rerun=True)
    age_weight_dict = age_weight_df.drop_duplicates(
        ['age_group_id', 'age_group_weight_value']
    ).set_index('age_group_id')['age_group_weight_value'].to_dict()
    df['weight'] = df['age_group_id'].map(age_weight_dict)
    report_if_merge_fail(df, 'weight', 'age_group_id')
    pop_run_id = CONF.get_id('pop_run')
    df = add_population(df, pop_run_id=pop_run_id, **cache_kwargs)
    null_pop = df['population'].isnull()
    null_pop_df = df[null_pop].drop('population', axis=1)
    if len(null_pop_df) > 0:
        null_pop_df = add_population(null_pop_df, pop_run_id=pop_run_id,
                                     location_set_id=location_set, **cache_kwargs)
        df = pd.concat([df[~null_pop], null_pop_df], ignore_index=True)
    report_if_merge_fail(df, "population", ['sex_id', 'age_group_id',
                                            'year_id', 'location_id'])
    if type(draw_cols) != list:
        draw_cols = [draw_cols]
    for draw_col in draw_cols:
        df[draw_col] = (df[draw_col] / df['population']) * df['weight']
    group_cols = [x for x in dem_cols if x != 'age_group_id']
    df = df.groupby(group_cols, as_index=False)[draw_cols].sum()
    df['age_group_id'] = 27

    return df


def aggregate_ages(df):
    age_ids = [22, 1, 158, 23, 159, 24, 25, 26, 21, 28, 157, 42, 162, 420]
    age_map = getcache_age_aggregate_to_detail_map(
        force_rerun=False, block_rerun=True, cache_results=False)
    age_map = age_map.loc[age_map.agg_age_group_id.isin(age_ids)]
    df = df.copy()
    df = df.merge(age_map, how='inner', on='age_group_id')
    df['age_group_id'] = df['agg_age_group_id']
    assert df['age_group_id'].notnull().all()
    return df.groupby(DEM_COLS, as_index=False)[DRAW_COLS].sum()


def aggregate_sexes(df):
    df = df.copy()
    group_cols = [x for x in DEM_COLS if x != 'sex_id']
    all_sex_df = df.groupby(group_cols, as_index=False)[DRAW_COLS].sum()
    all_sex_df['sex_id'] = 3
    return all_sex_df


def main(int_cause, description, end_product, cause_id, year_id):
    df = McauseResult(
        int_cause=int_cause, end_product=end_product, process='calculate_counts',
        year_id=year_id, cause_id=cause_id, description=description, conf=CONF
    ).read_results()
    print_log_message("Creating location aggregates")
    lhh = get_current_location_hierarchy(
        location_set_version_id=CONF.get_id('location_set_version'),
        force_rerun=False, block_rerun=True
    )
    loc_agg_df = aggregate_locs(df, lhh)
    df = pd.concat([df, loc_agg_df], ignore_index=True, sort=True)
    report_duplicates(df, DEM_COLS)

    print_log_message("Creating sex aggregates")
    all_sex_df = aggregate_sexes(df)
    df = pd.concat([all_sex_df, df], ignore_index=True, sort=True)

    if end_product in ['mortality', 'incidence', 'attributable_burden']:
        print_log_message("Creating age aggregates")
        agg_age_df = aggregate_ages(df)
        asr_df = create_age_standardized_rates(df, DRAW_COLS, DEM_COLS)
        df = pd.concat([df, asr_df, agg_age_df], ignore_index=True, sort=True)

    report_duplicates(df, DEM_COLS)
    df[DEM_COLS] = df[DEM_COLS].astype(int)
    McauseResult(
        int_cause=int_cause, end_product=end_product, process='age_loc_aggregation',
        year_id=year_id, cause_id=cause_id, description=description, conf=CONF
    ).write_results(df)


if __name__ == "__main__":
    description = str(sys.argv[1])
    end_product = str(sys.argv[2])
    int_cause = str(sys.argv[3])
    task_id = os.environ.get('SGE_TASK_ID')
    if task_id:
        print(f'Running as array job, task_id: {task_id}')
        base_out_dir = McauseResult(
            int_cause=int_cause, end_product=end_product,
            process='run_model', description=description, conf=CONF
        ).parent_model_dir
        task_row = pd.read_csv("FILEPATH").iloc[int(task_id) - 1]
        cause_id = int(task_row['cause_id'])
        year_id = int(task_row['year_id'])
    else:
        cause_id = int(sys.argv[4])
        year_id = int(sys.argv[5])
    print_log_message(f'Running year {year_id}, cause {cause_id}')
    main(int_cause, description, end_product, cause_id, year_id)
