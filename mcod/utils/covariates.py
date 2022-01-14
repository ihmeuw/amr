from __future__ import division
from builtins import range
import pandas as pd
import numpy as np
from sklearn import preprocessing
from db_tools import ezfuncs
from db_queries import get_covariate_estimates, get_best_model_versions
from cod_prep.claude.configurator import Configurator
from cod_prep.utils import (
    report_if_merge_fail, report_duplicates, drop_unmodeled_sex_causes, drop_unmodeled_asc,
    get_function_results, print_log_message
)
from cod_prep.downloaders import (
    add_cause_metadata, add_age_metadata, add_location_metadata,
    get_current_location_hierarchy, create_age_bins, add_population
)

CONF = Configurator('standard')

def get_int_cause_model_covariates(int_cause):
    covariates_df = pd.read_csv("FILEPATH").query('int_cause == @int_cause')
    assert len(covariates_df) == 1
    return covariates_df['covariates'].str.split(', ').iloc[0]


def get_covariate_id(covariate):
    query = """SELECT * FROM shared.covariate WHERE covariate_name_short = '{}'""".format(covariate)
    df = ezfuncs.query(query, conn_def=CONF.get_database_setup('db'))
    assert len(df) == 1
    return int(df['covariate_id'].iloc[0])


def get_covariate_id_cols(covariate_id):
    query = """SELECT * FROM shared.covariate WHERE covariate_id = {}""".format(covariate_id)
    df = ezfuncs.query(query, conn_def=CONF.get_database_setup('db'))
    assert len(df) == 1
    id_cols = ['location_id', 'year_id']
    if df['by_age'].iloc[0] == 1:
        id_cols += ['age_group_id']
    if df['by_sex'].iloc[0] == 1:
        id_cols += ['sex_id']
    return id_cols


def add_weight_to_covariate(cov_df, covariate_name_short, **cache_kwargs):
    if covariate_name_short is None:
        covariate_name_short = cov_df.loc[0, 'covariate_name_short']

    cov_to_weight_map = {
        'alcohol_lpc': 'population', 'SEV_metab_fpg': 'population',
        'SEV_wash_sanitation': 'population', 'smok_prev': 'population',
        'prev_obesity': 'population', 'idu_prop_agestandardized': 'none',
        'idu_prop_byage': 'population'
    }
    weight = cov_to_weight_map[covariate_name_short]

    if weight == 'none':
        cov_df['weight'] = 1
    elif weight == 'population':
        cov_df = add_population(cov_df, pop_run_id=CONF.get_id('pop_run'),
                                location_set_id=CONF.get_id('location_set'),
                                decomp_step=CONF.get_id('decomp_step'), **cache_kwargs)
        report_if_merge_fail(cov_df, check_col='population',
                             merge_cols=['age_group_id', 'location_id', 'year_id', 'sex_id'])
        cov_df.rename({'population': 'weight'}, axis='columns', inplace=True)
    else:
        raise NotImplementedError

    return cov_df


def aggregate_covariate(cov_df, agg_col, covariate_name_short=None,
                        agg_age_group_ids=None, id_cols=None, **cache_kwargs):
    if covariate_name_short is None:
        covariate_name_short = cov_df.loc[0, 'covariate_name_short']
    if id_cols is None:
        covariate_id = get_covariate_id(covariate_name_short)
        id_cols = get_covariate_id_cols(covariate_id)

    assert agg_col in ['sex_id', 'age_group_id']

    cov_df = add_weight_to_covariate(cov_df, covariate_name_short, **cache_kwargs)

    cov_df['weighted_cov'] = cov_df['mean_value'] * cov_df['weight']

    if agg_col == 'age_group_id':
        cov_df = create_age_bins(cov_df, agg_age_group_ids, dropna=False)
    elif agg_col == 'sex_id':
        assert set(cov_df.sex_id) == {1, 2}
        cov_df['sex_id'] = 3

    cov_df = cov_df.groupby(id_cols, as_index=False)[['weighted_cov', 'weight']].sum()

    cov_to_agg_type_map = {
        'alcohol_lpc': 'average', 'SEV_metab_fpg': 'average', 'SEV_wash_sanitation': 'average',
        'smok_prev': 'average', 'prev_obesity': 'average', 'idu_prop_byage': 'average',
        'idu_prop_agestandardized': 'average'
    }
    agg_type = cov_to_agg_type_map[covariate_name_short]

    if agg_type == 'average':
        cov_df.eval("agg_cov = weighted_cov / weight", inplace=True)
    elif agg_type == 'sum':
        cov_df['agg_cov'] = cov_df['weighted_cov']
    else:
        raise NotImplementedError

    cov_df.drop(['weighted_cov', 'weight'], axis='columns', inplace=True)
    cov_df.rename({'agg_cov': 'mean_value'}, axis='columns', inplace=True)
    assert set(cov_df.columns) == set(id_cols + ['mean_value'])
    assert cov_df.notnull().values.all()

    return cov_df


def merge_covariate(df, covariate_name_short, scale=False, **get_cov_kwargs):
    covariate_id = get_covariate_id(covariate_name_short)
    id_cols = get_covariate_id_cols(covariate_id)
    get_cov_kwargs.update(
        {'location_id': list(df.location_id.unique()),
         'year_id': list(df.year_id.unique())}
    )
    cov_df = get_cov(covariate_id=covariate_id, **get_cov_kwargs)[id_cols + ['mean_value']]
    cache_kwargs = {'force_rerun': False, 'block_rerun': True}

    if 'sex_id' in id_cols:
        if not set(df.sex_id).issubset(set(cov_df.sex_id)):
            assert set(df.sex_id) == {3}
            assert set(cov_df.sex_id) == {1, 2}
            cov_df = aggregate_covariate(
                cov_df, 'sex_id', covariate_name_short=covariate_name_short,
                id_cols=id_cols, **cache_kwargs
            )

    if 'age_group_id' in id_cols:
        df_age_group_ids = df.age_group_id.unique().tolist()
        cov_age_group_ids = cov_df.age_group_id.unique().tolist()

        if not set(df_age_group_ids).issubset(set(cov_age_group_ids)):
            print_log_message("Aggregating covariates to match incoming dataframe.")
            cov_df = add_age_metadata(
                cov_df, ['age_group_days_start', 'age_group_days_end'], **cache_kwargs
            )
            df_ag = add_age_metadata(
                df.copy(), ['age_group_days_start', 'age_group_days_end'], **cache_kwargs
            )
            too_young = cov_df['age_group_days_end'] <= df_ag.age_group_days_start.min()
            too_old = cov_df['age_group_days_start'] >= df_ag.age_group_days_end.max()
            cov_df = cov_df[~(too_young | too_old)]
            cov_df = cov_df.drop(['age_group_days_start', 'age_group_days_end'], axis=1)
            cov_df = aggregate_covariate(
                cov_df, 'age_group_id', covariate_name_short=covariate_name_short,
                agg_age_group_ids=df_age_group_ids, id_cols=id_cols, **cache_kwargs
            )

    cov_df = cov_df.rename(columns={'mean_value': covariate_name_short})
    report_duplicates(cov_df, id_cols)
    print_log_message("RUNNING MERGE COVARIATES")
    if covariate_name_short == 'haqi':
        cov_df.eval(f"{covariate_name_short} = {covariate_name_short} / 100", inplace=True)
    if (covariate_name_short == 'LDI_pc') or ('ldi_pc' in covariate_name_short):
        print_log_message('using natural log of LDI per capita')
        cov_df[covariate_name_short] = np.log(cov_df[covariate_name_short])

    if scale:
        scaler = preprocessing.MinMaxScaler()
        cov_df[[covariate_name_short]] = scaler.fit_transform(cov_df[[covariate_name_short]])

    df = df.merge(cov_df, on=id_cols, how='left')
    report_if_merge_fail(df, covariate_name_short, id_cols)
    return df


def drop_unmodeled_asc_aggregate_ages(df, cause_meta_df, age_meta_df, detail_ages):
    simple_age_map = (
        age_meta_df.loc[lambda d: d["age_group_id"].isin(detail_ages), :]
        .set_index("simple_age", verify_integrity=True)
    )
    df = add_cause_metadata(
        df, add_cols=['cause_start', 'cause_end', 'male', 'female', 'yld_only'],
        cause_meta_df=cause_meta_df.assign(
            cause_start=lambda d: d['yll_age_start'].map(simple_age_map['age_group_years_start']),
            cause_end=lambda d: d['yll_age_end'].map(simple_age_map['age_group_years_end']),
        )
    )
    df = add_age_metadata(
        df, add_cols=['age_group_years_start', 'age_group_years_end'],
        age_meta_df=age_meta_df
    )
    df = drop_unmodeled_sex_causes(df)

    unmodeled_age_old = df["age_group_years_start"] >= df["cause_end"]
    unmodeled_age_young = df["age_group_years_end"] <= df["cause_start"]
    df = df[~(unmodeled_age_old | unmodeled_age_young | (df['yld_only'] == 1))]
    df.drop(['cause_start', 'cause_end', 'male', 'female', 'yld_only',
             'age_group_years_start', 'age_group_years_end'], axis=1, inplace=True)
    return df


def validate_redistribution_restrictions(df, cause_meta_df):
    valid_ages = [0, 0.01, 0.1, 1] + list(range(5, 95, 5))
    bad_age_starts = set(df.age_start.dropna()) - set(valid_ages)
    assert bad_age_starts == set(
    ), "Age starts in redstribution restrictions are invalid {}".format(bad_age_starts)
    bad_age_ends = set(df.age_start.dropna()) - set(valid_ages)
    assert bad_age_ends == set(), \
        "Age ends in redstribution restrictions are invalid {}".format(bad_age_ends)
    valid_locations = set(get_current_location_hierarchy(
        location_set_version_id=CONF.get_id('location_set_version'),
        force_rerun=False, block_rerun=True
    ).location_name).union(set(["NONE"]))
    all_locations = set(df.super_region.dropna()).union(
        set(df.region.dropna())).union(
        set(df.country.dropna())).union(
        set(df.subnational_level1.dropna())).union(
        set(df.subnational_level2.dropna()))
    bad_locations = set(all_locations) - set(valid_locations)
    assert bad_locations == set(
    ), "Locations in redstribution restrictions are invalid {}".format(bad_locations)
    valid_years = list(range(1900, 2100))
    bad_years = set(df.year_start.dropna()).union(
        set(df.year_end.dropna())) - set(valid_years)
    assert bad_years == set(), \
        "years in redstribution restrictions are invalid {}".format(bad_years)
    valid_dev_status = ["D0", "D1"]
    bad_dev_status = set(df.dev_status.dropna()) - set(valid_dev_status)
    assert bad_dev_status == set(
    ), "Dev_status in redstribution restrictions are invalid {}".format(bad_dev_status)
    bad_causes = set(df.acause) - set(cause_meta_df.acause)
    assert bad_causes == set(
    ), "Acauses in redstribution restrictions don't exist " \
        "in cause hierarchy {}".format(bad_causes)


def enforce_redistribution_restrictions(int_cause, cause_meta_df, age_meta_df, template):
    if int_cause not in ['x59', 'y34']:
        df = pd.read_csv("FILEPATH")
    else:
        df = pd.read_csv("FILEPATH")
    validate_redistribution_restrictions(df, cause_meta_df)
    df['super_region'] = df['super_region'].replace({'NONE': None})
    assert df[
        ['region', 'country', 'subnational_level2', 'subnational_level1',
         'year_start', 'year_end', 'dev_status', 'sex']
    ].isnull().values.all(), 'write out other needed restrictions'
    df = df[['acause', 'age_start', 'age_end', 'super_region']]
    df['age_start'] = df['age_start'].fillna(0)
    df['age_end'] = df['age_end'].fillna(95)
    ignore_acauses = CONF.get_id('ignore_redistribution_overrides')
    df = df.loc[~(df['acause'].isin(ignore_acauses))]
    df = add_cause_metadata(df, add_cols='cause_id', merge_col='acause',
                            cause_meta_df=cause_meta_df)
    orig_cols = template.columns
    template = template.merge(df, on='cause_id', how='left')
    if not template.acause.isnull().values.all():
        template = add_location_metadata(
            template, location_set_version_id=CONF.get_id('location_set_version'),
            force_rerun=False, block_rerun=True, add_cols='super_region_name'
        )
        template = add_age_metadata(
            template, add_cols=['age_group_years_start', 'age_group_years_end'],
            age_meta_df=age_meta_df
        )
        too_young = template['age_start'] >= template['age_group_years_end']
        too_old = template['age_end'] < template['age_group_years_start']
        template['super_region'] = template['super_region'].fillna(template['super_region_name'])
        loc_violation = template['super_region_name'] != template['super_region']
        template = template[~(too_young | too_old | loc_violation)]
    template = template[orig_cols]
    return template


def get_cov(covariate_id=None, model_version_id=None, location_set_id=None, gbd_round_id=None,
            decomp_step=None, location_id='all', year_id='all', **cache_kwargs):
    if location_set_id is None:
        location_set_id = CONF.get_id('location_set')
    if gbd_round_id is None:
        gbd_round_id = CONF.get_id('gbd_round')
    if decomp_step is None:
        decomp_step = CONF.get_id('decomp_step')
    if model_version_id is None:
        model_version_id = get_best_model_versions(
            entity='covariate', ids=covariate_id, gbd_round_id=gbd_round_id,
            decomp_step=decomp_step, status='best'
        ).loc[0, 'model_version_id']
    cache_name = f"cov_{covariate_id}_mvid_{model_version_id}_lsid_{location_set_id}"
    function = get_covariate_estimates
    args = [covariate_id]
    kwargs = {
        'location_set_id': location_set_id,
        'gbd_round_id': gbd_round_id,
        'decomp_step': decomp_step,
        'model_version_id': model_version_id,
        'location_id': location_id,
        'year_id': year_id
    }
    df = get_function_results(
        function,
        args,
        kwargs,
        cache_name,
        **cache_kwargs
    )
    if isinstance(df, pd.DataFrame):
        assert covariate_id == df.loc[0, 'covariate_id'],\
            f"Covariate {covariate_id} and model version {model_version_id} do not match."
    return df
