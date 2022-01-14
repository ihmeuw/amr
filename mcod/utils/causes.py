import pandas as pd
from db_queries import get_cause_metadata
from cod_prep.claude.configurator import Configurator
from cod_prep.utils import (
    report_if_merge_fail, get_function_results, report_duplicates
)
from cod_prep.downloaders import (
    add_cause_metadata, get_best_cause_hierarchy_version,
    get_current_cause_hierarchy, get_cod_ages, get_age_weights, add_population
)

CONF = Configurator('standard')


def _get_sepsis_cause_hierarchy(int_cause, cause_set_version_id):
    assert 'sepsis' in int_cause, "this hierarchy is just for sepsis"
    df = pd.read_excel(
        "FILEPATH",
        header=1
    )
    nested_name_dict = df[
        ['Nested code ', 'acause']
    ].dropna().drop_duplicates().set_index('Nested code ')['acause'].to_dict()

    df = df.loc[df['level'] != 'XX']
    df = df.loc[~((df['cause_id'] == 605) & (
        df['parent nested code'] == 2400))]
    df['cause_id'] = df['cause_id'].replace({371: 995, 488: 490})
    df['parent nested code'] = df['parent nested code'].fillna(
        df['Nested code '])
    assert not df['cause_id'].duplicated().values.any()

    chh = get_current_cause_hierarchy(
        cause_set_version_id=cause_set_version_id,
        force_rerun=False, block_rerun=True
    )
    chh = chh.loc[
        (chh['yld_only'] != 1) & (chh['parent_id']!= 294)
        & (chh['is_estimate'] == 1)
    ]
    chh['level_2'] = chh['cause_id']
    chh.loc[chh['secret_cause'] == 1, 'level_2'] = chh['parent_id']

    level_dict = df.set_index('cause_id')['parent nested code'].to_dict()
    chh['level_1'] = chh['cause_id'].map(level_dict)

    chh.loc[chh['cause_id'] == 1003, 'level_1'] = 1300
    chh.loc[chh['cause_id'] == 1004, 'level_1'] = 1800
    chh.loc[chh['cause_id'].isin([1018, 1019]), 'level_1'] = 1000
    chh.loc[chh['cause_id'] == 1020, 'level_1'] = 1400
    chh.loc[chh['cause_id'].isin([1006, 1007]), 'level_1'] = 1700
    chh.loc[chh['cause_id'].isin(
        [1005] + list(range(1008, 1014))), 'level_1'] = 1600
    report_if_merge_fail(chh, 'level_1', 'cause_id')
    chh['level_1'] = chh['level_1'].astype(int)

    chh.loc[chh['cause_id'].isin([343, 407, 1003, 960]), 'level_1'] = 1000
    chh.loc[chh['cause_id'] == 505, 'level_1'] = 1800
    chh.loc[chh['cause_id'] == 720, 'level_1'] = 2600
    chh.loc[chh['cause_id'] == 940, 'level_1'] = 2500
    check_df = chh.drop_duplicates(['level_1', 'level_2'])
    report_duplicates(check_df, 'level_2')

    chh['nested_name'] = chh['level_1'].map(nested_name_dict)
    report_if_merge_fail(chh, 'nested_name', 'level_1')

    chh = chh[['cause_id', 'level_1', 'level_2',
               'acause', 'nested_name']].drop_duplicates()

    return chh


def _get_int_cause_hierarchy(int_cause, cause_set_version_id):
    chh = get_current_cause_hierarchy(
        cause_set_version_id=cause_set_version_id, force_rerun=False, block_rerun=True
    )
    chh['level_2'] = chh['cause_id']
    chh.loc[chh['secret_cause'] == 1, 'level_2'] = chh['parent_id']
    path_to_parent_df = chh['path_to_top_parent'].str.split(',', expand=True)
    chh['level_1'] = path_to_parent_df[2]
    assert chh.loc[chh['level_1'].isnull(), 'level'].isin([0, 1]).values.all()
    chh = chh.loc[chh['level'] >= 2]
    chh = chh[['cause_id', 'level_1', 'level_2', 'acause']].drop_duplicates()
    return chh


def get_int_cause_hierarchy(int_cause, cause_set_id=4, gbd_round_id=None,
                            cause_set_version_id=None, **cache_kwargs):
    if not gbd_round_id:
        gbd_round_id = CONF.get_id('gbd_round')

    if not cause_set_version_id:
        cause_set_version_id = CONF.get_id('cause_set_version')

    if 'sepsis' in int_cause:
        function = _get_sepsis_cause_hierarchy
    else:
        function = _get_int_cause_hierarchy
    cache_name = f"{int_cause}_nested_cause_hierarchy_v{cause_set_version_id}"

    args = [int_cause, cause_set_version_id]
    kwargs = {}

    df = get_function_results(
        function,
        args,
        kwargs,
        cache_name,
        **cache_kwargs
    )

    return df


def get_most_detailed_inj_causes(int_cause, cause_set_version_id=None,
                                 cause_set_id=3, gbd_round_id=None,
                                 **cache_kwargs):
    cause_df = get_current_cause_hierarchy(
        cause_set_version_id=cause_set_version_id,
        cause_set_id=cause_set_id, gbd_round_id=None, **cache_kwargs)
    cause_df = cause_df.loc[
        (cause_df['yld_only'] != 1) & (cause_df['most_detailed'] == 1) &
        (cause_df["secret_cause"] != 1)
    ]
    inj_causes = list(
        cause_df.loc[cause_df['acause'].str.contains('inj')].cause_id.unique())
    restricted_targets = [729, 945]
    if int_cause == "x59":
        restricted_targets += [721, 723, 725, 726, 727, 854, 941]
    return list(set(inj_causes) - set(restricted_targets))


def agg_secret_causes(df, **cache_kwargs):
    df = add_cause_metadata(
        df, ['secret_cause'], cause_set_version_id=CONF.get_id(
            'cause_set_version'), **cache_kwargs
    )
    chh = get_current_cause_hierarchy(
        cause_set_id=4, gbd_round_id=None, **cache_kwargs
    )
    sch = chh.loc[(chh["acause"].str.contains("inj")) & (
        chh["secret_cause"] == 1)]
    secret_id = sch.set_index('cause_id')['parent_id'].to_dict()
    parent_name = sch.set_index('parent_id')['acause_parent'].to_dict()
    df['parent_id'] = df['cause_id'].map(secret_id)
    df['acause_parent'] = df['parent_id'].map(parent_name)
    df.loc[(df["secret_cause"] == 1) & (
        df["cause_id"] != 743), "cause_id"] = df["parent_id"]
    df.loc[(df["secret_cause"] == 1) & (df["cause_id"] != 743),
           "acause"] = df["acause_parent"]
    df.drop(["secret_cause", "parent_id", "acause_parent"],
            axis=1, inplace=True)

    return df


def get_unique_patterns(df, col):
    df[col] = df[col].str.split("_")
    df.loc[df[col].isnull(), col] = ""
    df[col] = [list(set(x)) for x in df[col]]
    df[col] = df[col].apply(lambda x: [_f for _f in x if _f])
    df[col] = df[col].astype(str)
    df[col] = df[col].apply(lambda x: ' '.join(sorted(x.split()))).str.replace(
        "]", "").str.replace("[", ""). str.replace(
        ",", "").str.replace("'", "")

    return df


def get_injuries_cause_ids(int_cause, cache_kwargs):

    d = get_most_detailed_inj_causes(int_cause,
                                     cause_set_version_id=CONF.get_id(
                                         'reporting_cause_set_version'),
                                     **cache_kwargs)

    p = list(get_level_cause_dict(cache_kwargs, int_cause,
                                  cause_set_version_id=CONF.get_id(
                                      'reporting_cause_set_version')).values())
    p = [x for sub_list in p for x in sub_list]

    causes = get_cause_metadata(gbd_round_id=5, cause_set_version_id=4)
    inj_causes = causes.loc[causes.cause_id.isin(d + p)].cause_id.unique()

    return inj_causes


def combine_cause_sets(end_product,**cache_kwargs):

    subset_string = 'yld_only != 1'
    if end_product in ['incidence', 'attributable_burden']:
        subset_string += ' & yll_only != 1'

    chh = get_current_cause_hierarchy(
        cause_set_version_id=CONF.get_id('computation_cause_set_version'),
        cause_set_id=CONF.get_id('computation_cause_set'), **cache_kwargs
    ).query(f'{subset_string}')[
        ['level', 'cause_id', 'parent_id', 'most_detailed']]
    levels = list(chh.level.unique())
    chh_rep = get_current_cause_hierarchy(
        cause_set_version_id=CONF.get_id('reporting_cause_set_version'),
        cause_set_id=CONF.get_id('reporting_cause_set'), **cache_kwargs
    )
    reporting_aggs = list(chh_rep.loc[
        chh_rep['acause'].str.contains('reporting'), 'cause_id'].unique())
    reporting_agg_df = get_current_cause_hierarchy(
        cause_set_version_id=CONF.get_id('reporting_cause_agg_set_version'),
        cause_set_id=CONF.get_id('reporting_cause_agg_set'), **cache_kwargs
    ).query(f'{subset_string} & parent_id in {reporting_aggs}')[
        ['cause_id', 'parent_id', 'most_detailed']
    ]
    reporting_agg_df = reporting_agg_df.query("cause_id != parent_id")
    reporting_agg_df = reporting_agg_df.query("cause_id not in (1009, 1010)")
    reporting_agg_df = reporting_agg_df.append(
        {'cause_id': 1008, 'parent_id': 1022, 'most_detailed': 0},
        ignore_index=True
    )
    reporting_agg_df = reporting_agg_df.merge(
        chh[['cause_id', 'level']], how='left', on='cause_id'
    )
    assert (reporting_agg_df.groupby('parent_id')['level'].nunique() == 1).all()
    reporting_agg_df = reporting_agg_df.append(
        reporting_agg_df[['parent_id', 'level']].drop_duplicates()
        .rename(columns={'parent_id': 'cause_id'})
        .assign(level=lambda d: d['level'] - 1, most_detailed=0)
    )
    chh = chh.append(reporting_agg_df)
    return chh, levels


def get_level_parent_child_cause_dict(end_product, **cache_kwargs):
    chh, levels = combine_cause_sets(end_product, **cache_kwargs)
    level_dict = {}
    for level in levels:
        child_dict = {}
        level_df = chh.query(f"level == {level} & most_detailed == 0")
        for cause_id in level_df.cause_id.unique():
            children = list(chh.query(
                f"parent_id == {cause_id} & cause_id != {cause_id}"
            ).cause_id.unique())
            child_dict.update({cause_id: children})
        level_dict.update({level: child_dict})
    return level_dict


def get_child_causes(parent_cause_id, end_product, **cache_kwargs):
    cause_meta_df, levels = combine_cause_sets(end_product, **cache_kwargs)
    child_cause_list = list(cause_meta_df.query(
        f'parent_id == {parent_cause_id}'
    )['cause_id'].unique())
    if parent_cause_id in child_cause_list:
        child_cause_list.remove(parent_cause_id)
    return child_cause_list


def get_infsyn_hierarchy():
    return pd.read_csv("FILEPATH")


def get_child_to_available_parent_syndrome(available_parents, infsyn=None):
    if infsyn is None:
        infsyn = get_infsyn_hierarchy()
    infsyn['parent_candidates'] = infsyn['path_to_top_parent'].apply(
        lambda x: [s for s in x.split(',') if s in available_parents])
    infsyn = infsyn.loc[~infsyn.parent_candidates.apply(lambda x: len(x) == 0)]
    infsyn['parent'] = infsyn['parent_candidates'].apply(lambda x: x[-1])
    return infsyn.set_index('infectious_syndrome')['parent'].to_dict()


def get_all_related_syndromes(syndrome, infsyn=None):
    if infsyn is None:
        infsyn = get_infsyn_hierarchy()
    return infsyn.loc[
        infsyn.path_to_top_parent.apply(lambda x: syndrome in x.split(',')),
        'infectious_syndrome'
    ].unique().tolist()
