import pandas as pd
from cod_prep.claude.configurator import Configurator
from cod_prep.downloaders.causes import (
    get_current_cause_hierarchy,
    get_all_related_causes
)
from cod_prep.utils.misc import report_if_merge_fail
import warnings

CONF = Configurator()


def get_bug_drugs():
    bug_drug_universe = pd.read_csv("FILEPATH")
    possible_bug_drug_combos = list(
        bug_drug_universe[['pathogen', 'abx_class']].to_records(index=False))
    possible_bug_drug_combos = [tuple(bug_drug)
                                for bug_drug in possible_bug_drug_combos]
    return possible_bug_drug_combos


def get_pred_ex(codcorrect_version):
    if not codcorrect_version:
        codcorrect_version = CONF.get_id("codcorrect_version")
    return pd.read_hdf("FILEPATH")


def get_prepped_csbg_universe(burden_type, add_parents=False, **cache_kwargs):
    ch = get_current_cause_hierarchy(
        cause_set_version_id=CONF.get_id('computation_cause_set_version'),
        cause_set_id=CONF.get_id('computation_cause_set'),
        **cache_kwargs
    ).query("yld_only != 1")
    if burden_type == 'nonfatal':
        ch = ch.query("yll_only != 1")
    csbg_universe = pd.read_csv("FILEPATH")
    if burden_type == 'fatal':
        csbg_universe = csbg_universe.query("yld_only != 1")
    elif burden_type == 'nonfatal':
        csbg_universe = csbg_universe.query("yll_only != 1")
    csbg_universe = csbg_universe.drop(['yll_only', 'yld_only'], axis='columns')
    csbg_universe['cause_id'] = csbg_universe['acause'].apply(
        get_all_related_causes, cause_meta_df=ch
    )
    csbg_universe = csbg_universe.explode('cause_id').reset_index(drop=True)
    report_if_merge_fail(csbg_universe, 'cause_id', 'acause')
    csbg_universe.loc[(csbg_universe['abx_class'].isna()), 'abx_class'] = 'none_tested'
    csbg_universe.loc[(csbg_universe['infectious_syndrome'] == \
        'others_and_non_bacterial_infectious'), 'pathogen'] = '(none_estimated)'
    if add_parents:
        csbg_universe['path_to_top_parent'] = csbg_universe[
            'cause_id'
        ].map(ch.set_index('cause_id')['path_to_top_parent'].to_dict())
        assert csbg_universe.path_to_top_parent.notnull().values.all()
        csbg_universe['path_id'] = csbg_universe['path_to_top_parent'].str.split(',')
        csbg_universe = csbg_universe.explode('path_id').reset_index(drop=True)
        csbg_universe['path_id'] = csbg_universe['path_id'].astype(int)
        csbg_universe.loc[
            csbg_universe.acause != '_all', 'cause_id'
        ] = csbg_universe['path_id']
        csbg_universe.drop(
            ['path_to_top_parent', 'path_id'], axis='columns',
            inplace=True
        )
        csbg_universe.drop_duplicates(inplace=True)

        csbg_universe = csbg_universe.append(
            csbg_universe.assign(infectious_syndrome='all')
        ).drop_duplicates()

        csbg_universe['pathogen_aggregate'] = 0
        csbg_universe = csbg_universe.append(
            csbg_universe.loc[
                csbg_universe.abx_class.notnull()
            ].assign(pathogen='all', pathogen_aggregate=1)
        ).drop_duplicates()
        csbg_universe = csbg_universe.append(
            csbg_universe.loc[
                csbg_universe.abx_class.notnull() &
                csbg_universe.pathogen.isin(
                    ['enteropathogenic_escherichia_coli', 'enterotoxigenic_escherichia_coli']
                )
            ].assign(pathogen='escherichia_coli', pathogen_aggregate=1)
        ).drop_duplicates()
    return csbg_universe


def drop_unmodeled_age_causes(df):
    aggregate_age = df.age_group_id.isin([22, 27])
    too_old = df['simple_age'] > df['age_end']
    too_young = df['simple_age'] < df['age_start']
    unmodeled_ages = (too_old | too_young) & ~aggregate_age
    df = df[~unmodeled_ages]
    return df


def drop_unmodeled_sex_causes(df):
    return df.query("~((sex_id == 1 & male == 0) | (sex_id == 2 & female == 0))")


def drop_unmodeled_asc(df, cause_meta_df, age_meta_df, burden_type):
    assert burden_type in ['fatal', 'nonfatal']
    is_fatal = burden_type == 'fatal'
    col_prefix = {True: 'yll', False: 'yld'}
    age_start = f'{col_prefix[is_fatal]}_age_start'
    age_end = f'{col_prefix[is_fatal]}_age_end'
    only = f'{col_prefix[not is_fatal]}_only'
    cause_meta_df = cause_meta_df[[
        age_start, age_end, 'male', 'female', 'cause_id', only
    ]].rename(columns={age_start: 'age_start', age_end: 'age_end', only: 'only'})
    age_meta_df = age_meta_df[['simple_age', 'age_group_id']]
    return (
        df.merge(cause_meta_df, how="left", on="cause_id")
        .merge(age_meta_df, how="left")
        .pipe(drop_unmodeled_sex_causes)
        .pipe(drop_unmodeled_age_causes)
        .query("only != 1")
        .drop(columns=["simple_age", "age_start", "age_end", "male", "female", "only"])
    )
