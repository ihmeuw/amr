from builtins import str
from builtins import range
import os
import sys
import pandas as pd
import numpy as np
from get_draws.api import get_draws
from cod_prep.claude.configurator import Configurator
from cod_prep.downloaders import (
    get_current_location_hierarchy, create_age_bins, get_ages,
    get_current_cause_hierarchy, get_cause_map
)
from cod_prep.downloaders.causes import get_all_related_causes
from cod_prep.utils import print_log_message
from mcod_prep.mcod_mapping import MCoDMapper
from mcod_prep.utils.mcause_io import McauseResult
from mcod_prep.utils.causes import (
    get_all_related_syndromes, get_int_cause_hierarchy
)

DEM_COLS = ['cause_id', 'location_id', 'sex_id', 'year_id', 'age_group_id']
DRAW_COLS = ["draw_" + str(x) for x in range(0, 1000)]
CONF = Configurator('standard')
BLOCK_RERUN = {'force_rerun': False, 'block_rerun': True}
BAD_CAUSES = [657, 656, 668]


def convert_int_cols(df):
    df[DEM_COLS] = df[DEM_COLS].apply(pd.to_numeric, downcast='integer')
    return df


def get_location_list():
    lhh = get_current_location_hierarchy(
        location_set_version_id=CONF.get_id('location_set_version'), **BLOCK_RERUN
    ).query("most_detailed == 1")
    return list(lhh['location_id'].unique())


def get_codcorrect_deaths(cause_id, year_id, age_group_id, location_id, sex_id,
                          num_workers, gbd_round_id=None, release_id=None,
                          version_id=None):
    if gbd_round_id is None:
        gbd_round_id = CONF.get_id("gbd_round")
    if release_id is None:
        release_id = CONF.get_id("release")
    if version_id is None:
        version_id = CONF.get_id("codcorrect_version")
    get_draws_kwargs = {
        'gbd_id_type': 'cause_id', 'gbd_id': cause_id, 'source': 'codcorrect',
        'metric_id': 1, 'measure_id': 1, 'location_id': location_id, 'age_group_id': age_group_id,
        'year_id': year_id, 'gbd_round_id': gbd_round_id, 'version_id': version_id,
        'num_workers': num_workers, 'release_id': release_id, 'sex_id': sex_id
    }
    print_log_message("Reading draws of deaths from CoDCorrect")
    df = get_draws(**get_draws_kwargs)
    return df


def get_deaths_draws(cause_id, year_id, ages, locs, sexes, int_cause, end_product,
                     gbd_round_id=None, release_id=None, version_id=None,
                     get_all_ages=True):

    if gbd_round_id is None:
        gbd_round_id = CONF.get_id("gbd_round")
    if release_id is None:
        release_id = CONF.get_id("release")
    if version_id is None:
        version_id = CONF.get_id("codcorrect_version")
    is_syndrome = int_cause in MCoDMapper.infectious_syndromes
    get_codcorrect_kwargs = {
        'cause_id': cause_id, 'location_id': locs, 'age_group_id': ages,
        'year_id': year_id, 'gbd_round_id': gbd_round_id, 'version_id': version_id,
        'num_workers': 5, 'release_id': release_id, 'sex_id': sexes
    }
    if (end_product in ['incidence', 'attributable_burden']):
        mort_result = McauseResult(
            int_cause=int_cause, end_product='mortality', process='calculate_counts',
            year_id=year_id, cause_id=cause_id, status='latest',
            conf=CONF
        )
        print_log_message(
            f"Reading {int_cause} deaths from\n {mort_result.parent_model_dir}")
        df = mort_result.read_results()
        if end_product == 'attributable_burden':
            print_log_message("Reading draws of deaths from CoDCorrect")
            cc_df = get_codcorrect_deaths(**get_codcorrect_kwargs)
            df = df.merge(cc_df, how='left', on=DEM_COLS)
            print_log_message(f"Calculating non {int_cause}-related deaths")
            for draw in DRAW_COLS:
                df[draw] = df[f'{draw}_y'] - df[f'{draw}_x']
            assert (df[DRAW_COLS] >= 0).values.all()
    elif end_product == 'mortality' and is_syndrome:
        sepsis_result = McauseResult(
            int_cause='sepsis', end_product='mortality',
            process='cause_aggregation', year_id=year_id, cause_id=cause_id,
            status='latest', conf=CONF
        )
        print_log_message(
            f"Reading sepsis deaths from\n {sepsis_result.parent_model_dir}"
        )
        df = sepsis_result.read_results()
    elif (end_product in ['mortality', 'rdp']):
        print_log_message("Reading draws of deaths from CoDCorrect")
        df = get_codcorrect_deaths(**get_codcorrect_kwargs)
    missing_ages = set(ages) - set(df.age_group_id)
    if get_all_ages:
        assert missing_ages == set(), f'Missing ages: {missing_ages}'
    else:
        if missing_ages != set():
            print_log_message(
                f"Missing the following age groups from the denominator: {missing_ages}")
    assert df[DEM_COLS +
              DRAW_COLS].notnull().values.all(), "Error pulling draws of deaths"
    return df.groupby(DEM_COLS, as_index=False)[DRAW_COLS].sum()


def gbd2019_skin_fix(df):
    assert CONF.get_id('gbd_round') == 6
    yll_age_start = get_current_cause_hierarchy(
        cause_set_version_id=CONF.get_id('computation_cause_set_version'),
        cause_set_id=CONF.get_id('computation_cause_set'),
        **BLOCK_RERUN).query('cause_id == 653').reset_index().at[0, 'yll_age_start']
    drop_ages = get_ages(**BLOCK_RERUN).query(
        f'age_group_years_start < {yll_age_start} & age_group_years_end < {yll_age_start}'
    ).age_group_id.unique()
    df = df.loc[~(df['age_group_id'].isin(drop_ages))]
    return df


def convert_prop_to_count(prop_df, deaths_df, end_product):
    print_log_message("Calculating counts")
    if end_product in ['mortality', 'rdp']:
        df = prop_df.merge(deaths_df, how='left', on=DEM_COLS)
    elif end_product in ['incidence', 'attributable_burden']:
        prop_df[DRAW_COLS] = np.reciprocal(prop_df[DRAW_COLS])
        df = prop_df.merge(deaths_df, how='inner', on=DEM_COLS)
    for draw in DRAW_COLS:
        df[draw] = df[f'{draw}_x'] * df[f'{draw}_y']
    df = df[DEM_COLS + DRAW_COLS]
    assert df.notnull().values.all(), "Error calculating deaths"
    return df


def squeeze(int_cause, end_product, cause_id, prop_df):

    chh = get_current_cause_hierarchy(cause_set_version_id=CONF.get_id("computation_cause_set_version"),
                                      cause_set_id=CONF.get_id("computation_cause_set"), force_rerun=False, block_rerun=True)

    GBD_result_cause = ['tb', 'intest_typhoid', 'intest_paratyph',
                        'intest_ints', 'std_chlamydia', 'std_gonnorhea']
    zeroed_causes = []
    for cause in GBD_result_cause:
        GBD_causes = get_all_related_causes(cause, chh)
        zeroed_causes = zeroed_causes + GBD_causes

    single_syn_cause_ids = {
        332: ['cns_infectious', 'blood_stream_infectious'],
        503: ['cardiac_infectious'],
        529: ['peritoneal_and_intra_abdomen_infectious'],
        322: ['respiratory_infectious_comm'],
        980: ['skin_infectious'],
        302: ['diarrhea'],
        595: ['uti_plus_comm'],
        368: ['blood_stream_infectious'],
        383: ['blood_stream_infectious']
    }
    single_syn_dictionary = {}
    for parent_cause_id, key in single_syn_cause_ids.items():
        child_causes = get_all_related_causes(parent_cause_id, chh)
        single_syn_dictionary.update(
            {child_cause: key for child_cause in child_causes}
        )

    all_sepsis_cause_ids = GBD_result_cause + \
        list(single_syn_cause_ids) + ['uri']
    all_sepsis_causes = []
    for cause in all_sepsis_cause_ids:
        sepsis_causes = get_all_related_causes(cause, chh)
        all_sepsis_causes = all_sepsis_causes + sepsis_causes

    comm_inf_syn = [
        syn for syn in MCoDMapper.infectious_syndromes if ("_comm" in syn)]
    comm_cause_list = []

    inf_map_10 = pd.read_csv("FILEPATH")
    inf_map_10['code_system_id'] = 1
    inf_map_9 = pd.read_csv("FILEPATH")
    inf_map_9['code_system_id'] = 6
    inf_map = pd.concat([inf_map_10, inf_map_9], ignore_index=True, sort=True)

    cause_map_10 = get_cause_map(1, force_rerun=False, block_rerun=True)
    cause_map_9 = get_cause_map(6, force_rerun=False, block_rerun=True)
    cause_map = pd.concat([cause_map_10, cause_map_9],
                          ignore_index=True, sort=True)
    cause_map.rename(columns={'value': 'icd_code'}, inplace=True)

    cause_meta_df = get_int_cause_hierarchy(
        int_cause, force_rerun=False, block_rerun=True, cache_dir='standard',
        cause_set_version_id=CONF.get_id('cause_set_version')
    )[['cause_id', 'level_1', 'level_2']]

    for syn in comm_inf_syn:
        syn = syn.replace("_comm", "")
        comm_syn = get_all_related_syndromes(syn)
        comm_syn = inf_map.loc[inf_map['infectious_syndrome'].isin(
            comm_syn), ['icd_code', 'code_system_id', 'infectious_syndrome']]

        comm_syn = comm_syn.merge(
            cause_map[['icd_code', 'cause_id', 'code_system_id']], how='left', on=['icd_code', 'code_system_id'])

        comm_syn = comm_syn[~comm_syn['icd_code'].isin(['U07.1', 'U07.2'])]
        assert comm_syn['cause_id'].notnull().values.all()

        comm_syn['int_cause'] = syn + '_comm'
        comm_cause_list.append(comm_syn)

    community_pairs = pd.concat(comm_cause_list, ignore_index=True, sort=True)
    community_pairs = community_pairs.merge(
        cause_meta_df, how='left', on='cause_id')
    community_pairs.loc[(community_pairs['level_2'].isna()) & (
        community_pairs['cause_id'] == 743), 'level_2'] = community_pairs['cause_id']

    assert community_pairs['level_2'].notnull().values.all()

    if int_cause == 'sepsis':
        if cause_id in all_sepsis_causes:
            draw_cols = [x for x in list(prop_df) if ("draw" in x)]
            prop_df[draw_cols] = 1
    else:
        if cause_id not in zeroed_causes:
            version_df = pd.read_csv("FILEPATH").query('int_cause != @int_cause')

            prop_df['int_cause'] = int_cause
            prediction_dfs = [prop_df]

            for index in range(len(version_df)):
                print(f"Appending infectious syndrome: {version_df.iloc[index,0]}")
                predict_int_cause = version_df.iloc[index, 0]
                predict_description = version_df.iloc[index, 1]

                predict_df = McauseResult(
                    int_cause=predict_int_cause, end_product=end_product, process='predict',
                    year_id=year_id, cause_id=cause_id, description=predict_description, conf=CONF
                ).read_results()

                if 'level_2' in predict_df.columns:
                    predict_df.rename(
                        columns={'level_2': 'cause_id'}, inplace=True)
                    convert_int_cols(predict_df)

                if (end_product in ['mortality', 'rdp']) & (cause_id in BAD_CAUSES):
                    print_log_message('Fix for GBD2019 skin diseases')
                    predict_df = gbd2019_skin_fix(predict_df)

                predict_df['int_cause'] = predict_int_cause

                prediction_dfs.append(predict_df)

            prediction_dfs = pd.concat(
                prediction_dfs, ignore_index=True, sort=True)

            prediction_dfs = prediction_dfs[
                [x for x in list(prediction_dfs)
                    if ("_id" in x)] + ["level_1","int_cause"] + [x for x in list(prediction_dfs) if ("draw" in x)
                ]
            ]

            all_possible = MCoDMapper.infectious_syndromes
            model_result_causes = [e for e in all_possible if e not in ('infectious_syndrome',
                                                                        'tb', 'typhoid_paratyphoid_ints', 'chlamydia_and_gonorrheae')]
            assert set(prediction_dfs.int_cause) == set(model_result_causes)

            predict_cols = [x for x in list(
                prediction_dfs) if ("id" in x)] + ["level_1"]
            draw_cols = [x for x in list(prediction_dfs) if ("draw" in x)]

            if cause_id in single_syn_dictionary:
                possible_syndromes = single_syn_dictionary[cause_id]
                prediction_dfs.loc[
                    ~prediction_dfs['int_cause'].isin(possible_syndromes),
                    draw_cols
                ] = 0

            for inf_syn in community_pairs['int_cause'].unique():
                CAI_causes = community_pairs.loc[community_pairs['int_cause'] == inf_syn, 'level_2'].unique(
                ).tolist()

                if cause_id not in CAI_causes:
                    prediction_dfs.loc[prediction_dfs['int_cause']
                                       == inf_syn, draw_cols] = 0

            prediction_dfs[draw_cols] = pd.DataFrame(
                prediction_dfs[draw_cols].to_numpy(
                ) / prediction_dfs.groupby(predict_cols)[draw_cols].transform(sum).to_numpy(),
                index=prediction_dfs.index
            )

            prediction_dfs = prediction_dfs.loc[prediction_dfs['int_cause'] == int_cause]
            prediction_dfs = prediction_dfs.drop(['int_cause'], axis=1)

            prop_df = prediction_dfs
        else:
            prop_df = prop_df[[x for x in list(prop_df) if (
                "_id" in x)] + ["level_1"] + [x for x in list(prop_df) if ("draw" in x)]]

            draw_cols = [x for x in list(prop_df) if ("draw" in x)]
            prop_df[draw_cols] = 0

    print_log_message("Saving squeeze")
    McauseResult(
        int_cause=int_cause, end_product=end_product, process='squeeze',
        year_id=year_id, cause_id=cause_id, description=description, conf=CONF
    ).write_results(prop_df)

    return prop_df


def main(description, cause_id, year_id, end_product, int_cause, no_squeeze):
    prop_df = McauseResult(
        int_cause=int_cause, end_product=end_product, process='predict',
        year_id=year_id, cause_id=cause_id, description=description, conf=CONF
    ).read_results()
    if 'level_2' in prop_df.columns:
        prop_df.rename(columns={'level_2': 'cause_id'}, inplace=True)
    convert_int_cols(prop_df)
    if (end_product in ['mortality', 'rdp']) & (cause_id in BAD_CAUSES):
        print_log_message('Fix for GBD2019 skin diseases')
        prop_df = gbd2019_skin_fix(prop_df)

    if not no_squeeze and int_cause in MCoDMapper.infectious_syndromes + ['sepsis'] and end_product == 'mortality':
        prop_df = squeeze(int_cause, end_product, cause_id, prop_df)

    ages = list(prop_df.age_group_id.unique())
    locs = list(prop_df.location_id.unique())
    sexes = list(prop_df.sex_id.unique())
    draws_df = get_deaths_draws(
        cause_id, year_id, ages, locs, sexes, int_cause, end_product)
    df = convert_prop_to_count(prop_df, draws_df, end_product)
    print_log_message("Saving output")
    McauseResult(
        int_cause=int_cause, end_product=end_product, process='calculate_counts',
        year_id=year_id, cause_id=cause_id, description=description, conf=CONF
    ).write_results(df)


if __name__ == '__main__':
    description = str(sys.argv[1])
    int_cause = str(sys.argv[2])
    end_product = str(sys.argv[3])
    no_squeeze = str(sys.argv[4])

    assert no_squeeze in ["True", "False"]
    no_squeeze = (no_squeeze == "True")

    task_id = os.environ.get('SGE_TASK_ID')
    if task_id:
        print(f'Running as array job, task_id: {task_id}')
        base_out_dir = McauseResult(
            int_cause=int_cause,
            end_product=end_product,
            process='run_model',
            description=description,
            conf=CONF
        ).parent_model_dir
        task_row = pd.read_csv("FILEPATH").iloc[int(task_id) - 1]
        cause_id = int(task_row['cause_id'])
        year_id = int(task_row['year_id'])
    else:
        cause_id = int(sys.argv[5])
        year_id = int(sys.argv[6])
    print_log_message(f'Running year {year_id}, cause {cause_id}')
    print_log_message(f'No Squeeze: {no_squeeze}')
    main(description, cause_id, year_id, end_product, int_cause, no_squeeze)
