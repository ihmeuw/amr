import sys
import os
import numpy as np
import pandas as pd
from pathlib import Path
from get_draws.api import get_draws
from cod_prep.downloaders import (
    get_current_location_hierarchy, get_cod_ages,
    get_current_cause_hierarchy, add_population,
    getcache_age_aggregate_to_detail_map
)
from cod_prep.utils import (
    report_if_merge_fail, print_log_message,
    report_duplicates
)
from cod_prep.claude.configurator import Configurator
from amr_prep.utils.amr_io import AmrResult
from amr_prep.utils.misc import drop_unmodeled_asc
from multiply.split_pathogen import PathogenSplitter


CACHE_KWARGS = {
    'force_rerun': False,
    'block_rerun': True,
    'cache_results': False
}

CONF = Configurator('standard')
GBD_ROUND = CONF.get_id("gbd_round")
RELEASE_ID = CONF.get_id("release")
LSV_ID = CONF.get_id('location_set_version')
CSV_ID = CONF.get_id('computation_cause_set_version')
CS_ID = CONF.get_id('computation_cause_set')
POP_RUN = CONF.get_id('pop_run')

CAUSE_FOR_OTHERS = [
    'diptheria',
    'encephalitis',
    'hepatitis',
    'hiv_other',
    'infectious',
    'malaria',
    'measles',
    'ntd_afrtryp',
    'ntd_chagas',
    'ntd_cysticer',
    'ntd_dengue',
    'ntd_nema_ascar',
    'ntd_other',
    'ntd_rabies',
    'ntd_schisto',
    'ntd_yellowfever',
    'otitis',
    'std_other',
    'uri',
    'varicella',
]

PROXY_CAUSES = {
    'cns_infectious': ['meningitis'],
    'cardiac_infectious': ['cvd_endo'],
    'respiratory_infectious': ['lri'],
    'skin_infectious': ['skin_bacterial', 'skin_cellulitis', 'skin_decubitus'],
    'typhoid_paratyphoid_ints': ['intest_typhoid', 'intest_paratyph', 'intest_ints'],
    'diarrhea': ['diarrhea'],
    'chlamydia_and_gonorrheae': ['std_chlamydia', 'std_gonnorhea'],
    'uti_plus': ['urinary_nephritis'],
    'blood_stream_infectious': ['maternal_sepsis', 'neonatal_sepsis'],
    'tb': ['tb_other', 'tb_mdr', 'tb_xdr'],
    'others_and_non_bacterial_infectious': CAUSE_FOR_OTHERS
}

CFR_SYNS = [
    'blood_stream_infectious',
    'peritoneal_and_intra_abdomen_infectious',
    'bone_joint_infection',
    'cardiac_infectious',
    'skin_infectious',
    'cns_infectious',
    'diarrhea',
    'respiratory_infectious',
    'uti_plus',
]


def draws_wrapper(source, causes, year_id, measure_id, metric_id, version_id,
                  num_workers):
    locs = get_current_location_hierarchy(
        location_set_version_id=LSV_ID, **CACHE_KWARGS
    )
    all_cntry = list(locs.loc[locs['level'] == 3, 'location_id'].unique())
    ages = get_cod_ages(
        gbd_round_id=GBD_ROUND, **CACHE_KWARGS
    )['age_group_id'].unique().tolist()
    draws = get_draws(
        gbd_id_type=["cause_id"],
        gbd_id=causes,
        source=source,
        year_id=year_id,
        measure_id=measure_id,
        metric_id=metric_id,
        gbd_round_id=GBD_ROUND,
        version_id=version_id,
        release_id=RELEASE_ID,
        location_id=all_cntry,
        sex_id=[1, 2],
        age_group_id=ages,
        num_workers=num_workers
    )
    return draws


def get_mi_ratios(num_workers, year_id, infectious_syndrome):
    proxy_causes = PROXY_CAUSES[infectious_syndrome]
    ch = get_current_cause_hierarchy(
        cause_set_id=CS_ID, cause_set_version_id=CSV_ID,
        **CACHE_KWARGS
    ).set_index("acause")['cause_id'].to_dict()
    proxy_cause_ids = [ch[acause] for acause in proxy_causes]

    print_log_message(f"Pulling incidence for causes {proxy_cause_ids}")
    incd = draws_wrapper(
        source='como',
        causes=proxy_cause_ids,
        year_id=year_id,
        measure_id=6,  # Incidence
        metric_id=3,  # Rates
        version_id=CONF.get_id("como_version"),
        num_workers=num_workers
    )

    print_log_message(f"Pulling deaths for causes {proxy_cause_ids}")
    mort = draws_wrapper(
        source='codcorrect',
        causes=proxy_cause_ids,
        year_id=year_id,
        measure_id=1,  # Deaths
        metric_id=1,  # Counts
        version_id=CONF.get_id("codcorrect_version"),
        num_workers=num_workers
    )

    print_log_message("Calculating MI ratios...")
    demo_cols = ['year_id', 'location_id', 'age_group_id', 'sex_id']
    incd = add_population(incd, pop_run_id=POP_RUN, **CACHE_KWARGS)
    report_if_merge_fail(incd, 'population', demo_cols)

    draw_cols = [col for col in incd if 'draw_' in col]
    incd.loc[:, draw_cols] = incd.loc[:, draw_cols].multiply(
        incd.loc[:, 'population'], axis="index")

    incd = incd[demo_cols + ['cause_id'] + draw_cols]
    mort = mort[demo_cols + ['cause_id'] + draw_cols]

    age_meta_df = get_cod_ages(gbd_round_id=GBD_ROUND, **CACHE_KWARGS)
    cause_meta_df = get_current_cause_hierarchy(
        cause_set_id=CS_ID, cause_set_version_id=CSV_ID,
        **CACHE_KWARGS
    )
    mort = drop_unmodeled_asc(mort, cause_meta_df, age_meta_df, 'fatal')
    incd = drop_unmodeled_asc(incd, cause_meta_df, age_meta_df, 'nonfatal')

    merge_cols = demo_cols + ['cause_id']
    mir = pd.merge(
        mort, incd, how='outer', on=merge_cols,
        suffixes=('_mortality', '_incidence'), indicator=True
    )
    mort_draws = [f'draw_{i}_mortality' for i in range(0, 1000)]
    incd_draws = [f'draw_{i}_incidence' for i in range(0, 1000)]

    if infectious_syndrome == 'others_and_non_bacterial_infectious':
        locs = get_current_location_hierarchy(
            location_set_version_id=LSV_ID, **CACHE_KWARGS
        )
        ssa_loc_ids = locs.loc[
            (locs['level'] == 3)
            & (locs['super_region_name'] == 'Sub-Saharan Africa'),
            'location_id'].unique().tolist()
        la_loc_ids = locs.loc[
            (locs['level'] == 3)
            & (locs['region_name'].str.contains('Latin America')),
            'location_id'].unique().tolist()
        mir = mir.loc[~((mir['cause_id'] == 300) & (mir['age_group_id'].isin([2, 3]))), ]
        mir = mir.loc[
            ~((mir['_merge'] == 'left_only')
              & (mir[incd_draws].isna().any(axis=1) | (mir[incd_draws] == 0).all(axis=1)))
        ]
        mir = mir.loc[~((mir['cause_id'] == 346) & ~(mir['location_id'].isin(la_loc_ids))), ]
        mir = mir.loc[
            ~((mir['cause_id'] == 358) & ~(mir['location_id'].isin(la_loc_ids + ssa_loc_ids))), ]
        mir = mir.loc[
            ~((mir[mort_draws].isna().any(axis=1) | (mir[mort_draws] == 0).all(axis=1))
              & (mir[incd_draws] == 0).all(axis=1)), ]
        mir.loc[mir[mort_draws].isna().any(axis=1), mort_draws] = 0
        mir = mir.groupby(demo_cols, as_index=False)[incd_draws + mort_draws].sum()

        mir[draw_cols] = pd.DataFrame(
            mir.filter(regex='.*_incidence').to_numpy()
            / mir.filter(regex='.*_mortality').to_numpy(),
            index=mir.index
        )
        assert mir[draw_cols].notnull().values.all()
        assert (mir[draw_cols] != np.inf).values.all()
        mir = mir[
            demo_cols + draw_cols
        ]
        report_duplicates(mir, demo_cols)
    else:
        assert (mir._merge != 'left_only').all()
        if infectious_syndrome == 'skin_infectious':
            mir = mir.loc[~mir.age_group_id.isin([2, 3])]
        mir['incidence_only'] = (mir._merge == 'right_only') | (mir[mort_draws] == 0).any(axis=1)
        assert mir.loc[mir.incidence_only, 'cause_id'].isin(
            [319, 320, 383, 934, 946, 947, 395]).all()
        assert ~mir.loc[
            mir.incidence_only & (mir.cause_id == 383), 'age_group_id'
        ].isin([2, 3]).any()
        assert mir.loc[mir._merge == 'both'].notnull().values.all()
        mir[mort_draws] = mir[mort_draws].fillna(0)
        assert mir.notnull().values.all()

        mir[draw_cols] = pd.DataFrame(
            mir.filter(regex='.*_incidence').to_numpy()
            / mir.filter(regex='.*_mortality').to_numpy(),
            index=mir.index
        )
        assert mir.loc[~mir.incidence_only, draw_cols].notnull().values.all()
        assert (mir.loc[~mir.incidence_only, draw_cols] != np.inf).values.all()
        if infectious_syndrome == 'skin_infectious':
            pn = mir.query("age_group_id == 4")
            mir = mir.append(pn.assign(age_group_id=2)).append(pn.assign(age_group_id=3))
        mir = mir[
            demo_cols + ['cause_id', 'incidence_only'] + draw_cols + incd_draws
        ]
        report_duplicates(mir, demo_cols + ['cause_id'])
    return mir


def get_cfr(year_id, infectious_syndrome):
    print_log_message("Reading in CFRs...")
    if infectious_syndrome == 'cardiac_infectious':
        input_dir = Path("FILEPATH")
    else:
        input_dir = Path("FILEPATH")

    cfr = pd.concat(
        [pd.read_csv(file).query(f"year_id == {year_id}")
         for file in input_dir.rglob("*.csv")],
        sort=False
    )
    assert cfr.notnull().values.all()
    cfr = cfr.rename(
        columns=lambda x: "draw_" + str(int(x.replace('draw_', "")) - 1) if 'draw' in x else x
    )

    print_log_message("Aggregating pathogens...")
    dist = PathogenSplitter.read_pathogen_models(
        conf=CONF, infectious_syndrome=infectious_syndrome,
        year_id=year_id, force_rerun=False, block_rerun=True,
        cache_results=False
    )
    for pathogen in list(set(dist.pathogen) - set(cfr.pathogen)):
        if infectious_syndrome == 'diarrhea':
            if 'escherichia_coli' in pathogen:
                use_pathogen = 'escherichia_coli'
            elif 'virus' in pathogen:
                use_pathogen = 'virus'
            else:
                use_pathogen = 'other'
        else:
            use_pathogen = 'other'
        print_log_message(f"No CFR for pathogen {pathogen}, using {use_pathogen}")
        cfr = cfr.append(
            cfr.query(f"pathogen == '{use_pathogen}'").assign(pathogen=pathogen)
        )
    if infectious_syndrome == 'diarrhea':
        cfr = cfr.loc[~cfr.pathogen.isin(['escherichia_coli', 'virus'])]

    id_cols = ['location_id', 'year_id', 'age_group_id', 'sex_id']
    if 'hosp' in cfr:
        id_cols += ['hosp']
        dist = dist.loc[dist['hosp'] != 'unknown']
        cfr = cfr.loc[cfr['hosp'] != 'unknown']
    cfr = cfr.drop_duplicates()
    draw_cols = [c for c in cfr if 'draw' in c]
    validate = 'one_to_one'
    if infectious_syndrome == 'diarrhea':
        age_map = getcache_age_aggregate_to_detail_map(**CACHE_KWARGS)
        age_map = age_map.loc[age_map.agg_age_group_id.isin(cfr.age_group_id.unique())]
        dist = dist.merge(
            age_map, how='left', on='age_group_id', validate='many_to_one'
        )
        report_if_merge_fail(dist, 'agg_age_group_id', 'age_group_id')
        id_cols.remove('age_group_id')
        dist['agg_sex_id'] = 3
        id_cols.remove('sex_id')
        id_cols += ['agg_age_group_id', 'agg_sex_id']
        cfr = cfr.rename(columns={'age_group_id': 'agg_age_group_id', 'sex_id': 'agg_sex_id'})
        validate = 'many_to_one'
    dist = dist.loc[dist['measure_id'] == 6, ]
    df = pd.merge(
        dist, cfr, how='outer', on=id_cols + ['pathogen'],
        validate=validate, indicator=True,
        suffixes=('_dist', '_cfr')
    )
    assert (df._merge != 'right_only').all()
    df = pd.merge(
        df, cfr.query("pathogen == 'other'"), how='left',
        on=id_cols, validate='many_to_one'
    )
    cfr_cols = [d + '_cfr' for d in draw_cols]
    dist_cols = [d + '_dist' for d in draw_cols]
    df[cfr_cols] = df[cfr_cols].fillna(
        df[draw_cols].rename(columns=lambda x: x + '_cfr')
    )
    assert df[cfr_cols].notnull().values.all()
    df[draw_cols] = pd.DataFrame(
        df[cfr_cols].to_numpy() * df[dist_cols].to_numpy(),
        index=df.index
    )
    assert df[draw_cols].notnull().values.all()
    if 'agg_age_group_id' in id_cols:
        id_cols.remove('agg_age_group_id')
        id_cols += ['age_group_id']
    if 'agg_sex_id' in id_cols:
        id_cols.remove('agg_sex_id')
        id_cols += ['sex_id']
    df = df.groupby(id_cols, as_index=False)[draw_cols].sum()
    assert (df[draw_cols] <= 1).values.all()
    assert (df[draw_cols] >= 0).values.all()
    df[draw_cols] = 1 / df[draw_cols]
    assert df[draw_cols].notnull().values.all()
    return df


def main(num_workers, year_id, infectious_syndrome):
    dfs = []
    if infectious_syndrome in PROXY_CAUSES:
        dfs.append(
            get_mi_ratios(num_workers, year_id, infectious_syndrome).assign(
                infectious_syndrome=infectious_syndrome
            )
        )
    if infectious_syndrome in CFR_SYNS:
        dfs.append(get_cfr(year_id, infectious_syndrome).assign(
            infectious_syndrome=infectious_syndrome
        ))
    df = pd.concat(dfs, sort=False)
    AmrResult(
        process='calculate_mi_ratios',
        burden_type='nonfatal',
        year_id=year_id,
        infectious_syndrome=infectious_syndrome
    ).write_results(df)
    print_log_message("Done!")


if __name__ == '__main__':
    num_workers = int(sys.argv[1])
    task_id = os.environ.get('SGE_TASK_ID')
    if task_id:
        print(f'Running as array job, task_id: {task_id}')
        task_row = pd.read_csv("FILEPATH").iloc[int(task_id) - 1]
        year_id = int(task_row['year_id'])
        infectious_syndrome = str(task_row['infectious_syndrome'])
    else:
        year_id = int(sys.argv[2])
        infectious_syndrome = str(sys.argv[3])
    main(num_workers, year_id, infectious_syndrome)
