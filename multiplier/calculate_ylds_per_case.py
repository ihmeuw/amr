import sys
import os
import numpy as np
import pandas as pd
from get_draws.api import get_draws
from cod_prep.downloaders import (
    get_current_location_hierarchy, get_cod_ages,
    get_current_cause_hierarchy
)
from cod_prep.utils import (
    print_log_message,
    report_duplicates
)
from cod_prep.claude.configurator import Configurator
from amr_prep.utils.amr_io import AmrResult
from amr_prep.utils.misc import drop_unmodeled_asc


CACHE_KWARGS = {
    'force_rerun': False,
    'block_rerun': True,
    'cache_results': False
}

CONF = Configurator('standard')
GBD_ROUND = CONF.get_id("gbd_round")
DECOMP_STEP = CONF.get_id("decomp_step")
LSV_ID = CONF.get_id('location_set_version')
CSV_ID = CONF.get_id('computation_cause_set_version')
CS_ID = CONF.get_id('computation_cause_set')
POP_RUN = CONF.get_id('pop_run')

cause_for_others = [
    'diptheria',
    'encephalitis',
    'hepatitis',
    'hiv_other',
    'infectious',
    'malaria',
    'measles',
    'ntd_afrtryp',
    'ntd_chagas',
    'ntd_dengue',
    'ntd_rabies',
    'ntd_yellowfever',
    'otitis',
    'skin_scabies',
    'std_herpes',
    'uri',
    'varicella'
]

PROXY_CAUSES = {
    'cns_infectious': ['meningitis'],
    'cardiac_infectious': ['cvd_endo'],
    'respiratory_infectious': ['lri'],
    'skin_infectious': ['skin_infect', 'skin_cellulitis', 'skin_bacterial', 'skin_decubitus'],
    'typhoid_paratyphoid_ints': ['intest_typhoid', 'intest_paratyph', 'intest_ints'],
    'diarrhea': ['diarrhea'],
    'chlamydia_and_gonorrheae': ['std_chlamydia', 'std_gonnorhea'],
    'uti_plus': ['urinary_nephritis'],
    'blood_stream_infectious': ['maternal_sepsis', 'neonatal_sepsis'],
    'tb': ['tb_other', 'tb_mdr', 'tb_xdr'],
    'peritoneal_and_intra_abdomen_infectious': ['digest_ileus'],
    'bone_joint_infection': ['skin_infect'],
    'others_and_non_bacterial_infectious': cause_for_others
}


def draws_wrapper(causes, year_id, measure_id, num_workers):
    locs = get_current_location_hierarchy(
        location_set_version_id=LSV_ID, **CACHE_KWARGS
    )
    all_cntry = list(locs.loc[locs['level'] == 3, 'location_id'].unique())
    ages = get_cod_ages(
        gbd_round_id=GBD_ROUND, **CACHE_KWARGS
    )['age_group_id'].unique().tolist()
    if 368 in causes:
        ages += [22]
    draws = get_draws(
        gbd_id_type=["cause_id"],
        gbd_id=causes,
        source='como',
        year_id=year_id,
        measure_id=measure_id,
        metric_id=3,  # Rates
        gbd_round_id=GBD_ROUND,
        version_id=CONF.get_id("como_version"),
        release_id=CONF.get_id("release"),
        location_id=all_cntry,
        sex_id=[1, 2],
        age_group_id=ages,
        num_workers=num_workers
    )
    return draws


def get_ylds_per_case(num_workers, year_id, infectious_syndrome):
    proxy_causes = PROXY_CAUSES[infectious_syndrome]
    ch = get_current_cause_hierarchy(
        cause_set_id=CS_ID, cause_set_version_id=CSV_ID,
        **CACHE_KWARGS
    ).set_index("acause")['cause_id'].to_dict()
    proxy_cause_ids = [ch[acause] for acause in proxy_causes]

    print_log_message(f"Pulling incidence for causes {proxy_cause_ids}")
    incd = draws_wrapper(
        causes=proxy_cause_ids,
        year_id=year_id,
        measure_id=6,  # Incidence
        num_workers=num_workers
    )

    print_log_message(f"Pulling YLDs for causes {proxy_cause_ids}")
    ylds = draws_wrapper(
        causes=proxy_cause_ids,
        year_id=year_id,
        measure_id=3,  # YLDs
        num_workers=num_workers
    )

    print_log_message("Calculating YLDs per case...")
    demo_cols = ['year_id', 'location_id', 'age_group_id', 'sex_id']
    draw_cols = [col for col in incd if 'draw_' in col]
    incd = incd[demo_cols + ['cause_id'] + draw_cols]
    ylds = ylds[demo_cols + ['cause_id'] + draw_cols]

    age_meta_df = get_cod_ages(gbd_round_id=GBD_ROUND, **CACHE_KWARGS)
    cause_meta_df = get_current_cause_hierarchy(
        cause_set_id=CS_ID, cause_set_version_id=CSV_ID,
        **CACHE_KWARGS
    )
    incd = drop_unmodeled_asc(incd, cause_meta_df, age_meta_df, 'nonfatal')
    ylds = drop_unmodeled_asc(ylds, cause_meta_df, age_meta_df, 'nonfatal')

    merge_cols = demo_cols + ['cause_id']
    df = pd.merge(
        ylds, incd, how='outer', on=merge_cols, validate='one_to_one',
        suffixes=('_yld', '_incd'), indicator=True
    )
    assert (df._merge == 'both').all()
    if infectious_syndrome == 'others_and_non_bacterial_infectious':
        yld_draw_cols = [f'draw_{i}_yld' for i in range(0, 1000)]
        incd_draw_cols = [f'draw_{i}_incd' for i in range(0, 1000)]
        df = df.groupby(demo_cols, as_index=False)[incd_draw_cols + yld_draw_cols].sum()

    yld_draws = df.filter(regex='.*_yld').to_numpy()
    incd_draws = df.filter(regex='.*_incd').to_numpy()
    df[draw_cols] = pd.DataFrame(
        np.true_divide(
            yld_draws, incd_draws, out=np.zeros_like(yld_draws),
            where=incd_draws != 0
        ), index=df.index
    )
    assert df.notnull().values.all()
    assert (df != np.inf).values.all()
    if infectious_syndrome != 'others_and_non_bacterial_infectious':
        df = df[demo_cols + ['cause_id'] + draw_cols]
    else:
        df = df[demo_cols + draw_cols]
    if infectious_syndrome == 'blood_stream_infectious':
        age_limits = tuple(
            cause_meta_df.query(f"cause_id == 368").iloc[0][
                ['yld_age_start', 'yld_age_end']
            ])
        add_ages = age_meta_df.query(
            f"simple_age < {age_limits[0]} | simple_age > {age_limits[1]}"
        ).age_group_id.unique().tolist()
        copy_df = df.query("age_group_id == 22 & cause_id == 368")
        df = df.append(pd.concat(
            [copy_df.assign(age_group_id=age) for age in add_ages]
        ))
        df = df.append(
            df.query("cause_id == 368 & sex_id == 2")
            .assign(sex_id=1)
        )
    if infectious_syndrome != 'others_and_non_bacterial_infectious':
        report_duplicates(df, demo_cols + ['cause_id'])
    else:
        report_duplicates(df, demo_cols)
    AmrResult(
        process='calculate_ylds_per_case',
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
    get_ylds_per_case(num_workers, year_id, infectious_syndrome)
