import pandas as pd
import sys
import os
import numpy as np
from cod_prep.claude.configurator import Configurator
from cod_prep.utils import print_log_message, report_if_merge_fail, report_duplicates
from cod_prep.downloaders import (
    getcache_age_aggregate_to_detail_map,
    get_current_location_hierarchy,
    prep_child_to_available_parent_loc,
    get_current_cause_hierarchy,
    add_age_metadata
)
from amr_prep.utils.amr_io import AmrResult
from amr_prep.utils.misc import get_prepped_csbg_universe

CONF = Configurator()
CSV_ID = CONF.get_id('computation_cause_set_version')
CS_ID = CONF.get_id('computation_cause_set')


class PathogenSplitter(object):

    def __init__(self, burden_type, year_id, cause_id, infectious_syndrome):
        self.dem_cols = [
            'location_id', 'sex_id', 'year_id', 'age_group_id', 'hosp'
        ]
        self.draw_cols = ["draw_" + str(x) for x in range(0, 1000)]
        self.cache_kwargs = {'force_rerun': False, 'block_rerun': True, 'cache_results': False}
        self.conf = Configurator()
        self.burden_type = burden_type
        self.year_id = year_id
        self.cause_id = cause_id
        self.infectious_syndrome = infectious_syndrome

    def convert_to_incidence(self, df, mir):
        print_log_message("Prepping for merge")
        merge_cols = [
            'location_id', 'year_id', 'age_group_id', 'sex_id',
            'infectious_syndrome'
        ]
        special_syndromes = self.conf.get_id("syndromes_w_one_path")
        if self.infectious_syndrome in special_syndromes:
            merge_cols += ['cause_id']
        elif 'cause_id' in mir:
            mir = mir.drop('cause_id', axis='columns')
        if 'hosp' in mir:
            if mir['hosp'].isna().values.all():
                mir.drop(columns=['hosp'], inplace=True)
            else:
                merge_cols += ['hosp']
        if not (set(df.age_group_id) <= set(mir.age_group_id)):
            age_map = getcache_age_aggregate_to_detail_map(**self.cache_kwargs)
            age_map = age_map.loc[age_map.agg_age_group_id.isin(mir.age_group_id.unique())]
            df = df.merge(
                age_map, how='left', on='age_group_id', validate='many_to_one'
            )
            report_if_merge_fail(df, 'agg_age_group_id', 'age_group_id')
            merge_cols.remove('age_group_id')
            merge_cols.append('agg_age_group_id')
            mir = mir.rename(columns={'age_group_id': 'agg_age_group_id'})
        if not (set(df.sex_id) <= set(mir.sex_id)):
            df['parent_sex_id'] = 3
            merge_cols.remove('sex_id')
            merge_cols.append('parent_sex_id')
            mir = mir.rename(columns={'sex_id': 'parent_sex_id'})

        print_log_message("Merging")
        df = df.merge(
            mir, how='outer', on=merge_cols, validate='many_to_one',
            indicator=True, suffixes=('', '_ratio')
        )
        assert (df._merge != 'left_only').all()
        print_log_message("Calculating incident cases")
        df[self.draw_cols] = pd.DataFrame(
            df[self.draw_cols].to_numpy() * df.filter(regex='draw_\\d{1,3}_ratio').to_numpy(),
            index=df.index
        )
        if self.infectious_syndrome in special_syndromes:
            df = df.loc[(df._merge != 'right_only') | df.incidence_only]
            assert set(df.loc[df.hosp.notnull(), 'hosp']) == {'all'}
            df['hosp'] = df['hosp'].fillna('all')
            df.loc[df.incidence_only, self.draw_cols] = df.filter(
                regex='draw_\\d{1,3}_incidence'
            ).rename(columns=lambda x: x.replace("_incidence", ""))
        else:
            df = df.loc[df._merge != 'right_only']
        df = df[
            self.dem_cols + self.draw_cols + ['infectious_syndrome', 'cause_id']
        ]
        assert df.notnull().values.all()
        assert (df != np.inf).values.all()
        return df

    @staticmethod
    def read_pathogen_models(conf, infectious_syndrome, year_id, **cache_kwargs):
        model_versions = pd.read_csv("FILEPATH")
        model_versions = model_versions.groupby('infectious_syndrome').apply(
            lambda d: dict(zip(d['model_version'], d['query_string']))).to_dict()
        if infectious_syndrome in ['cardiac_infectious']:
            pull_syndrome = 'blood_stream_infectious'
        else:
            pull_syndrome = infectious_syndrome

        bug = pd.DataFrame()
        for mv, query_string in model_versions[pull_syndrome].items():
            model_df = pd.read_csv("FILEPATH")
            if type(query_string) == str:
                model_df = add_age_metadata(
                    model_df, ['age_group_days_start', 'age_group_days_end'],
                    **cache_kwargs
                )
                report_if_merge_fail(model_df, 'age_group_days_start', 'age_group_id')
                model_df = model_df.query(query_string)
                model_df = model_df.drop(
                    ['age_group_days_start', 'age_group_days_end'],
                    axis='columns'
                )
            bug = bug.append(model_df)
        bug['infectious_syndrome'] = infectious_syndrome
        return bug

    def split_pathogen(self):
        print_log_message(f"Reading pathogen draws")
        bug = self.read_pathogen_models(
            self.conf, self.infectious_syndrome, self.year_id, **self.cache_kwargs)

        print_log_message("Multiplying by pathogen distribution")
        if self.burden_type == 'fatal':
            bug = bug.loc[bug.measure_id == 1]
        elif self.burden_type == 'nonfatal':
            bug = bug.loc[bug.measure_id == 6]

        if (self.infectious_syndrome == 'diarrhea'):
            pathogen_merge_cols = [
                'infectious_syndrome', 'location_id',
                'age_group_id', 'sex_id', 'year_id'
            ]
        else:
            age_map = getcache_age_aggregate_to_detail_map(**self.cache_kwargs)
            age_map = age_map.loc[age_map.agg_age_group_id.isin(bug.age_group_id.unique())]
            self.df = self.df.merge(age_map, how='left', on='age_group_id', validate='many_to_one')
            report_if_merge_fail(self.df, 'agg_age_group_id', 'age_group_id')
            self.df['parent_sex_id'] = 3
            lh = get_current_location_hierarchy(
                location_set_version_id=self.conf.get_id("location_set_version"),
                **self.cache_kwargs)
            loc_map = prep_child_to_available_parent_loc(
                lh, bug.location_id.unique().tolist(), min_level=3)
            self.df = self.df.merge(loc_map, how='left', on='location_id', validate='many_to_one')
            report_if_merge_fail(self.df, 'parent_location_id', 'location_id')
            bug = bug.rename(
                columns={
                    'age_group_id': 'agg_age_group_id',
                    'location_id': 'parent_location_id',
                    'sex_id': 'parent_sex_id'
                }
            )
            pathogen_merge_cols = [
                'agg_age_group_id', 'parent_sex_id', 'parent_location_id',
                'year_id', 'infectious_syndrome', 'hosp'
            ]

        assert np.allclose(bug.groupby(pathogen_merge_cols)[self.draw_cols].sum(), 1)
        report_duplicates(bug, pathogen_merge_cols + ['pathogen'])
        self.df = self.df.merge(
            bug, how='left',
            on=pathogen_merge_cols,
            suffixes=('', '_path')
        )

        report_if_merge_fail(self.df, self.draw_cols[0] + '_path', pathogen_merge_cols)
        self.df[self.draw_cols] = pd.DataFrame(
            self.df[self.draw_cols].to_numpy() * self.df.filter(
                regex='draw_\\d{1,3}_path').to_numpy(),
            index=self.df.index
        )
        self.df = self.df.drop([d + '_path' for d in self.draw_cols], axis='columns')

    def run_split(self):
        print_log_message(f"Reading syndrome draws")
        self.df = AmrResult(
            process='split_sepsis_syndrome',
            burden_type='fatal',
            year_id=self.year_id,
            cause_id=self.cause_id,
            infectious_syndrome=self.infectious_syndrome
        ).read_results()

        if self.burden_type == 'nonfatal':
            print_log_message("Reading MI ratios")
            mir = AmrResult(
                process='calculate_mi_ratios',
                burden_type=self.burden_type,
                year_id=self.year_id,
                infectious_syndrome=self.infectious_syndrome
            ).read_results()
            if self.infectious_syndrome == 'others_and_non_bacterial_infectious':
                self.df = self.convert_to_incidence(self.df, mir)
            elif self.infectious_syndrome == 'blood_stream_infectious':
                if self.cause_id == 368:
                    mir = mir.loc[mir.cause_id == 368]
                    self.df = self.convert_to_incidence(self.df, mir)
                else:
                    self.df = add_age_metadata(
                        self.df, 'age_group_days_end', **self.cache_kwargs
                    )
                    report_if_merge_fail(self.df, 'age_group_days_end', 'age_group_id')
                    neonatal = self.df.query("age_group_days_end <= 27")
                    non_neonatal = self.df.query("age_group_days_end > 27")
                    neonatal = self.convert_to_incidence(
                        neonatal, mir.loc[mir.cause_id == 383].copy()
                    )
                    non_neonatal = self.convert_to_incidence(
                        non_neonatal, mir.loc[mir.cause_id.isnull()].copy()
                    )
                    self.df = pd.concat([neonatal, non_neonatal], sort=False)
            elif self.infectious_syndrome in \
                    ['peritoneal_and_intra_abdomen_infectious', 'bone_joint_infection']:
                self.df = self.convert_to_incidence(self.df, mir)
            else:
                synd_cause_mir = {
                    'cns_infectious': ['meningitis'],
                    'cardiac_infectious': ['cvd_endo'],
                    'respiratory_infectious': ['lri'],
                    'skin_infectious': ['skin_bacterial', 'skin_cellulitis', 'skin_decubitis'],
                    'typhoid_paratyphoid_ints': [
                        'intest_typhoid', 'intest_paratyph', 'intest_ints'],
                    'diarrhea': ['diarrhea'],
                    'chlamydia_and_gonorrheae': ['std_chlamydia', 'std_gonnorhea'],
                    'uti_plus': ['urinary_nephritis'],
                    'tb': ['tb_other', 'tb_mdr', 'tb_xdr'],
                }
                synd_acauses = synd_cause_mir.get(self.infectious_syndrome)
                ch = get_current_cause_hierarchy(
                    cause_set_id=CS_ID, cause_set_version_id=CSV_ID,
                    **self.cache_kwargs
                )
                synd_cause_ids = ch.loc[ch['acause'].isin(synd_acauses), 'cause_id'].tolist()
                if self.cause_id in synd_cause_ids:
                    mir = mir.loc[mir.cause_id == self.cause_id, ]
                    self.df = self.convert_to_incidence(self.df, mir)
                else:
                    mir = mir.loc[mir.cause_id.isnull()]
                    self.df = self.convert_to_incidence(self.df, mir)

        if self.infectious_syndrome in self.conf.get_id("syndromes_w_one_path"):
            print_log_message("This syndrome/cause is 1 pathogen")
            csbg_universe = get_prepped_csbg_universe(
                burden_type=self.burden_type, add_parents=False
            )
            pathogen = csbg_universe.query(
                f"infectious_syndrome == '{self.infectious_syndrome}' & "
                f"cause_id == {self.cause_id}"
            )[['pathogen']].drop_duplicates()
            assert len(pathogen) == 1
            pathogen = pathogen['pathogen'].unique()[0]
            self.df['pathogen'] = pathogen
        elif self.infectious_syndrome == 'others_and_non_bacterial_infectious':
            self.df['pathogen'] = '(none_estimated)'
        else:
            self.split_pathogen()

        self.df = self.df[
            self.dem_cols + ['cause_id', 'infectious_syndrome', 'pathogen']
            + self.draw_cols
        ]
        assert self.df.notnull().values.all()

    def save_output(self):
        print_log_message("Saving outputs")
        for pathogen, pathogen_df in self.df.groupby('pathogen'):
            print_log_message(f"Working on pathogen {pathogen}")
            AmrResult(
                process='split_pathogen',
                burden_type=self.burden_type,
                year_id=self.year_id,
                cause_id=self.cause_id,
                infectious_syndrome=self.infectious_syndrome,
                pathogen=pathogen
            ).write_results(pathogen_df)


if __name__ == '__main__':
    burden_type = str(sys.argv[1])

    task_id = os.environ.get('SGE_TASK_ID')
    if task_id:
        print(f'Running as array job, task_id: {task_id}')
        task_row = pd.read_csv("FILEPATH").iloc[int(task_id) - 1]
        year_id = int(task_row['year_id'])
        cause_id = int(task_row['cause_id'])
        infectious_syndrome = str(task_row['infectious_syndrome'])
    else:
        year_id = int(sys.argv[2])
        cause_id = int(sys.argv[3])
        infectious_syndrome = str(sys.argv[4])

    splitter = PathogenSplitter(
        burden_type=burden_type,
        year_id=year_id, cause_id=cause_id,
        infectious_syndrome=infectious_syndrome
    )
    splitter.run_split()
    splitter.save_output()
