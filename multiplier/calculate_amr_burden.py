import pandas as pd
import sys
import os
from cod_prep.claude.configurator import Configurator
from amr_prep.utils.amr_io import AmrResult
from cod_prep.utils import (
    print_log_message,
    report_if_merge_fail,
    report_duplicates
)
from cod_prep.downloaders.ages import add_age_metadata
from amr_prep.utils.misc import get_pred_ex, get_prepped_csbg_universe

CONF = Configurator()


class AmrBurdenCalculator(object):

    def __init__(self, burden_type, year_id, cause_id, infectious_syndrome,
                 pathogen):
        self.burden_type = burden_type
        self.year_id = year_id
        self.cause_id = cause_id
        self.infectious_syndrome = infectious_syndrome
        self.pathogen = pathogen
        self.conf = Configurator()
        self.cache_kwargs = {'force_rerun': False, 'block_rerun': True, 'cache_results': False}
        self.dem_cols = ['location_id', 'age_group_id', 'sex_id', 'year_id']
        self.draw_cols = ["draw_" + str(x) for x in range(0, 1000)]
        self.csbg = get_prepped_csbg_universe(burden_type=self.burden_type)

    def calculate_ylls(self, df):
        df = df.copy()
        pred_ex = get_pred_ex(codcorrect_version=self.conf.get_id("codcorrect_version"))
        pred_ex_id_cols = ['location_id', 'age_group_id', 'sex_id', 'year_id']
        df = df.merge(
            pred_ex, how='left', on=pred_ex_id_cols,
            validate='many_to_one'
        )
        report_if_merge_fail(df, 'pred_ex', pred_ex_id_cols)
        draw_cols = self.draw_cols
        df[draw_cols] = pd.DataFrame(
            df.filter(regex='draw_\\d{1,3}$').to_numpy()
            * df[['pred_ex']].to_numpy(),
            index=df.index
        )
        df = df.drop('pred_ex', axis='columns')
        return df

    def calculate_ylds(self, df, ypc):
        df = df.copy()
        merge_cols = [
            'location_id', 'year_id', 'age_group_id', 'sex_id'
        ]
        special_syndromes = self.conf.get_id("syndromes_w_one_path")
        if self.infectious_syndrome in special_syndromes:
            merge_cols += ['cause_id']
        elif self.infectious_syndrome != 'others_and_non_bacterial_infectious':
            ypc = ypc.drop('cause_id', axis='columns')

        df = df.merge(
            ypc, how='left', on=merge_cols, validate='many_to_one',
            indicator=True, suffixes=('', '_ypc')
        )
        assert (df._merge == 'both').all()
        df[self.draw_cols] = pd.DataFrame(
            df[self.draw_cols].to_numpy() * df.filter(regex='draw_\\d{1,3}_ypc').to_numpy(),
            index=df.index
        )
        assert df[self.draw_cols].notnull().values.all()
        df = df.drop([f'draw_{x}_ypc' for x in range(0, 1000)] + ['_merge'], axis='columns')
        return df

    def multiply_prop(self, df, prop):
        assert (prop.age_group_id == 22).all()
        df['agg_age_group_id'] = 22
        df['parent_sex_id'] = 3
        prop = prop.rename(columns={'age_group_id': 'agg_age_group_id', 'sex_id': 'parent_sex_id'})
        prop_merge_cols = [
            'agg_age_group_id', 'parent_sex_id', 'location_id', 'year_id', 'pathogen',
            'measure_id'
        ]
        if self.pathogen == 'mycobacterium_tuberculosis':
            cause_to_abx = {
                934: ['isoniazid', 'rifampicin', 'all_susceptible'],
                946: ['mdr'],
                947: ['xdr']
            }
            prop = prop.loc[prop.abx_class.isin(cause_to_abx[self.cause_id])]
        report_duplicates(prop, prop_merge_cols + ['counterfactual', 'abx_set', 'abx_class'])
        assert set(prop.location_id) == set(df.location_id)
        df = df.set_index(sorted([c for c in df if 'draw' not in c]))
        prop = prop.set_index(sorted([c for c in prop if 'draw' not in c]))
        df = df.mul(prop).reset_index(drop=False)
        return df

    def run_calculator(self):
        print_log_message(
            f"Reading in {self.burden_type} counts for pathogen {self.pathogen}, "
            f"year {self.year_id}, cause {self.cause_id}, infectious_syndrome "
            f"{self.infectious_syndrome}"
        )
        df = AmrResult(
            process='split_pathogen',
            burden_type=self.burden_type,
            year_id=self.year_id,
            cause_id=self.cause_id,
            infectious_syndrome=self.infectious_syndrome,
            pathogen=self.pathogen
        ).read_results()

        print_log_message("Converting to year-based metrics")
        if self.burden_type == 'fatal':
            df = df.assign(measure_id=1).append(
                self.calculate_ylls(df).assign(measure_id=4)
            )
        elif self.burden_type == 'nonfatal':
            ypc = AmrResult(
                process='calculate_ylds_per_case', burden_type=self.burden_type,
                year_id=self.year_id, infectious_syndrome=self.infectious_syndrome
            ).read_results()
            if self.infectious_syndrome == 'blood_stream_infectious':
                df = add_age_metadata(
                    df, 'age_group_days_end', **self.cache_kwargs
                )
                report_if_merge_fail(df, 'age_group_days_end', 'age_group_id')
                neonatal = df.query("age_group_days_end <= 27")
                non_neonatal = df.query("age_group_days_end > 27")
                neonatal = self.calculate_ylds(
                    neonatal, ypc.loc[ypc.cause_id == 383].copy()
                )
                non_neonatal = self.calculate_ylds(
                    non_neonatal, ypc.loc[ypc.cause_id == 368].copy()
                )
                df = df.assign(measure_id=6).append(
                    pd.concat([neonatal, non_neonatal], sort=False)
                    .assign(measure_id=3)
                )
            else:
                if self.infectious_syndrome == 'skin_infectious':
                    if self.cause_id in [656, 657, 665]:
                        ypc = ypc.query(f"cause_id == {self.cause_id}")
                    else:
                        ypc = ypc.query(f"cause_id == 980")
                df = df.assign(measure_id=6).append(
                    self.calculate_ylds(df, ypc).assign(measure_id=3)
                )

        print_log_message("Calculating AMR burden for 2 counterfactuals")
        abxs = self.csbg.loc[(self.csbg['pathogen'] == self.pathogen),
                             'abx_class'].unique().tolist()
        if ('none_tested' not in abxs):
            props = AmrResult(
                process='calculate_amr_props',
                burden_type=self.burden_type,
                year_id=self.year_id,
                pathogen=self.pathogen
            ).read_results()
            print_log_message("Got the amr props, multiplying them to split pathogen results")
            df = self.multiply_prop(df, props)
        else:
            print_log_message("No AMR burden for this pathogen, report all pathogen deaths as is")
            df['abx_set'] = 'all'
            df['abx_class'] = '(none_tested)'
            df['counterfactual'] = 'no_infection'

        df = df[
            self.dem_cols + [
                'cause_id', 'infectious_syndrome',
                'pathogen', 'abx_set', 'abx_class', 'hosp',
                'measure_id', 'counterfactual'
            ] + self.draw_cols
        ]
        assert df.notnull().values.all()

        print_log_message("Saving output")
        AmrResult(
            'calculate_amr_burden',
            self.burden_type,
            self.year_id,
            self.cause_id,
            self.infectious_syndrome,
            self.pathogen
        ).write_results(df)
        return df


if __name__ == '__main__':
    burden_type = str(sys.argv[1])
    task_id = os.environ.get('SGE_TASK_ID')
    if task_id:
        print(f'Running as array job, task_id: {task_id}')
        task_row = pd.read_csv("FILEPATH").iloc[int(task_id) - 1]
        year_id = int(task_row['year_id'])
        cause_id = int(task_row['cause_id'])
        infectious_syndrome = str(task_row['infectious_syndrome'])
        pathogen = str(task_row['pathogen'])
    else:
        year_id = int(sys.argv[2])
        cause_id = int(sys.argv[3])
        infectious_syndrome = str(sys.argv[4])
        pathogen = str(sys.argv[5])
    calc = AmrBurdenCalculator(
        burden_type=burden_type, year_id=year_id, cause_id=cause_id,
        infectious_syndrome=infectious_syndrome, pathogen=pathogen
    )
    calc.run_calculator()
