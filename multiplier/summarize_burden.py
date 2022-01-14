import pandas as pd
import sys
import os
import numpy as np
from mcod_prep.age_loc_aggregation import aggregate_locs
from cod_prep.downloaders import (
    get_current_location_hierarchy,
    getcache_age_aggregate_to_detail_map,
    add_population,
    get_pop,
    get_cod_ages,
    get_age_weights
)
from cod_prep.claude.configurator import Configurator
from amr_prep.utils.amr_io import AmrResult
from cod_prep.utils import print_log_message, report_duplicates, report_if_merge_fail
from mcod_prep.compile_burden import summarize_draws
from get_draws.api import get_draws
from amr_prep.utils.misc import get_prepped_csbg_universe

CONF = Configurator()


class AmrSummarizer(object):

    def __init__(self, burden_type, year_id, cause_id, infectious_syndrome,
                 pathogen, num_workers):
        self.burden_type = burden_type
        self.pathogen = pathogen
        self.year_id = year_id
        self.cause_id = cause_id
        self.infectious_syndrome = infectious_syndrome
        self.num_workers = num_workers
        self.num_draws = 1000
        self.id_cols = [
            'location_id', 'age_group_id', 'sex_id', 'year_id',
            'cause_id', 'infectious_syndrome', 'pathogen', 'abx_set',
            'abx_class', 'hosp', 'measure_id', 'counterfactual'
        ]
        self.draw_cols = ["draw_" + str(x) for x in range(0, 1000)]
        self.cache_kwargs = {'force_rerun': False, 'block_rerun': True, 'cache_results': False}
        self.conf = Configurator()
        self.csbg = get_prepped_csbg_universe(
            self.burden_type, add_parents=True, **self.cache_kwargs
        )

    def simple_aggregate(self, df, col_to_aggregate, detail_to_agg_map,
                         id_cols=None, value_cols=None):
        print_log_message(f"Aggregating {col_to_aggregate}")
        df = df.copy()
        id_cols = id_cols or self.id_cols
        value_cols = value_cols or self.draw_cols
        df = df.merge(detail_to_agg_map, how='inner', on=col_to_aggregate)
        df[col_to_aggregate] = df['agg_' + col_to_aggregate]
        assert df[col_to_aggregate].notnull().all()
        df = df.groupby(id_cols, as_index=False)[value_cols].sum()
        return df

    def aggregate_location(self, df):
        lhh = get_current_location_hierarchy(
            self.conf.get_id("location_set_version"),
            **self.cache_kwargs
        )
        return aggregate_locs(
            df, lhh, draw_cols=self.draw_cols, id_cols=self.id_cols,
            num_draws=len(self.draw_cols)
        )

    def aggregate_age(self, df, agg_age_group_ids, **agg_kwargs):
        age_map = getcache_age_aggregate_to_detail_map(**self.cache_kwargs)
        age_map = age_map.loc[age_map.agg_age_group_id.isin(agg_age_group_ids)]
        return self.simple_aggregate(df, 'age_group_id', age_map, **agg_kwargs)

    def aggregate_sex(self, df):
        sex_map = pd.DataFrame([[1, 3], [2, 3]], columns=['sex_id', 'agg_sex_id'])
        return self.simple_aggregate(df, 'sex_id', sex_map)

    def aggregate_hosp(self, df):
        hosp_map = pd.DataFrame(
            [['hospital', 'all'], ['community', 'all']],
            columns=['hosp', 'agg_hosp']
        )
        return self.simple_aggregate(df, 'hosp', hosp_map)

    def distinguish_mdr_xdr(self, df):
        pathogen_to_prefixes = {
            'salmonella_typhi': 'typh_paratyph_',
            'salmonella_paratyphi': 'typh_paratyph_',
            'mycobacterium_tuberculosis': 'tb_'
        }
        df = df.reset_index()
        df.loc[(df.abx_class.isin(['mdr', 'xdr'])), 'abx_class'] = df['pathogen'].map(
            pathogen_to_prefixes
        ).fillna('') + df['abx_class']
        return df

    def calculate_dalys(self, df):
        ylls = pd.concat([
            AmrResult(
                process='aggregate_cause_syndrome',
                burden_type='fatal',
                year_id=self.year_id,
                cause_id=self.cause_id,
                infectious_syndrome=self.infectious_syndrome,
                pathogen=pathogen
            ).read_results().query("measure_id == 4") for pathogen in self.read_pathogens
        ], sort=False)
        ylls = self.distinguish_mdr_xdr(ylls)
        ylls = ylls.assign(pathogen=self.pathogen)
        dalys = df.query("measure_id == 3").append(ylls)\
            .assign(measure_id=2)\
            .groupby(self.id_cols, as_index=False)[self.draw_cols].sum()
        return dalys

    def get_denominator(self, df):
        ages = df.age_group_id.unique().astype(int).tolist()
        locations = df.location_id.unique().astype(int).tolist()
        sexes = df.sex_id.unique().astype(int).tolist()

        get_draws_params = {
            'fatal': [{'source': 'codcorrect', 'measure_id': [1, 4], 'metric_id': [1]}],
            'nonfatal': [
                {'source': 'como', 'measure_id': [3, 6], 'metric_id': [3]},
                {'source': 'dalynator', 'measure_id': [2], 'metric_id': [1]}
            ]
        }
        ccs = []
        for run in get_draws_params[self.burden_type]:
            version_id = self.conf.get_id(f"{run['source']}_version")
            print_log_message(f"Reading draws from {run['source']}...")
            cc = get_draws(
                gbd_id_type=["cause_id"],
                gbd_id=[self.cause_id],
                source=run['source'],
                year_id=self.year_id,
                measure_id=run['measure_id'],
                metric_id=run['metric_id'],
                gbd_round_id=self.conf.get_id("gbd_round"),
                version_id=version_id,
                decomp_step=self.conf.get_id("decomp_step"),
                location_id=locations,
                sex_id=sexes,
                age_group_id=ages,
                num_workers=self.num_workers
            )
            if run['metric_id'] == [3]:
                print_log_message("Converting from rates to counts...")
                pop_cols = ['year_id', 'location_id', 'age_group_id', 'sex_id']
                cc = add_population(cc, pop_run_id=self.conf.get_id("pop_run"), **self.cache_kwargs)
                report_if_merge_fail(cc, 'population', pop_cols)
                cc.loc[:, self.draw_cols] = cc.loc[:, self.draw_cols].multiply(
                    cc.loc[:, 'population'], axis="index")
            ccs.append(cc)

        cc = pd.concat(ccs, sort=False)

        merge_cols = [
            'location_id', 'year_id', 'age_group_id', 'sex_id', 'cause_id',
            'measure_id'
        ]
        cc = cc[merge_cols + self.draw_cols]
        df = df.merge(
            cc, how='left', on=merge_cols,
            indicator=True, suffixes=('', '_denom'), validate='many_to_one'
        )
        assert (df._merge == 'both').all()
        df = df.drop(['_merge'], axis='columns')
        self.denom_cols = [f'draw_{x}_denom' for x in range(0, self.num_draws)]
        self.draw_cols += self.denom_cols
        return df

    def calculate_props(self, df):
        df = df.copy()
        self.draw_cols = ["draw_" + str(x) for x in range(0, 1000)]
        df[self.draw_cols] = pd.DataFrame(np.clip(
            np.true_divide(
                df[self.draw_cols].to_numpy(), df[self.denom_cols].to_numpy(),
                out=np.zeros_like(df[self.draw_cols]),
                where=df[self.denom_cols].to_numpy() != 0
            ), None, 1
        ), index=df.index)
        df = df.drop(self.denom_cols, axis='columns')
        return df

    def calculate_rates(self, df):
        df = df.copy()

        pop = get_pop(pop_run_id=self.conf.get_id("pop_run"), **self.cache_kwargs)
        need_agg_ages = list(set(df.age_group_id) - set(pop.age_group_id))
        if len(need_agg_ages) > 0:
            pop = pop.append(self.aggregate_age(
                pop, agg_age_group_ids=need_agg_ages,
                id_cols=['location_id', 'age_group_id', 'sex_id', 'year_id'],
                value_cols=['population']
            ))
        df = add_population(df, pop_df=pop)
        report_if_merge_fail(
            df, 'population', ['location_id', 'age_group_id', 'sex_id', 'year_id'])
        df[self.draw_cols] = pd.DataFrame(
            df[self.draw_cols].to_numpy() / df[['population']].to_numpy(),
            index=df.index
        )

        age_weights = get_age_weights(
            gbd_round_id=self.conf.get_id("gbd_round"),
            **self.cache_kwargs
        )
        detailed_ages = get_cod_ages(
            gbd_round_id=self.conf.get_id("gbd_round"),
            **self.cache_kwargs
        ).age_group_id.unique()
        age_weights = age_weights.loc[age_weights.age_group_id.isin(detailed_ages)]
        assert np.isclose(age_weights.age_group_weight_value.sum(), 1)
        age_std = pd.merge(
            df, age_weights, how='inner', on='age_group_id',
            validate='many_to_one'
        )
        age_std = age_std.reset_index(drop=True)
        age_std[self.draw_cols] = pd.DataFrame(
            age_std[self.draw_cols].to_numpy() * age_std[['age_group_weight_value']].to_numpy(),
            index=age_std.index
        )
        age_std = age_std.assign(age_group_id=27)
        age_std = age_std.groupby(self.id_cols)[self.draw_cols].sum().reset_index()

        df = df.drop(['population'], axis='columns')
        df = df.append(age_std, sort=False)
        return df

    def is_pathogen_aggregate(self):
        return self.csbg.query(
            f"infectious_syndrome == '{self.infectious_syndrome}' & "
            f"pathogen == '{self.pathogen}'"
        )['pathogen_aggregate'].max()

    def summarize(self):
        print_log_message(f"Reading in the draws files for pathogen {self.pathogen}")
        if self.is_pathogen_aggregate() == 0:
            self.read_pathogens = [self.pathogen]
        else:
            self.read_pathogens = self.csbg.loc[self.csbg.abx_class.notnull()].query(
                f"cause_id == {self.cause_id} & "
                f"infectious_syndrome == '{self.infectious_syndrome}' & "
                f"pathogen_aggregate == 0"
            ).pathogen.unique().tolist()
            if self.pathogen == 'escherichia_coli':
                self.read_pathogens = [
                    path for path in self.read_pathogens if 'escherichia_coli' in path]

        df = pd.concat([AmrResult(
            process='aggregate_cause_syndrome',
            burden_type=self.burden_type,
            year_id=self.year_id,
            cause_id=self.cause_id,
            infectious_syndrome=self.infectious_syndrome,
            pathogen=pathogen
        ).read_results() for pathogen in self.read_pathogens], sort=False)
        df = self.distinguish_mdr_xdr(df)
        df = df.assign(pathogen=self.pathogen)
        df = df.groupby(self.id_cols, as_index=False)[self.draw_cols].sum()

        if self.burden_type == 'nonfatal':
            print_log_message("Calculating DALYs")
            df = df.append(self.calculate_dalys(df))

        print_log_message("Aggregating drugs")
        if self.infectious_syndrome in ['tb', 'all']:
            all_resistant = df.loc[
                (df['counterfactual'] == 'no_infection') & (
                    df['abx_class'].isin(
                        ['tb_mdr', 'tb_xdr', 'isoniazid', 'rifampicin']
                    )
                )
            ].assign(abx_class='all_resistant')\
             .groupby(self.id_cols, as_index=False)[self.draw_cols].sum()
            df = df.append(all_resistant)
        resistant_and_susceptible = df.loc[
            (df['counterfactual'] == 'no_infection') & (
                df['abx_class'].isin(['all_resistant', 'all_susceptible']))] \
            .assign(abx_class='all_resistant_and_susceptible') \
            .groupby(self.id_cols, as_index=False)[self.draw_cols].sum()
        all_resistant = df.loc[df['counterfactual'] == 'susceptible_infection', ] \
            .assign(abx_class='all_resistant') \
            .groupby(self.id_cols, as_index=False)[self.draw_cols].sum()
        df = pd.concat([df, all_resistant, resistant_and_susceptible])
        if self.pathogen == 'all':
            all_total = df.loc[
                (df['counterfactual'] == 'no_infection')
                & (df['abx_class'].isin(['(none_tested)', 'all_resistant_and_susceptible']))
                & (df['pathogen'] == 'all'), ] \
                .assign(abx_class='all_total') \
                .groupby(self.id_cols, as_index=False)[self.draw_cols].sum()
            df = df.append(all_total)

        df = df.append(self.aggregate_hosp(df.loc[df.hosp != 'all']))
        df = df.groupby(self.id_cols, as_index=False)[self.draw_cols].sum()

        print_log_message("Aggregating location, age, and sex")
        df = df.append(self.aggregate_location(df))
        df = df.append(self.aggregate_age(df, agg_age_group_ids=[42, 192, 22]))
        df = df.append(self.aggregate_sex(df))
        report_duplicates(df, self.id_cols)

        print_log_message("Calculating rates")
        df = df.assign(metric_id=1)
        self.id_cols += ['metric_id']
        report_duplicates(df, self.id_cols)
        df = df.append(
            self.calculate_rates(df.loc[df.metric_id == 1]).assign(metric_id=3)
        )
        report_duplicates(df, self.id_cols)
        assert df.notnull().values.all()

        print_log_message("Saving results with draws...")
        AmrResult(
            process='summarize_burden',
            burden_type=self.burden_type,
            year_id=self.year_id, cause_id=self.cause_id,
            infectious_syndrome=self.infectious_syndrome,
            pathogen=self.pathogen,
            save_draws=True
        ).write_results(df)

        print_log_message("Summarizing draws...")
        df = summarize_draws(df, draw_cols=self.draw_cols, prefix='amr_')
        df = df.reset_index(drop=True)
        df.loc[
            df.counterfactual == 'susceptible_infection',
            ['amr_mean', 'amr_upper', 'amr_lower']
        ] = np.clip(df[['amr_mean', 'amr_upper', 'amr_lower']], 0, None)
        assert df.notnull().values.all()

        print_log_message("Saving results without draws...")
        AmrResult(
            process='summarize_burden',
            burden_type=self.burden_type,
            year_id=self.year_id, cause_id=self.cause_id,
            infectious_syndrome=self.infectious_syndrome,
            pathogen=self.pathogen
        ).write_results(df)
        return df


if __name__ == '__main__':
    burden_type = str(sys.argv[1])
    num_workers = int(sys.argv[2])
    task_id = os.environ.get('SGE_TASK_ID')
    if task_id:
        print(f'Running as array job, task_id: {task_id}')
        task_row = pd.read_csv("FILEPATH").iloc[int(task_id) - 1]
        year_id = int(task_row['year_id'])
        cause_id = int(task_row['cause_id'])
        infectious_syndrome = str(task_row['infectious_syndrome'])
        pathogen = str(task_row['pathogen'])
    else:
        year_id = int(sys.argv[3])
        cause_id = int(sys.argv[4])
        infectious_syndrome = str(sys.argv[5])
        pathogen = str(sys.argv[6])

    summarizer = AmrSummarizer(
        burden_type=burden_type,
        year_id=year_id,
        cause_id=cause_id,
        infectious_syndrome=infectious_syndrome,
        pathogen=pathogen,
        num_workers=num_workers
    )
    summarizer.summarize()
