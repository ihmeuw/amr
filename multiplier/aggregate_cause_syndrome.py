import pandas as pd
import os
import sys
from cod_prep.downloaders import (
    get_all_related_causes,
    get_current_cause_hierarchy
)
from cod_prep.claude.configurator import Configurator
from amr_prep.utils.amr_io import read_results_wrapper, AmrResult
from cod_prep.utils import print_log_message
from amr_prep.utils.misc import get_prepped_csbg_universe
from itertools import repeat
CONF = Configurator()


class CauseSyndromeAggregator(object):

    def __init__(self, burden_type, year_id, cause_id, infectious_syndrome,
                 pathogen, num_workers=1):
        self.burden_type = burden_type
        self.year_id = year_id
        self.cause_id = cause_id
        self.infectious_syndrome = infectious_syndrome
        self.pathogen = pathogen
        self.num_workers = num_workers

        self.conf = Configurator()
        self.cache_kwargs = {'force_rerun': False, 'block_rerun': True, 'cache_results': False}
        self.id_cols = [
            'location_id', 'age_group_id', 'sex_id', 'year_id',
            'cause_id', 'infectious_syndrome', 'pathogen', 'abx_set',
            'abx_class', 'hosp', 'measure_id', 'counterfactual'
        ]
        self.draw_cols = ["draw_" + str(x) for x in range(0, 1000)]

        self.chh = get_current_cause_hierarchy(
            cause_set_version_id=self.conf.get_id("computation_cause_set_version"),
            cause_set_id=self.conf.get_id("computation_cause_set"),
            **self.cache_kwargs
        ).query("yld_only != 1")
        if self.burden_type == 'nonfatal':
            self.chh = self.chh.query("yll_only != 1")
        self.csbg_universe = get_prepped_csbg_universe(
            burden_type=self.burden_type,
            add_parents=True,
            **self.cache_kwargs
        )

    def get_child_causes(self):
        all_related = get_all_related_causes(self.cause_id, self.chh)
        detailed_causes = self.chh.query("most_detailed == 1").cause_id.unique()
        self.child_causes = list(set(all_related) & set(detailed_causes))

        if self.infectious_syndrome in self.conf.get_id("syndromes_w_one_path"):
            assoc_causes = self.csbg_universe.loc[
                (self.csbg_universe.infectious_syndrome == self.infectious_syndrome) & (
                    self.csbg_universe.pathogen == self.pathogen),
                'cause_id'
            ].unique().tolist()
            self.child_causes = list(set(assoc_causes) & set(self.child_causes))

    def get_child_syndromes(self):
        if self.infectious_syndrome != 'all':
            self.child_syndromes = [self.infectious_syndrome]
        else:
            self.child_syndromes = self.csbg_universe.loc[
                (self.csbg_universe.cause_id == self.cause_id) & (
                    self.csbg_universe.pathogen == self.pathogen) & (
                        self.csbg_universe.infectious_syndrome != 'all'),
                'infectious_syndrome'
            ].unique().tolist()
            assert len(self.child_syndromes) > 0

    def aggregate(self):
        self.get_child_causes()
        self.get_child_syndromes()
        if (self.child_causes != [self.cause_id]) & (self.infectious_syndrome == 'all'):
            self.child_causes = [self.cause_id]
        assert ((self.child_causes == [self.cause_id]) ^ (
            self.child_syndromes == [self.infectious_syndrome]))

        if self.child_syndromes == [self.infectious_syndrome]:
            arguments = list(zip(
                repeat(self.burden_type),
                self.child_causes,
                repeat(self.year_id),
                repeat(self.infectious_syndrome),
                repeat(self.pathogen)
            ))
            pieces = self.child_causes
        elif self.child_causes == [self.cause_id]:
            if pathogen == 'escherichia_coli':
                self.child_syndromes.remove('diarrhea')
            arguments = list(zip(
                repeat(self.burden_type),
                repeat(self.cause_id),
                repeat(self.year_id),
                self.child_syndromes,
                repeat(self.pathogen)
            ))
            pieces = self.child_syndromes

        print_log_message(f"Working on files {pieces}")
        chunks = [
            arguments[i:i + 20]
            for i in range(0, len(arguments), 20)
        ]
        df = pd.DataFrame()
        for chunk in chunks:
            df_list = []
            for arg in chunk:
                df_list.append(read_results_wrapper(*arg))
            print_log_message("Concatenating and grouping")
            df = df.append(
                pd.concat(df_list, axis='index', ignore_index=True, sort=True)
            )
            assert df.notnull().values.all()
            df = df.assign(
                cause_id=self.cause_id,
                infectious_syndrome=self.infectious_syndrome
            ).groupby(self.id_cols, as_index=False)[self.draw_cols].sum()

        print_log_message("Writing output!")
        AmrResult(
            process='aggregate_cause_syndrome',
            burden_type=self.burden_type,
            year_id=self.year_id,
            cause_id=self.cause_id,
            infectious_syndrome=self.infectious_syndrome,
            pathogen=self.pathogen,
        ).write_results(df)

        print_log_message("Finished!")


if __name__ == '__main__':
    burden_type = str(sys.argv[1])
    num_workers = int(sys.argv[2])
    status = int(sys.argv[3])
    task_id = os.environ.get('SGE_TASK_ID')
    if task_id:
        print(f'Running as array job, task_id: {task_id}')
        task_row = pd.read_csv("FILEPATH").iloc[int(task_id) - 1]
        year_id = int(task_row['year_id'])
        cause_id = int(task_row['cause_id'])
        infectious_syndrome = str(task_row['infectious_syndrome'])
        pathogen = str(task_row['pathogen'])
    else:
        year_id = int(sys.argv[4])
        cause_id = int(sys.argv[5])
        infectious_syndrome = str(sys.argv[6])
        pathogen = str(sys.argv[7])
    agg = CauseSyndromeAggregator(
        burden_type=burden_type, year_id=year_id, cause_id=cause_id,
        infectious_syndrome=infectious_syndrome, pathogen=pathogen,
        num_workers=num_workers
    )
    agg.aggregate()
