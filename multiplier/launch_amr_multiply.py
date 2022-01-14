import pandas as pd
import json
import requests
from pandas.io.json import json_normalize
from cod_prep.utils import print_log_message, cod_timestamp
from amr_prep.utils.amr_io import AmrResult
from amr_prep.utils.misc import get_prepped_csbg_universe
from mcod_prep.utils.mcod_cluster_tools import submit_mcod
from cod_prep.claude.configurator import Configurator
import argparse
from cod_prep.downloaders import (
    get_current_cause_hierarchy
)


class AmrMultiplier(object):

    def __init__(self, process, burden_type, year_id=None, cause_id=None,
                 infectious_syndrome=None, pathogen=None, abx_class=None,
                 run_params_table=None, queue='all.q', task_throttle=5000):
        self.process = process
        self.burden_type = burden_type
        assert self.burden_type in AmrResult.processes[self.process]['valid_burden_types'],\
            f"You cannot launch {self.process} for {self.burden_type} burden"
        self.year_id = self.coerce_to_list(year_id)
        self.cause_id = self.coerce_to_list(cause_id)
        self.infectious_syndrome = self.coerce_to_list(infectious_syndrome)
        self.pathogen = self.coerce_to_list(pathogen)
        self.abx_class = self.coerce_to_list(abx_class)
        self.run_params_table = run_params_table
        self.queue = queue
        self.task_throttle = task_throttle
        self.conf = Configurator()
        assert self.conf.config_type == 'amr'
        self.cache_kwargs = {'force_rerun': True, 'block_rerun': False}
        self.validate_and_set_run_parameters()

    def coerce_to_list(self, x):
        if x is not None and type(x) != list:
            x = [x]
        return x

    def validate_and_set_run_parameters(self):
        assert self.process in AmrResult.processes
        assert self.burden_type in AmrResult.burden_types

        process_to_valid_variables = {
            process: AmrResult.processes[process]['valid_params'] + ['year_id']
            for process in AmrResult.processes.keys()
        }
        process_to_valid_variables['split_sepsis_syndrome'] = list(
            set(process_to_valid_variables['split_sepsis_syndrome']) - {'infectious_syndrome'}
        )
        process_to_valid_variables['split_pathogen'] = list(
            set(process_to_valid_variables['split_pathogen']) - {'pathogen'}
        )

        run_params_dict = {
            'year_id': self.year_id,
            'cause_id': self.cause_id,
            'infectious_syndrome': self.infectious_syndrome,
            'pathogen': self.pathogen,
            'abx_class': self.abx_class
        }

        if self.run_params_table is not None:
            print_log_message(
                "You passed run_params_table, ignoring all other filters"
            )
            assert set(self.run_params_table) <= set(
                process_to_valid_variables[self.process])
            assert 'year_id' in self.run_params_table
            assert self.run_params_table.notnull().values.all()
        else:
            for param, value in run_params_dict.items():
                if param not in process_to_valid_variables[self.process]:
                    assert value is None, param + ' cannot be passed for ' + self.process + \
                        '. Please refer to process_to_valid_variables dictionary'

        ch = get_current_cause_hierarchy(
            cause_set_version_id=self.conf.get_id('computation_cause_set_version'),
            cause_set_id=self.conf.get_id('computation_cause_set'),
            **self.cache_kwargs
        )
        self.csbg_universe = get_prepped_csbg_universe(
            burden_type=self.burden_type,
            add_parents=self.process in ["aggregate_cause_syndrome", "summarize_burden"]
        )

        if self.process not in ["aggregate_cause_syndrome", "summarize_burden"]:
            self.csbg_universe = self.csbg_universe.loc[
                self.csbg_universe.cause_id.isin(
                    ch.query("most_detailed == 1").cause_id
                )
            ]
        elif self.process == 'aggregate_cause_syndrome':
            self.csbg_universe = self.csbg_universe.loc[
                self.csbg_universe.cause_id.isin(
                    ch.query("most_detailed != 1").cause_id
                ) | (self.csbg_universe.infectious_syndrome == 'all')
            ]
            self.csbg_universe = self.csbg_universe.loc[
                self.csbg_universe.pathogen_aggregate != 1
            ]

        tasks = self.csbg_universe.copy()
        if self.run_params_table is not None:
            tasks = tasks.merge(
                self.run_params_table, how='right',
                indicator=True
            )
            assert (tasks._merge == 'both').all(), f"You asked to run the "\
                f"following tasks which are not recognized for {self.process}"\
                f" : {tasks.loc[tasks._merge != 'both']}"
        else:
            for param, value in run_params_dict.items():
                if value is not None and param in tasks:
                    assert set(value) <= set(tasks[param]), f"You asked to run the "\
                        f"following values of {param} which are not "\
                        f"recognized for {self.process}: {set(value) - set(tasks[param])}"
            for param, value in run_params_dict.items():
                if value is not None and param in tasks:
                    tasks = tasks.loc[tasks[param].isin(value)]

        if 'year_id' not in tasks:
            assert self.year_id is not None
            tasks = pd.concat(
                [tasks.assign(year_id=yr) for yr in self.year_id],
                sort=False
            )

        if self.process in [
            'calculate_amr_props'
        ]:
            tasks = tasks.loc[tasks.abx_class != 'none_tested', ]

        tasks = tasks.loc[tasks[
            AmrResult.processes[self.process]['valid_params']
        ].notnull().all(axis=1)]
        tasks = tasks[process_to_valid_variables[self.process]].drop_duplicates()
        assert tasks.notnull().values.all()
        assert not tasks.duplicated().any()
        assert len(tasks) > 0, "Your filters returned no valid run parameters!"
        self.tasks = tasks

        self.expected_files = self.tasks.merge(
            self.csbg_universe[AmrResult.processes[self.process]['valid_params']]
            .drop_duplicates(),
            how='left'
        )
        assert self.expected_files.notnull().values.all()

    def set_launch_dirs(self):
        self.parent_dir = "FILEPATH"
        if self.process != 'aggregate_cause_syndrome':
            self.array_dir = "FILEPATH"
        else:
            self.array_dir = "FILEPATH"

    def launch(self):
        self.set_launch_dirs()
        process_resource_dict = {
            'split_sepsis_syndrome': {
                'num_workers': 10, 'memory': '20G', 'runtime': "00:30:00",
                'params': [10]
            },
            'calculate_mi_ratios': {
                'num_workers': 20, 'memory': '40G', 'runtime': "01:00:00",
                'params': [20]
            },
            'split_pathogen': {
                'num_workers': 1, 'memory': '20G', 'runtime': "24:00:00",
                'params': [self.burden_type]
            },
            'calculate_ylds_per_case': {
                'num_workers': 20, 'memory': '40G', 'runtime': "24:00:00",
                'params': [20]
            },
            'calculate_amr_props': {
                'num_workers': 5, 'memory': '50G', 'runtime': "24:00:00",
                'params': [self.burden_type]
            },
            'calculate_amr_burden': {
                'num_workers': 1, 'memory': '20G', 'runtime': "04:00:00",
                'params': [self.burden_type]
            },
            'aggregate_cause_syndrome': {
                'num_workers': 1, 'memory': '425G', 'runtime': "12:45:00",
                'params': [self.burden_type, 1]
            },
            'summarize_burden': {
                'num_workers': 10, 'memory': '475G', 'runtime': "15:30:00",
                'params': [self.burden_type, 10]
            }
        }
        if self.process != 'aggregate_cause_syndrome':
            print_log_message(f"Clearing {len(self.expected_files)} files")
            for index, row in self.expected_files.iterrows():
                AmrResult(
                    self.process, self.burden_type, **row.to_dict()
                ).clear_results()

            print_log_message(f"Submitting array with {len(self.tasks)} tasks")
            jobname = f"{self.process}_{self.burden_type}"
            worker = "FILEPATH"
            self.tasks['task_id'] = list(range(1, len(self.tasks) + 1))
            self.tasks.to_csv("FILEPATH")
            resources = process_resource_dict[self.process]
            jid = submit_mcod(
                jobname, 'python', worker, cores=resources['num_workers'],
                memory=resources['memory'], params=resources['params'],
                verbose=True, logging=True, runtime=resources['runtime'],
                queue=self.queue, num_tasks=len(self.tasks),
                task_throttle=self.task_throttle
            )
            with open("FILEPATH") as f:
                json.dump([int(jid)], f, indent=2)
        else:
            print_log_message(f"Clearing {len(self.expected_files)} files")
            for index, row in self.expected_files.iterrows():
                AmrResult(
                    self.process, self.burden_type, **row.to_dict()
                ).clear_results()
            print_log_message(f"Submitting array jobs with {len(self.tasks)} tasks")
            ch = get_current_cause_hierarchy(
                cause_set_version_id=self.conf.get_id('computation_cause_set_version'),
                cause_set_id=self.conf.get_id('computation_cause_set'),
                **self.cache_kwargs
            ).query("yld_only != 1")
            self.tasks['most_detailed'] = self.tasks['cause_id'].map(
                ch.set_index('cause_id')['most_detailed'].to_dict()
            )
            assert self.tasks.most_detailed.notnull().all()
            self.tasks['status'] = 1
            self.tasks.loc[
                (self.tasks.most_detailed != 1)
                & (self.tasks.infectious_syndrome == 'all'),
                'status'
            ] = 2
            holds = []
            for status in [1, 2]:
                tasks = self.tasks.query(f"status == {status}")
                if (status == 2) & (len(tasks) == 0):
                    print_log_message('No status 2 jobs to run for aggregation given parameters')
                    break
                jobname = f"{self.process}_{self.burden_type}_{status}"
                worker = "FILEPATH"
                tasks['task_id'] = list(range(1, len(tasks) + 1))
                tasks.to_csv("FILEPATH")
                resources = process_resource_dict[self.process]
                jid = submit_mcod(
                    jobname, 'python', worker, cores=resources['num_workers'],
                    memory=resources['memory'], params=resources['params'] + [status],
                    verbose=True, logging=True, runtime=resources['runtime'],
                    queue=self.queue, num_tasks=len(tasks), holds=holds,
                    task_throttle=self.task_throttle
                )
                holds += [jid]
            with open("FILEPATH") as f:
                json.dump([int(jid) for jid in holds], f, indent=2)

    def check_outputs_exist(self, save_path=None):
        self.expected_files['completed'] = self.expected_files.apply(
            lambda x: AmrResult(
                self.process, self.burden_type, **x.to_dict()
            ).results_path.exists(),
            axis=1
        )
        if self.process != 'aggregate_cause_syndrome':
            self.tasks['task_id'] = list(range(1, len(self.tasks) + 1))
        self.expected_files = self.expected_files.merge(
            self.tasks, how='left', validate='many_to_one')
        missing_files = self.expected_files.loc[~self.expected_files.completed]
        if len(missing_files) > 0:
            if save_path is not None:
                missing_files.to_excel("FILEPATH")
            print_log_message(
                f"Outputs for the following year/cause/pathogens do not exist: \n"
                f"{missing_files}"
            )
        else:
            print_log_message("All outputs exist!")

    def check_last_run(self, save_path=None):
        self.set_launch_dirs()
        if self.process != 'aggregate_cause_syndrome':
            tasks = pd.read_csv("FILEPATH")
            with open("FILEPATH") as f:
                jids = json.load(f)
            tasks['job_id'] = jids[0]
        else:
            tasks = [
                pd.read_csv("FILEPATH")
                for status in [1, 2]
            ]
            with open("FILEPATH") as f:
                jids = json.load(f)
            tasks[0]['job_id'] = jids[0]
            tasks[1]['job_id'] = jids[1]
            tasks = pd.concat(tasks, sort=False)
        qpid = pd.DataFrame()
        for job_id in jids:
            j = requests.get(
                "ADDRESS",
                params=[('job_number', job_id), ('limit', 40000)]
            ).json()
            qpid = qpid.append(json_normalize(j['data']))
        df = pd.merge(
            tasks, qpid, how='left', left_on=['job_id', 'task_id'],
            right_on=['job_number', 'task_number'],
            validate='one_to_one'
        )
        if (df.exit_status != 0).any():
            if save_path is not None:
                df.loc[df.exit_status != 0].to_excel("FILEPATH")
            print_log_message(
                f"Exit status other than zero for the following job/tasks: \n"
                f"{df.loc[df['exit_status'] != 0, list(tasks.columns) + ['exit_status']]}"
            )
        else:
            print_log_message("All jobs had exit status 0!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Launch various processes that together combine all'
        ' component models and estimate the fatal burden of AMR'
    )
    parser.add_argument('process', type=str, choices=list(AmrResult.processes.keys()))
    parser.add_argument('burden_type', type=str, choices=AmrResult.burden_types)
    parser.add_argument('--year_id', nargs='*', type=int, default=[2019])
    parser.add_argument('--cause_id', nargs='*', type=int)
    parser.add_argument('--infectious_syndrome', nargs='*', type=str)
    parser.add_argument('--pathogen', nargs='*', type=str)
    parser.add_argument('--abx_class', nargs='*', type=str)
    parser.add_argument('--queue', type=str, default='all.q')
    parser.add_argument('--task_throttle', type=int, default=5000)
    parser.add_argument('--check_outputs', action='store_true')
    parser.add_argument('--check_last_run', action='store_true')
    parser.add_argument('--failed_jobs_save_path', type=str)
    args = parser.parse_args()
    kwargs = vars(args)
    check_outputs = kwargs['check_outputs']
    check_last_run = kwargs['check_last_run']
    save_path = kwargs['failed_jobs_save_path']
    for kwarg in ['check_outputs', 'check_last_run', 'failed_jobs_save_path']:
        kwargs.pop(kwarg)
    multiplier = AmrMultiplier(**kwargs)
    if check_outputs or check_last_run:
        if check_outputs:
            multiplier.check_outputs_exist(save_path=save_path)
        if check_last_run:
            multiplier.check_last_run(save_path=save_path)
    else:
        multiplier.launch()
