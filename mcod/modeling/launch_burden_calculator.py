from filecmp import cmp
from builtins import str
from builtins import object
import os
import datetime
from db_queries import get_ids
import pandas as pd
import numpy as np
import argparse
import itertools
from collections import defaultdict
from mcod_prep.utils.mcod_cluster_tools import submit_mcod
from mcod_prep.utils.causes import (get_int_cause_hierarchy, get_level_parent_child_cause_dict,
                                    get_child_causes, get_most_detailed_inj_causes)
from mcod_prep.utils.covariates import get_int_cause_model_covariates, get_cov, get_covariate_id
from mcod_prep.mcod_mapping import MCoDMapper
from cod_prep.claude.configurator import Configurator
from mcod_prep.utils.mcause_io import makedirs_safely
from cod_prep.downloaders import (
    get_current_cause_hierarchy, get_age_weights,
    get_pop, getcache_age_aggregate_to_detail_map, 
    get_cause_map
)
from mcod_prep.utils.mcause_io import McauseResult

class BurdenCalculator(object):
    conf = Configurator('standard')
    cache_options = {
        'force_rerun': True, 'block_rerun': False, 'cache_dir': "standard",
        'cache_results': True, 'verbose': True}
    full_time_series = list(range(1990, conf.get_id('year_end') + 1))
    block_rerun = {'force_rerun': False, 'block_rerun': True}

    def __init__(self, run_filters):
        self.run_filters = self.validate_run_filters(run_filters)
        self.int_causes = self.run_filters['int_causes']
        self.processes = self.run_filters['processes']
        self.description = self.run_filters['model_description']
        self.end_product = self.run_filters['end_product']
        self.years = self.run_filters['year_id']
        self.custom = self.run_filters['custom']
        self.run_lasso = self.run_filters['no_lasso']
        self.no_squeeze = self.run_filters['no_squeeze']
        self.subnationals = self.run_filters['subnationals']

    def validate_run_filters(self, run_filters):
        if "run_model" in run_filters['processes']:
            assert set(run_filters['processes']) == set(["run_model"]), \
                "Modeling should be run independently for now."
        if run_filters['year_id'] == ['all']:
            run_filters.update({'year_id': self.full_time_series})
        else:
            run_filters.update({'year_id': [int(x) for x in run_filters['year_id']]})
        assert all(isinstance(x, int) for x in run_filters['year_id'])
        ordered_processes = []
        no_more_processes = False
        all_processes = McauseResult.processes
        if 'remove_predictions_template' in all_processes:
            all_processes.remove('remove_predictions_template')
        if run_filters['end_product'] == 'rdp':
            all_processes.remove('calculate_props')
        for process in all_processes:
            if process in run_filters['processes']:
                assert not no_more_processes, "Mind the gap!"
                ordered_processes.append(process)
            else:
                if len(ordered_processes) > 0:
                    no_more_processes = True
        if ('run_model' in run_filters['processes']) and (run_filters['no_lasso']):
            print("Appending date to model description")
            if run_filters['model_description'] == '':
                description = '{:%Y_%m_%d}'.format(datetime.datetime.now())
            else:
                description = '{:%Y_%m_%d}'.format(
                    datetime.datetime.now()
                ) + '_' + run_filters['model_description']
            run_filters.update({'model_description': description})
        if 'sepsis_fatal_model' in run_filters['int_causes']:
            assert run_filters['model_description'] == 'GBD2017', \
                'using this keyword for modeling sepsis mortality is deprecated'
        if set(run_filters['int_causes']) <= set(McauseResult.infectious_syndromes + ['sepsis']):
            assert self.conf.config_type == 'amr', \
                "please change FILEPATH/setting.yaml to 'amr'"
        else:
            assert self.conf.config_type == 'mcod'
        return run_filters

    def cache_resources(self, cache_functions_to_run_with_args):
        for cache_function, kwargs in cache_functions_to_run_with_args:
            function_name = cache_function.__name__
            cache_exists = cache_function(only_check_cache_exists=True, verbose=True, **kwargs)
            if cache_exists:
                print(f"No need to recache method {function_name} with args: {kwargs}")
            else:
                print(f"Running {function_name} with args: {kwargs}")
                kwargs.update(self.cache_options)
                cache_function(**kwargs)

    def get_int_cause_age_bins(self, int_cause):
        if self.end_product == 'rdp':
            if int_cause in ['pulmonary_embolism', 'right_hf', 'left_hf', 'unsp_hf', 'aki',
                             'arterial_embolism', 'hypertension', 'atherosclerosis', 'amyloidosis',
                             'crf', 'alc_hepatic_failure']:
                return [39, 195, 211, 222, 229, 47, 21]
            elif int_cause in ['x59', 'y34']:
                return [39, 24, 224, 229, 47, 268, 294]
            elif int_cause in ['sepsis', 'explicit_sepsis', 'hepatic_failure']:
                return [28, 5, 23, 24, 41, 234]
            elif int_cause in [
                'arf', 'pneumonitis', 'unsp_cns', 'renal_failure', 'arrhythmia', 'peritonitis',
                    'fe_acid_base', 'cerebral_palsy', 'plegia', 'cachexia', 'empyema']:
                return [28, 5, 23, 24, 224, 229, 47, 30, 160]
            elif int_cause in ['gi_bleeding', 'pneumothorax', 'cardiac_arrest', 'osteomyelitis']:
                return [1, 189, 375, 217, 224, 229, 292, 160]
            else:
                AssertionError, f'no age bins specified for {int_cause}'
        elif int_cause == 'pulmonary_embolism_non_septic':
            return [39, 198, 227, 234]
        elif "respiratory" in int_cause:
            return [42, 179, 191, 26]
        elif "uti" in int_cause:
            return [175, 221]
        elif "bone" in int_cause: 
            return [172, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, 31, 32, 235]
        else:
            return self.conf.get_id('cod_ages')

    def get_detailed_causes(self, end_product, int_cause,
                            model_ages=None, custom=False,
                            **cache_kwargs):
        if end_product == 'rdp':
            if int_cause in ['x59', 'y34']:
                cause_list = get_most_detailed_inj_causes(
                    int_cause, cause_set_version_id=self.conf.get_id(
                        'reporting_cause_set_version'),
                    **cache_kwargs)
            else:
                if 'by_age' in self.description:
                    assert model_ages is not None
                    base_out_dir = McauseResult(
                        int_cause=int_cause, end_product=end_product,
                        process='run_model', description=self.description,
                        age_group_id=model_ages[0], conf=self.conf
                    ).results_path
                else:
                    base_out_dir = McauseResult(
                        int_cause=int_cause, end_product=end_product,
                        process='run_model', description=self.description,
                        conf=self.conf
                    ).results_path
                cause_file = f"{base_out_dir}/cause_list.csv"
                if 'by_age' in self.description:
                    assert all(
                        [cmp(f"{base_out_dir.parent}/{age}/cause_list.csv", cause_file)
                         for age in model_ages[1:]]
                    ), "cause list files differ across age groups"
                cause_list = list(pd.read_csv(cause_file).keep_causes.unique())
                bad_targets = [500, 543, 544, 729, 945, 843]
                cause_list = list(set(cause_list) - set(bad_targets))
        else:
            subset_string = 'yld_only != 1 & most_detailed == 1'
            if end_product in ['incidence', 'attributable_burden']:
                subset_string += ' & yll_only != 1'
            if custom:
                return [1, 2, 3]
            else:
                cause_list = list(get_current_cause_hierarchy(
                    cause_set_version_id=self.conf.get_id(
                        'computation_cause_set_version'),
                    cause_set_id=self.conf.get_id('computation_cause_set'), **cache_kwargs
                ).query(f'{subset_string}').cause_id.unique())
        return cause_list

    def get_causes_to_run(self, int_cause, model_ages):
        if self.run_filters['cause_id'] is not None:
            return self.run_filters['cause_id']
        else:
            return self.get_detailed_causes(
                self.end_product, int_cause,
                model_ages, self.custom, **self.block_rerun
            )

    def get_all_reporting_causes(self):
        if self.run_filters['cause_id'] is not None:
            return self.run_filters['cause_id']
        else:
            subset = 'yld_only != 1'
            if self.end_product in ['incidence', 'attributable_burden']:
                subset += ' & yll_only != 1'
            return list(get_current_cause_hierarchy(
                cause_set_version_id=self.conf.get_id('reporting_cause_set_version'),
                cause_set_id=self.conf.get_id('reporting_cause_set'), **self.block_rerun
            ).query(f'{subset}')['cause_id'].unique())

    def set_process_cause_lists(self, int_cause, model_ages):
        self.process_cause_lists = {}
        for process in set(McauseResult.processes) - {'run_model'}:
            if process in [
                'predict', 'calculate_counts', 'age_loc_aggregation',
                'remove_predictions_template'
            ]:
                cause_list = self.get_causes_to_run(int_cause, model_ages)
            elif process == 'cause_aggregation':
                if self.end_product == 'rdp':
                    cause_list = ['_all_int_cause']
                elif self.run_filters['cause_id'] is not None:
                    cause_list = self.run_filters['cause_id']
                else:
                    level_cause_dict = get_level_parent_child_cause_dict(
                        self.end_product
                    )
                    cause_list = []
                    for level in list(level_cause_dict.values()):
                        cause_list += list(level.keys())
            elif process == 'calculate_props':
                cause_list = self.get_all_reporting_causes()
            elif process == 'compile':
                if self.end_product == 'rdp':
                    cause_list = []
                else:
                    cause_list = self.get_all_reporting_causes()
            self.process_cause_lists[process] = cause_list

    def launch_compile(self, int_cause, holds=[]):
        if self.end_product == 'rdp':
            cause_list = self.process_cause_lists['age_loc_aggregation']
            worker = os.path.join("FILEPATH", 'compile_rdp.py')
            jobname = "compile_rd_props_" + int_cause
            params = [self.description, int_cause] + \
                ['--year_id '.join([''] + [str(x) + ' ' for x in self.years]).strip()] + \
                ['--target '.join([''] + [str(x) + ' ' for x in cause_list]).strip()]
            memory = self.calculate_memory(len(cause_list))
            McauseResult(
                int_cause=int_cause, end_product=self.end_product, process='compile',
                description=self.description
            ).clear_results()
            submit_mcod(jobname, 'python', worker, cores=8, verbose=True, logging=True,
                        memory=f'{memory}G', params=params, holds=holds, runtime="02:00:00")
        else:
            worker = os.path.join("FILEPATH", 'compile_burden.py')
            for cause_id in self.process_cause_lists['compile']:
                for year_id in self.years:
                    jobname = f"compile_summary_files_{int_cause}_{cause_id}_{year_id}"
                    params = [self.description, self.end_product, int_cause, cause_id, year_id]
                    if self.custom:
                        raise NotImplementedError
                        params += ['--custom']
                    submit_mcod(jobname, 'python', worker, cores=1, verbose=True, logging=True,
                                memory=f'52G', params=params, holds=holds)

    def calculate_memory(self, n):
        return int((np.log(n) * 30) + 15)

    def launch_cause_aggregation_with_holds(self, holds, int_cause, year,
                                            kwarg_dict, worker):
        if not self.custom:
            level_cause_dict = get_level_parent_child_cause_dict(self.end_product)
        else:
            raise NotImplementedError
        holds_dict = defaultdict(lambda: None)
        for level in sorted(level_cause_dict.keys(), reverse=True):
            for parent_cause, children in level_cause_dict[level].items():
                jobname = f'cause_agg_{int_cause}_{year}_{parent_cause}'
                params = [self.description, year, int_cause, self.end_product,
                          f'--parent_cause_id {parent_cause}']
                params += ['--child_cause '.join([''] + [str(x) + ' ' for x in children]).strip()]
                memory = self.calculate_memory(len(children))
                hold_ids = [holds_dict[x] for x in children if holds_dict[x] is not None] + holds
                if len(children) > 20 and len(children) <= 30:
                    kwarg_dict.update({'runtime': '01:00:00'})
                elif len(children) > 30:
                    kwarg_dict.update({'runtime': '02:00:00'})
                kwarg_dict.update({'params': params, 'memory': f'{memory}G', 'holds': hold_ids})
                jid = submit_mcod(jobname, 'python', worker, 1, **kwarg_dict)
                holds_dict[parent_cause] = jid
        return holds_dict.values()

    def launch_cause_aggregator(self, year, int_cause, holds):
        worker = os.path.join("FILEPATH", 'cause_aggregation.py')
        kwarg_dict = {'verbose': True, 'logging': True, 'holds': holds, 'runtime': '01:00:00'}
        if self.end_product == 'rdp':
            jobname = f"cause_agg_{int_cause}_{year}"
            params = [self.description, year, int_cause,
                      self.end_product, '--parent_cause_id _all_int_cause']
            params += ['--child_cause '.join([''] + [
                str(x) + ' ' for x in self.process_cause_lists['age_loc_aggregation']]).strip()
            ]
            memory = self.calculate_memory(len(
                self.process_cause_lists['age_loc_aggregation']
            ))
            jid = submit_mcod(jobname, 'python', worker, cores=2,
                              memory=f'{memory}G', params=params, **kwarg_dict)
            return jid
        else:
            if self.run_filters['cause_id'] is None:
                self.launch_cause_aggregation_with_holds(holds, int_cause, year,
                                                         kwarg_dict, worker)
            else:
                print("Submitting cause aggregation jobs without cause-specific holds")
                for cause_id in self.run_filters['cause_id']:
                    child_causes = get_child_causes(cause_id, self.end_product, **self.block_rerun)
                    params = [self.description, year, int_cause,
                              self.end_product, f'--parent_cause_id {cause_id}']
                    params += [
                        '--child_cause '.join([''] + [str(x) + ' ' for x in child_causes]).strip()]
                    jobname = f'cause_agg_{int_cause}_{year}_{cause_id}'
                    if len(child_causes) > 20 and len(child_causes) <= 30:
                        kwarg_dict.update({'runtime': '01:00:00'})
                    elif len(child_causes) > 30:
                        kwarg_dict.update({'runtime': '02:00:00'})
                    memory = f"{self.calculate_memory(len(child_causes))}G"
                    submit_mcod(jobname, 'python', worker, 1,
                                memory=memory, params=params, **kwarg_dict)
            return None

    def _create_year_cause_array_job_csv(self, cause_list, filepath):
        rows = itertools.product(*[self.years, cause_list])
        template = pd.DataFrame.from_records(rows, columns=['year_id', 'cause_id'])
        template.to_csv(f'{filepath}.csv', index=False)
        return len(template)

    def remove_predictions_template(self, cause_list, int_cause, model_ages):
        for year_id in self.years:
            for cause_id in cause_list:
                if 'by_age' in self.description:
                    for age_group_id in model_ages:
                        McauseResult(
                            int_cause=int_cause, end_product=self.end_product,
                            process='predictions_template', year_id=year_id,
                            cause_id=cause_id, age_group_id=age_group_id,
                            description=self.description, conf=self.conf
                        ).clear_results()
                else:
                    McauseResult(
                        int_cause=int_cause, end_product=self.end_product,
                        process='predictions_template', year_id=year_id,
                        cause_id=cause_id, description=self.description,
                        conf=self.conf
                    ).clear_results()

    def launch_save_predictions_template(self, int_cause, model_ages):
        worker = os.path.join("FILEPATH", "prep_predictions_template.py")
        jobname = f"save_predict_template_{int_cause}"
        params = [int_cause, self.end_product, self.description] + \
            ['--age_group_id '.join([''] + [str(x) + ' ' for x in model_ages]).strip()]
        if self.end_product == 'rdp':
            params += ['--subnationals']
        if self.end_product != 'rdp' and self.subnationals == True: 
            params += ['--subnationals']
        self.remove_predictions_template(
            self.process_cause_lists['predict'], int_cause, model_ages
        )
        num_tasks = self._create_year_cause_array_job_csv(
            self.process_cause_lists['predict'],
            f'{self.parent_model_dir}/save_predict_template_array_job_template')
        jid = submit_mcod(jobname, 'python', worker, cores=1, params=params, verbose=True,
                          memory='10G', logging=True, runtime="00:30:00", num_tasks=num_tasks)
        return jid

    def launch_predictions(self, int_cause, model_ages, holds=[]):
        worker = os.path.join("FILEPATH", "run_predictionsworker.R")
        jobname = f"predict_{int_cause}"
        params = [self.description, int_cause, self.parent_model_dir]
        if 'simple_glm' in self.description:
            for year_id in self.years:
                year_jobname = jobname + f'_{year_id}'
                year_params = params + [year_id]
                jid = submit_mcod(year_jobname, 'r', worker, cores=5, verbose=True,
                                  params=year_params, memory='80G', logging=True,
                                  runtime="05:00:00", holds=holds)
        else:
            num_tasks = self._create_year_cause_array_job_csv(
                self.process_cause_lists['predict'],
                f'{self.parent_model_dir}/predictions_array_job_template')
            jid = submit_mcod(jobname, 'r', worker, cores=2, params=params, verbose=True,
                              memory='20G', logging=True, runtime="00:50:00",
                              num_tasks=num_tasks, holds=holds)
        return jid

    def launch_counts_calculator(self, int_cause, holds=[]):
        worker = os.path.join("FILEPATH", 'calculate_counts.py')
        jobname = f'counts_{int_cause}'
        num_tasks = self._create_year_cause_array_job_csv(
            self.process_cause_lists['calculate_counts'],
            f'{self.parent_model_dir}/counts_array_job_template'
        )
        jid = submit_mcod(
            jobname, 'python', worker, cores=2, memory='120G',
            params=[self.description, int_cause, self.end_product, self.no_squeeze], verbose=True,
            logging=True, runtime="02:00:00", holds=holds, num_tasks=num_tasks
        )
        return jid

    def launch_age_loc_agg(self, int_cause, holds=[]):
        num_tasks = self._create_year_cause_array_job_csv(
            self.process_cause_lists['age_loc_aggregation'],
            f'{self.parent_model_dir}/age_loc_agg_array_job_template')
        worker = os.path.join("FILEPATH", 'age_loc_aggregation.py')
        jobname = f'agg_{int_cause}'
        jid = submit_mcod(
            jobname, 'python', worker, cores=1, memory='5G',
            params=[self.description, self.end_product, int_cause],
            verbose=True, logging=True, runtime="01:00:00", holds=holds, num_tasks=num_tasks
        )
        return jid

    def launch_back_calculate_props(self, int_cause, holds=[]):
        worker = os.path.join("FILEPATH", 'back_calculate_proportions.py')
        jobname = f'props_{int_cause}'
        num_tasks = self._create_year_cause_array_job_csv(
            self.process_cause_lists['calculate_props'],
            f'{self.parent_model_dir}/props_array_job_template')
        mem = 8
        if self.end_product == 'mortality':
            mem = 50
        jid = submit_mcod(
            jobname, 'python', worker, cores=5, memory=f'{mem}G',
            params=[self.description, int_cause, self.end_product], verbose=True,
            logging=True, runtime="06:00:00", holds=holds, num_tasks=num_tasks
        )
        return jid

    def launch_model(self, int_cause, model_ages):
        jobname = int_cause + "_run_model_" + self.description + '_' + self.end_product
        params = [self.description, int_cause]
        if self.end_product == 'rdp':
            worker = os.path.join("FILEPATH", 'run_model_rdp.py')
            params += ['--age_group_id '.join([''] + [str(x) + ' ' for x in model_ages]).strip()]
            if not self.run_lasso:
                params += ['--run_lasso']
        else:
            worker = os.path.join("FILEPATH", 'run_model.py')
            params += [self.end_product] + \
                ['--age_group_id '.join([''] + [str(x) + ' ' for x in model_ages]).strip()]
        submit_mcod(jobname, 'python', worker, 1, '35G', params=params,
                    verbose=True, logging=True, runtime="00:30:00")

    def launch(self):
        cache_functions_to_run_with_args = [
            (get_current_cause_hierarchy, {
                'cause_set_version_id': self.conf.get_id('computation_cause_set_version'),
                'cause_set_id': self.conf.get_id('computation_cause_set')}),
            (get_current_cause_hierarchy, {
                'cause_set_version_id': self.conf.get_id('reporting_cause_set_version'),
                'cause_set_id': self.conf.get_id('reporting_cause_set')}),
            (get_pop, {'pop_run_id': self.conf.get_id('pop_run')}),
            (get_age_weights, {}),
            (getcache_age_aggregate_to_detail_map, {}),
            (get_cause_map, {'code_system_id': 1}), 
            (get_cause_map, {'code_system_id': 6})
        ]
        self.cache_resources(cache_functions_to_run_with_args)

        if self.description == 'best_models': 
            assert self.processes != 'run_model'
            CONF = Configurator('standard')
            version_df = pd.read_csv("FILEPATH")
            self.int_causes = version_df['int_cause'].tolist()
            update_description = True
        else: 
            update_description = False

        for int_cause in self.int_causes:
            predict_jobs = []
            counts_jobs = []
            age_loc_agg_jobs = []
            cause_agg_jobs = []
            props_jobs = []
            model_ages = self.get_int_cause_age_bins(int_cause)
            covariates = get_int_cause_model_covariates(int_cause)
            cache_functions_with_args = [
                (get_int_cause_hierarchy, {
                    'int_cause': int_cause,
                    'cause_set_version_id': self.conf.get_id('cause_set_version')})
            ]
            if len(set(['run_model', 'predict']) - set(self.processes)) < 2:
                for cov in covariates:
                    cache_functions_with_args.append(
                        (get_cov, {'covariate_id': get_covariate_id(cov)})
                    )
            self.cache_resources(cache_functions_with_args)
            if 'run_model' in self.processes:
                self.launch_model(int_cause, model_ages)
            else:
                if update_description: 
                    self.description = ''.join(version_df.loc[version_df['int_cause'] == int_cause, 'model_description'].tolist())

                self.parent_model_dir = McauseResult(
                    int_cause=int_cause,
                    end_product=self.end_product,
                    process='run_model',
                    description=self.description,
                    age_group_id=model_ages[0],
                    conf=self.conf
                ).parent_model_dir
                assert self.parent_model_dir.exists()
                self.set_process_cause_lists(int_cause, model_ages)
                for process in self.processes:
                    for year_id in self.years:
                        for cause_id in self.process_cause_lists[process]:
                            McauseResult(
                                int_cause=int_cause, end_product=self.end_product,
                                process=process, year_id=year_id, cause_id=cause_id,
                                description=self.description, conf=self.conf
                            ).clear_results()
                if 'predict' in self.processes:
                    template_jid = self.launch_save_predictions_template(
                        int_cause, model_ages)
                    jid = self.launch_predictions(
                        int_cause, model_ages, holds=[template_jid])
                    predict_jobs.append(jid)
                if 'calculate_counts' in self.processes:
                    jid = self.launch_counts_calculator(int_cause, holds=predict_jobs)
                    counts_jobs.append(jid)
                if 'age_loc_aggregation' in self.processes:
                    jid = self.launch_age_loc_agg(int_cause, holds=counts_jobs)
                    age_loc_agg_jobs.append(jid)
                if 'remove_predictions_template' in self.processes:
                    self.remove_predictions_template(
                        self.process_cause_lists['remove_predictions_template'],
                        int_cause, model_ages
                    )
                for year_id in self.years:
                    if 'cause_aggregation' in self.processes:
                        jid = self.launch_cause_aggregator(
                            year_id, int_cause, holds=age_loc_agg_jobs)
                        cause_agg_jobs.append(jid)
                if 'calculate_props' in self.processes:
                    jid = self.launch_back_calculate_props(
                        int_cause, holds=cause_agg_jobs)
                    props_jobs.append(jid)
                if 'compile' in self.processes:
                    if self.end_product == 'rdp':
                        compile_holds = cause_agg_jobs
                    else:
                        compile_holds = props_jobs
                    self.launch_compile(int_cause, holds=compile_holds)

    def check_output_exists(self):
        missing_processes = {}
        for process in self.processes:
            missing_int_causes = {}

            if self.description == 'best_models': 
                assert self.processes != 'run_model'
                CONF = Configurator('standard')
                version_df = pd.read_csv("FILEPATH")
                self.int_causes = version_df['int_cause'].tolist()
                update_description = True
            else: 
                update_description = False

            for int_cause in self.int_causes:
                if update_description: 
                    self.description = ''.join(version_df.loc[version_df['int_cause'] == int_cause, 'model_description'].tolist())
                model_ages = self.get_int_cause_age_bins(int_cause)
                self.set_process_cause_lists(int_cause, model_ages)
                year_cause_dict = {}
                for year in self.years:
                    missing_causes = []
                    for cause in self.process_cause_lists[process]:
                        result = McauseResult(
                            int_cause=int_cause, end_product=self.end_product,
                            process=process, year_id=year, cause_id=cause,
                            description=self.description, conf=self.conf
                        )
                        if result.results_path.exists():
                            size = os.path.getsize(str(result.results_path))
                            if size <= 0:
                                missing_causes.append(cause)
                                year_cause_dict.update({year: missing_causes})
                                missing_int_causes.update({int_cause: year_cause_dict})
                                missing_processes.update({process: missing_int_causes})
                        else:
                            missing_causes.append(cause)
                            year_cause_dict.update({year: missing_causes})
                            missing_int_causes.update({int_cause: year_cause_dict})
                            missing_processes.update({process: missing_int_causes})
        if len(missing_processes) > 0:
            print("The following year/causes were not available: {}".format(missing_processes))
        else:
            print("All year/causes are available!")
        return missing_processes


if __name__ == '__main__':
    processes = McauseResult.processes
    product_list = McauseResult.product_list
    int_causes = McauseResult.valid_intermediate_causes + McauseResult.infectious_syndromes
    parser = argparse.ArgumentParser(description='Launch intermediate cause analyses')
    parser.add_argument(
        'model_description',
        help='version tracking for model runs. A few sneaky features:'
        '1) if you pass an empty string, date is appended 2) add "by_age" to run separate models by age group '
        '3) add "simple_glm" for a fixed effect on cause_id'
        '4) use "best_models" to run all int_causes and their best model descriptions', type=str)
    parser.add_argument(
        'processes', help='processes to be run, e.g. age_loc_aggregation',
        type=str, nargs='+', choices=processes)
    parser.add_argument(
        'end_product', help='how the end result will be used',
        choices=product_list, type=str)
    parser.add_argument(
        '--int_causes', help='intermediate cause(s) of interest'
        'Note: if description is "best_models", all int_causes will run regardless of those selected',
        nargs='+', choices=int_causes, type=str, required=True)
    parser.add_argument(
        '--cause_id', help='causes for which to run post-modeling processes',
        nargs='*', type=int)
    parser.add_argument(
        '--year_id', help='years for which to run post-modeling processes',
        nargs='*', default=[2019])
    parser.add_argument(
        '--custom', help='Use custom cause aggregation', action='store_true')
    parser.add_argument(
        '--no_lasso', help='Do not run LASSO; argument is ignored unless end_product is "rdp"',
        action='store_false')
    parser.add_argument(
        '--check', action='store_true', help='Do not launch, just check for outputs')
    parser.add_argument(
        '--no_squeeze', action='store_true', help='For calculate_counts only'
        'Stops from squeezing prediction results to the sepsis envelope')
    parser.add_argument(
        '--subnationals', action='store_true',
        help='Predict for subnationals, default false unless rdp')

    args = vars(parser.parse_args())
    check = args['check']
    args.pop('check')
    launcher = BurdenCalculator(args)
    print(f"You've submitted the following arguments: {args}")
    if check:
        launcher.check_output_exists()
    else:
        launcher.launch()