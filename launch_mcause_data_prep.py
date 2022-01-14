"""Launch formatting, mapping, and redistribution steps for multilple cause data."""
import pandas as pd
import argparse
from builtins import object
from mcod_prep.utils.mcod_cluster_tools import submit_mcod
from mcod_prep.utils.nids import get_datasets, add_nid_metadata
from mcod_prep.utils.mcause_io import delete_phase_output, check_output_exists, makedirs_safely
from mcod_prep.mcod_mapping import MCoDMapper
from cod_prep.downloaders import (
    get_map_version, get_cause_map, get_current_cause_hierarchy, get_ages,
    get_current_location_hierarchy, get_remove_decimal
)
from cod_prep.claude.configurator import Configurator
from cod_prep.utils import print_log_message, cod_timestamp


class MCauseLauncher(object):
    conf = Configurator()
    assert conf.config_type in ['amr']
    cache_options = {
        'force_rerun': True,
        'block_rerun': False,
        'cache_dir': "standard",
        'cache_results': True,
        'verbose': True
    }
    source_memory_dict = {
        'SOURCE': "5G"
    }

    def __init__(self, run_filters):
        self.run_filters = run_filters
        self.location_set_version_id = self.conf.get_id('location_set_version')
        self.cause_set_version_id = self.conf.get_id('cause_set_version')
        self.mcod_code = self.conf.get_directory('mcod_code')
        self.datasets = self.prep_run_filters()
        self.log_base_dir = "FILEPATH"

    def prep_run_filters(self):
        datasets_kwargs = {
            'force_rerun': True,
            'block_rerun': False,
            'cache_results': True,
            'is_active': True
        }
        datasets_kwargs.update({k: v for k, v in self.run_filters.items() if k not in
                                ['intermediate_causes', 'phases']})
        index_cols = ['nid', 'extract_type_id']
        datasets = get_datasets(**datasets_kwargs).drop_duplicates(index_cols).set_index(
            index_cols)[['year_id', 'code_system_id', 'source', 'data_type_id']]
        datasets = datasets.assign(
            code_map_version_id=datasets['code_system_id'].apply(
                lambda x: get_map_version(x, 'YLL', 'best')),
            remove_decimal=datasets['code_system_id'].apply(
                lambda x: get_remove_decimal(x))
        )
        return datasets

    def cache_resources(self, cache_functions_to_run_with_args):
        """Cache metadata files."""
        for cache_function, kwargs in cache_functions_to_run_with_args:
            function_name = cache_function.__name__
            cache_exists = cache_function(only_check_cache_exists=True, verbose=True, **kwargs)
            if cache_exists:
                print(f"No need to recache method {function_name} with args: {kwargs}")
            else:
                print(f"Running {function_name} with args: {kwargs}")
                kwargs.update(self.cache_options)
                cache_function(**kwargs)

    def launch_format_map(self, year, source, int_cause, code_system_id,
                          code_map_version_id, nid, extract_type_id, data_type_id):
        """Submit qsub for format_map phase."""
        delete_phase_output('format_map', nid, extract_type_id, sub_dirs=int_cause)
        worker = f"{self.mcod_code}/run_phase_format_map.py"
        params = [int(year), source, int_cause, int(code_system_id), int(code_map_version_id),
                  int(nid), int(extract_type_id), int(data_type_id)]
        jobname = f'mcause_format_map_{source}_{year}_{int_cause}'
        try:
            memory = self.source_memory_dict[source]
        except KeyError:
            print(f"{source} is not in source_memory_dict. Trying with 5G.")
            memory = '5G'

        if source in ['SOURCE']:
            runtime = '06:00:00'
        else:
            runtime = '00:45:00'
        if source in ['SOURCE']:
            jdrive = True
        else:
            jdrive = False
        jid = submit_mcod(
            jobname, 'python', worker, cores=1, memory=memory, params=params,
            verbose=True, logging=True, jdrive=jdrive, runtime=runtime,
            log_base_dir=self.log_base_dir
        )
        return jid

    def launch_redistribution(self, nid, extract_type_id, code_system_id, code_map_version_id,
                              remove_decimal, data_type_id, int_cause, holds=[]):
        """Submit qsub for redistribution."""
        delete_phase_output('redistribution', nid, extract_type_id, sub_dirs=int_cause)
        worker = "FILEPATH/run_phase_redistribution.py"
        jobname = f"{int_cause}_redistribution_{nid}_{extract_type_id}"
        params = [nid, extract_type_id, self.cause_set_version_id,
                  self.location_set_version_id, code_system_id, code_map_version_id,
                  remove_decimal, int(data_type_id), int_cause]
        submit_mcod(jobname, 'python', worker, cores=1, memory='10G', params=params,
                    holds=holds, verbose=True, logging=True, log_base_dir=self.log_base_dir)

    def launch(self):
        cache_functions_to_run_with_args = [
            (get_current_cause_hierarchy, {'cause_set_version_id': self.cause_set_version_id}),
            (get_ages, {}),
            (get_current_location_hierarchy,
                {'location_set_version_id': self.location_set_version_id})
        ]
        for code_map_version_id in list(self.datasets.code_map_version_id.unique()):
            cache_functions_to_run_with_args.append(
                (get_cause_map, {'code_map_version_id': code_map_version_id})
            )
        self.cache_resources(cache_functions_to_run_with_args)

        makedirs_safely(self.log_base_dir)

        format_map_jobs = []
        for row in self.datasets.itertuples():
            nid, extract_type_id = row.Index
            for int_cause in self.run_filters['intermediate_causes']:
                if 'format_map' in self.run_filters['phases']:
                    jid = self.launch_format_map(
                        row.year_id, row.source, int_cause, row.code_system_id,
                        row.code_map_version_id, nid, extract_type_id, row.data_type_id
                    )
                    format_map_jobs.append(jid)
                if 'redistribution' in self.run_filters['phases']:
                    self.launch_redistribution(
                        nid, extract_type_id, row.code_system_id, row.code_map_version_id,
                        row.remove_decimal, row.data_type_id, int_cause, holds=format_map_jobs
                    )

    def check_output_exists(self):
        failures = pd.DataFrame()
        for row in self.datasets.itertuples():
            nid, extract_type_id = row.Index
            for int_cause in self.run_filters['intermediate_causes']:
                for phase in self.run_filters['phases']:
                    if not check_output_exists(phase, nid, extract_type_id, sub_dirs=int_cause):
                        failures = failures.append({
                            'nid': nid,
                            'extract_type_id': extract_type_id,
                            'int_cause': int_cause,
                            'phase': phase
                        }, ignore_index=True)
        if len(failures) == 0:
            print_log_message("Outputs all exist!")
        else:
            failures = add_nid_metadata(failures, ['source'])
            print_log_message(f"The following nid/etid/phase/int causes failed: \n {failures}")
        return failures


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='NID/etids for which to run formatting, mapping, and redistribution')
    parser.add_argument('intermediate_causes', help='intermediate cause(s) of interest',
                        nargs='+', choices=MCoDMapper.possible_int_causes, type=str)
    parser.add_argument('phases', help='data processing phases', type=str,
                        nargs='+', choices=['format_map', 'redistribution'])
    parser.add_argument('--iso3', nargs='*')
    parser.add_argument('--code_system_id', nargs='*', type=int)
    parser.add_argument('--data_type_id', nargs='*', type=int)
    parser.add_argument('--year_id', type=int, nargs='*')
    parser.add_argument('--source', nargs='*')
    parser.add_argument('--nid', nargs='*', type=int)
    parser.add_argument('--check', action='store_true')
    args = vars(parser.parse_args())
    check = args['check']
    args.pop('check')
    launcher = MCauseLauncher(args)
    print(f"You've submitted the following arguments: {args}")
    if check:
        launcher.check_output_exists()
    else:
        launcher.launch()
