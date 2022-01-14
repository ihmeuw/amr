from cod_prep.claude.configurator import Configurator
from pathlib import Path
from cod_prep.claude.claude_io import makedirs_safely
from cod_prep.utils import print_log_message, cod_timestamp
from amr_prep.utils.misc import get_prepped_csbg_universe
import pandas as pd
import glob
import os
from cod_prep.utils import wrap

CONF = Configurator()


def read_results_wrapper(burden_type, cause_id, year_id, infectious_syndrome,
                         pathogen):
    return AmrResult(
        process='calculate_amr_burden',
        burden_type=burden_type,
        year_id=year_id,
        cause_id=cause_id,
        infectious_syndrome=infectious_syndrome,
        pathogen=pathogen
    ).read_results()


class AmrResult():
    burden_types = ["fatal", "nonfatal"]
    processes = {
        "split_sepsis_syndrome": {
            'valid_burden_types': ['fatal'],
            'valid_params': ['cause_id', 'infectious_syndrome']
        },
        "calculate_mi_ratios": {
            'valid_burden_types': ['nonfatal'],
            'valid_params': ['infectious_syndrome'],
        },
        "split_pathogen": {
            'valid_burden_types': ['fatal', 'nonfatal'],
            'valid_params': ['cause_id', 'infectious_syndrome', 'pathogen']
        },
        "calculate_ylds_per_case": {
            'valid_burden_types': ['nonfatal'],
            'valid_params': ['infectious_syndrome']
        },
        "calculate_amr_props": {
            'valid_burden_types': ['fatal', 'nonfatal'],
            'valid_params': ['pathogen']
        },
        "calculate_amr_burden": {
            'valid_burden_types': ['fatal', 'nonfatal'],
            'valid_params': ['cause_id', 'infectious_syndrome', 'pathogen']
        },
        "aggregate_cause_syndrome": {
            'valid_burden_types': ['fatal', 'nonfatal'],
            'valid_params': ['cause_id', 'infectious_syndrome', 'pathogen']
        },
        "summarize_burden": {
            'valid_burden_types': ['fatal', 'nonfatal'],
            'valid_params': ['cause_id', 'infectious_syndrome', 'pathogen']
        },
    }

    conf = Configurator()

    def __init__(self, process: str, burden_type: str, year_id: int,
                 cause_id: int = None, infectious_syndrome: str = None,
                 pathogen: str = None, abx_class: str = None,
                 save_draws: bool = False, version: str = None):
        self.process = process
        self.burden_type = burden_type
        self.year_id = year_id
        if (save_draws) & (process == 'summarize_burden'):
            self.draws_suffix = "_draws"
        else:
            self.draws_suffix = ""
        self.optional_params = {
            'cause_id': cause_id or "",
            'infectious_syndrome': infectious_syndrome or "",
            'pathogen': pathogen or "",
            'abx_class': abx_class or "",
        }
        assert self.conf.config_type == 'amr'
        self.version = version
        self.construct_results_path()

    def construct_results_path(self):
        if self.process not in AmrResult.processes:
            raise NotImplementedError
        if self.burden_type not in AmrResult.burden_types:
            raise NotImplementedError

        assert self.burden_type in AmrResult.processes[self.process]['valid_burden_types']
        for param, val in self.optional_params.items():
            if param in AmrResult.processes[self.process]['valid_params']:
                assert val != "", f"You must pass a valid {param} for process {self.process}"
                if param == 'cause_id':
                    self.optional_params[param] = int(val)
            else:
                assert val == "", f"You cannot pass {param} for process {self.process}"

        if self.process == 'aggregate_cause_syndrome':
            subfolder = "FILEPATH"
        else:
            subfolder = "FILEPATH"
        if self.version is not None:
            subfolder = "FILEPATH"
        self.results_path = "FILEPATH"

    def make_results_directory(self):
        makedirs_safely(str(self.results_path.parent))

    def clear_results(self):
        if self.version is not None:
            raise AssertionError('Cannot clear an archive file')
        if self.results_path.exists():
            self.results_path.unlink()

    def write_results(self, df):
        if self.version is not None:
            raise AssertionError('Cannot write to an archive file')
        self.make_results_directory()
        df.to_hdf(self.results_path, key='model', format='fixed', mode='w')

    def read_results(self):
        return pd.read_hdf(self.results_path, key='model')


def get_amr_results(process, burden_type, year_id, cause_id=None, infectious_syndrome=None,
                    pathogen=None, abx_class=None, location_id=None, age_group_id=None,
                    sex_id=None, measure_id=None, metric_id=None, hosp=None,
                    counterfactual=None, exec_function=None, exec_function_args=None,
                    filter_continuously=True, draws=False, version=None):
    exec_function_args = exec_function_args or {}
    if process != 'summarize_burden':
        raise NotImplementedError
    assert burden_type in AmrResult.burden_types

    csbg = get_prepped_csbg_universe(burden_type=burden_type, add_parents=True)
    csbg = csbg.loc[csbg.notnull().all(axis=1)]
    csbg_filters = {
        'cause_id': cause_id,
        'infectious_syndrome': infectious_syndrome,
        'pathogen': pathogen
    }
    for var, values in csbg_filters.items():
        if values is not None:
            csbg = csbg.loc[csbg[var].isin(wrap(values))]

    csbg = csbg[
        AmrResult.processes['summarize_burden']['valid_params']
    ].drop_duplicates()

    print_log_message(f"Reading {len(csbg) * len(wrap(year_id))} files")

    filters = {
        'abx_class': abx_class,
        'age_group_id': age_group_id,
        'sex_id': sex_id,
        'location_id': location_id,
        'measure_id': measure_id,
        'metric_id': metric_id,
        'hosp': hosp,
        'counterfactual': counterfactual
    }

    dfs = []
    i = 0
    for year in wrap(year_id):
        for index, row in csbg.iterrows():
            if i % 50 == 0:
                print(f"Reading file {i}")
            df = AmrResult(
                process, burden_type,
                year, save_draws=draws,
                version=version, **row.to_dict()
            ).read_results()
            if exec_function is not None:
                df = exec_function(df, **exec_function_args)
            if filter_continuously:
                for var, values in filters.items():
                    if values is not None:
                        df = df.loc[df[var].isin(wrap(values))]
            dfs.append(df)
            i += 1
    df = pd.concat(dfs, sort=False)

    if not filter_continuously:
        for var, values in filters.items():
            if values is not None:
                df = df.loc[df[var].isin(wrap(values))]
    return df


def save_amr_data(df: pd.DataFrame, phase: str, ldir: str):
    assert 'nid' in df, "You must have a real NID"
    assert df['nid'].notnull().values.all(), "You must have a real NID"
    assert df.nid.dtype == 'int64', "You must have a real NID"
    assert (df.nid != -1).all(), "You must have a real NID"

    if 'SOURCE_NAME' not in ldir:
        assert ('FILEPATH' in ldir) or (
            'FILEPATH' in ldir),\
            "You specified a directory not in the (FILEPATH) share folder"
    assert phase in ['unmapped', 'mapped']
    ldir = Path(ldir)
    assert ldir.exists()
    mapping_dir = ldir / "FILEPATH"
    archive_dir = mapping_dir / "FILEPATH"
    makedirs_safely(str(mapping_dir))
    makedirs_safely(str(archive_dir))
    file_name = "FILENAME.csv"
    archive_file_name = "FILENAME.csv"

    df.to_csv(mapping_dir / file_name, index=False)
    df.to_csv(archive_dir / archive_file_name, index=False)


def find_mapped_filepaths(target_sources_dict):
    target_source_path_dict = {}
    for key, value in target_sources_dict.items():
        for path in Path(value).rglob('*'):
            filepath = str(path)
            if "FILEPATH" in filepath:
                target_source_path_dict.update({key: path})

    return target_source_path_dict


def get_amr_data(sources=None, location_id=None, year_id=None, specimen=None,
                 pathogen=None, abx_class=None, infectious_syndrome=None, active=True):
    arguments = locals()

    nid_metadata = pd.read_csv("FILEPATH")

    if active:
        active_cols = [x for x in nid_metadata.columns if x not in [
            'source', 'file_path', 'code_system_id', 'DUA', 'IHME_Oxford']]
        nid_metadata = nid_metadata.loc[~(nid_metadata[active_cols] == 0).all(axis=1)]

    print_log_message('Getting the right filepaths')
    if sources is not None:
        if type(sources) is not list:
            sources = [sources]
        nid_metadata = nid_metadata.loc[nid_metadata['source'].isin(sources), ]
    wanted_source_dirs = nid_metadata.set_index('source')['file_path'].to_dict()

    dfs = []
    print_log_message('Reading in data')
    for key, value in wanted_source_dirs.items():
        mappath = "FILEPATH"
        if os.path.isfile(mappath):
            df = pd.read_csv(mappath)
            df['source'] = key
            dfs.append(df)
        else:
            print_log_message(f'{key} is not mapped yet, not gathering data')
    df = pd.concat(dfs)

    print_log_message('subsetting based on additional args')
    del arguments['sources']
    del arguments['active']
    for key, value in arguments.items():
        if value is not None:
            if type(value) is not list:
                value = [value]
            df = df.loc[df[key].isin(value), ]
    return df
