"""
Read in split groups and pass them through redistribution.
"""
import sys
import os
import pandas as pd
import time
from mcod_prep.utils.nids import get_datasets
from cod_prep.claude.redistribution import GarbageRedistributor
from cod_prep.claude.configurator import Configurator
from cod_prep.downloaders.nids import get_value_from_nid
from cod_prep.utils import print_log_message

CONF = Configurator('standard')
PACKAGE_DIR = "FILEPATH"
SG_DIR = "FILEPATH"


def read_cause_map(code_system_id):
    """Read in cause map csv."""
    df = pd.read_csv("FILEPATH")
    return df


def read_split_group(nid, extract_type_id, sg, int_cause):
    """Read in split group dataframe."""
    indir = "FILEPATH"
    df = pd.read_csv("FILEPATH", dtype={'cause': 'object'})
    return df


def write_split_group(df, nid, extract_type_id, sg, int_cause):
    """Write completed split group."""
    indir = "FILEPATH"
    df.to_csv("FILEPATH", index=False)


def run_pipeline(df, nid, extract_type_id, cause_map, code_system_id, sg,
                 data_type_id, int_cause, write_diagnostics=True):
    proportion_ids_to_add = [
        int_cause, 'died', 'sepsis', 'code_system_id', 'pathogen_from_cause'
    ]

    redistributor = GarbageRedistributor(
        code_system_id, package_dir=PACKAGE_DIR, first_and_last_only=False,
        col_meta={c: {'col_type': 'demographic'} for c in proportion_ids_to_add}
    )
    redistributor.proportion_ids += proportion_ids_to_add
    df = redistributor.get_computed_dataframe(df, cause_map)

    if write_diagnostics:
        outdir = "FILEPATH"

        signature_metadata = redistributor.get_signature_metadata()
        signature_metadata.to_csv("FILEPATH", index=False)

        proportion_metadata = redistributor.get_proportion_metadata()
        proportion_metadata.to_csv("FILEPATH", index=False)

        magic_table = redistributor.get_diagnostic_dataframe()
        magic_table.to_csv("FILEPATH", index=False)

    return df


def main(nid, extract_type_id, split_group, csid, int_cause, data_type_id):
    start_time = time.time()
    df = read_split_group(nid, extract_type_id, split_group, int_cause)
    cause_map = read_cause_map(code_system_id)
    df = run_pipeline(df, nid, extract_type_id, cause_map,
                      code_system_id, split_group, data_type_id, int_cause)
    write_split_group(df, nid, extract_type_id, split_group, int_cause)
    run_time = time.time() - start_time
    print_log_message(f"Total run time: {run_time}")


if __name__ == "__main__":
    nid = int(sys.argv[1])
    extract_type_id = int(sys.argv[2])
    int_cause = str(sys.argv[3])
    code_system_id = int(sys.argv[4])
    data_type_id = int(sys.argv[5])
    task_id = os.environ.get('SGE_TASK_ID')
    if task_id:
        print(f'Running as array job, task_id: {task_id}')
        rd_process_dir = CONF.get_directory('rd_process_data')
        split_group = pd.read_csv(
            "FILEPATH"
        ).iloc[int(task_id) - 1]['split_group']
        print(f'Working on split group: {split_group}')
    main(nid, extract_type_id, split_group, code_system_id, int_cause, data_type_id)
