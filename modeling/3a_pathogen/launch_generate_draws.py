import sys
import pandas as pd
from pathlib import Path
from cod_prep.utils import (
    print_log_message, wait_for_job_ids
)
from mcod_prep.utils.mcod_cluster_tools import submit_mcod
from cod_prep.claude.configurator import Configurator

CONF = Configurator()

SYNDROMES = [
    'blood_stream_infectious',
    'cns_infectious',
    'respiratory_infectious',
    'peritoneal_and_intra_abdomen_infectious',
    'skin_infectious',
    'uti_plus',
    'bone_joint_infection'
]


def main(diag_dir):
    best_mods = pd.read_csv(CONF.get_resource("best_pathogen_model_versions"))
    best_mods = best_mods.loc[best_mods.infectious_syndrome.isin(SYNDROMES)]
    best_mods = best_mods[['infectious_syndrome', 'model_version']].to_records(index=False)
    worker = "FILEPATH/generate_draws.py"
    jobs = []
    for mod in best_mods:
        params = [mod[0], mod[1], diag_dir]
        jobname = f"get_draws_{mod[0]}_{mod[1]}"
        jid = submit_mcod(
            jobname, language='python', worker=worker,
            cores=1, memory="100G", params=params,
            runtime="01:00:00", logging=True,
            log_base_dir=diag_dir
        )
        jobs.append(jid)
    print_log_message("Finished submitting!")
    print_log_message("Waiting...")
    wait_for_job_ids(jobs)
    print_log_message("Jobs complete")

    parent_dir = Path(diag_dir)

    df = pd.DataFrame()
    for file in parent_dir.iterdir():
        if '.csv' in str(file) and 'all_syndromes' not in str(file):
            df = df.append(pd.read_csv(file))
    df.to_csv("FILEPATH")


if __name__ == '__main__':
    diag_dir = str(sys.argv[1])
    main(diag_dir)
