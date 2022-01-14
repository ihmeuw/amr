import pandas as pd
from cod_prep.utils import (
    report_duplicates, wait_for_job_ids,
    print_log_message
)
from mcod_prep.utils.mcod_cluster_tools import submit_mcod
from cod_prep.claude.configurator import Configurator
import getpass
from cod_prep.downloaders import get_current_location_hierarchy
user = getpass.getuser()

CONF = Configurator()

IN_DIR = ADDRESS

corrs = pd.read_csv(FILEPATH)
corrs_salmonella = pd.read_csv(FILEPATH)
corrs = corrs.append(corrs_salmonella)

bugs = [ele for ele in reversed(corrs['pathogen'].value_counts().index.tolist())]

locs = get_current_location_hierarchy(
    location_set_version_id=CONF.get_id('location_set_version'),
)

countries = list(locs.loc[locs['level'] == 3, 'location_id'])

log_base_dir = FILEPATH
for bug in bugs:
    print_log_message(f"Submitting jobs for {bug} models by location...")
    for loc in countries:
        worker = f"{CONF.get_directory('amr_repo')}/modelling/"\
            "4_prevalence_resistance/simulate_resistance_profiles.py"
        params = [bug, loc]
        jobname = f"model_{bug}_{loc}"
        submit_mcod(
            jobname, language='python', worker=worker,
            cores=1, memory="10G", params=params,
            runtime="10:00:00", logging=True,
            log_base_dir=log_base_dir
        )
    print_log_message("Finished submitting!")
