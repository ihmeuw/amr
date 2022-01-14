import pandas as pd
from cod_prep.utils import (
    report_duplicates, wait_for_job_ids,
    print_log_message
)
from mcod_prep.utils.mcod_cluster_tools import submit_mcod
from cod_prep.claude.configurator import Configurator
from modelling.pathogen.run_model import PathogenNetwork

CONF = Configurator()
HOLDS = []


def main():
    config_file = pd.read_excel("FILEPATH", sheet_name='run')
    report_duplicates(config_file, ['model_version', 'infectious_syndrome'])
    worker = "FILEPATH/run_model.py"
    print_log_message(f"Submitting jobs for {len(config_file)} models...")
    for model in config_file.to_dict('rows'):
        params = [
            model['model_version'],
            model['infectious_syndrome'],
        ]
        jobname = f"model_{model['model_version']}_{model['infectious_syndrome']}"
        submit_mcod(
            jobname, language='python', worker=worker,
            cores=1, memory="100G", params=params,
            runtime="24:00:00", logging=True,
            log_base_dir="FILEPATH",
            holds=HOLDS
        )
    print_log_message("Finished submitting!")


if __name__ == '__main__':
    main()
