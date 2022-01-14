from builtins import str
from builtins import range
import os
import sys
import pandas as pd
import argparse
from cod_prep.utils import print_log_message, report_duplicates
from cod_prep.downloaders import get_current_cause_hierarchy
from cod_prep.claude.configurator import Configurator
from mcod_prep.age_loc_aggregation import create_age_standardized_rates
from mcod_prep.utils.mcause_io import McauseResult

CONF = Configurator('standard')
DEM_COLS = ['cause_id', 'location_id', 'sex_id', 'year_id', 'age_group_id']
DRAW_COLS = ["draw_" + str(x) for x in range(0, 1000)]
CACHE_KWARGS = {'force_rerun': False, 'block_rerun': True}


def main(description, year_id, int_cause, end_product, parent_cause_id, custom, child_causes):
    print_log_message("Reading in age/sex/location aggregated files")
    cause_dfs = []
    for cause_id in child_causes:
        print(f"Working on child cause {cause_id}")
        df = McauseResult(
            int_cause=int_cause, end_product=end_product, process='age_loc_aggregation',
            year_id=year_id, cause_id=cause_id, description=description, conf=CONF
        ).read_results().query('age_group_id != 27')
        cause_dfs.append(df)
    df = pd.concat(cause_dfs, ignore_index=True, sort=True)
    df['cause_id'] = parent_cause_id
    df = df.groupby(DEM_COLS, as_index=False)[DRAW_COLS].sum()
    if end_product in ['mortality', 'incidence']:
        print_log_message("Creating age standardized rates")
        age_std_df = create_age_standardized_rates(df, DRAW_COLS, DEM_COLS)
        df = pd.concat([age_std_df, df], ignore_index=True, sort=True)
        df[DEM_COLS] = df[DEM_COLS].astype(int)
    report_duplicates(df, DEM_COLS)
    print_log_message("Saving output")
    McauseResult(
        int_cause=int_cause, end_product=end_product, process='cause_aggregation',
        year_id=year_id, cause_id=parent_cause_id, description=description, conf=CONF
    ).write_results(df)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Aggregate child causes to parent')
    parser.add_argument('description', type=str)
    parser.add_argument('year_id', type=int)
    parser.add_argument('int_cause', type=str)
    parser.add_argument('end_product', type=str)
    parser.add_argument('--parent_cause_id', default=None)
    parser.add_argument('--custom', action='store_true')
    parser.add_argument('--child_causes', action='append', type=int)
    args = parser.parse_args()
    main(**vars(args))
