from builtins import str
from builtins import range
import os
import sys
import pandas as pd
import numpy as np
from cod_prep.claude.configurator import Configurator
from cod_prep.utils import print_log_message, report_if_merge_fail
from mcod_prep.calculate_counts import get_deaths_draws
from mcod_prep.age_loc_aggregation import create_age_standardized_rates
from cod_prep.downloaders import get_cod_ages
from mcod_prep.utils.mcause_io import McauseResult

DEM_COLS = ['cause_id', 'location_id', 'sex_id', 'year_id', 'age_group_id']
DRAW_COLS = ["draw_" + str(x) for x in range(0, 1000)]
CONF = Configurator('standard')
BLOCK_RERUN = {'force_rerun': False, 'block_rerun': True}


def calculate_prop(df, draw):
    over_one = df[f"{draw}_x"] > df[f"{draw}_y"]
    outside_tol = ~np.isclose(df[f"{draw}_x"], df[f"{draw}_y"], rtol=0.01, atol=10)
    bad_rows = over_one & outside_tol
    assert not bad_rows.any(), f"Found the following bad rows with proportion "\
        f"greater than 1: \n {df.loc[bad_rows]}"
    return np.clip((df[f"{draw}_x"] / df[f"{draw}_y"]).fillna(0), None, 1)


def convert_count_to_prop(counts_df, deaths_df):
    print_log_message("Calculating proportions")
    df = counts_df.merge(deaths_df, how='left', on=DEM_COLS)
    report_if_merge_fail(df, DRAW_COLS[0] + '_y', DEM_COLS)
    for draw in DRAW_COLS:
        df[draw] = calculate_prop(df, draw)
    df = df[DEM_COLS + DRAW_COLS]
    print_log_message("Running checks")
    assert df.notnull().values.all(), \
        f"Error calculating proportions \n{df.columns[df.isnull().any()]}"
    assert (df[DRAW_COLS] <= 1).all().all(), \
        f"Proportion should not exceed 1"
    return df


def main(description, cause_id, year_id, end_product, int_cause):
    counts_df = McauseResult(
        int_cause=int_cause, end_product=end_product, process='cause_aggregation',
        year_id=year_id, cause_id=cause_id, description=description, conf=CONF
    ).read_results()
    counts_df[DEM_COLS] = counts_df[DEM_COLS].astype(int)
    counts_df = counts_df.query('age_group_id != 27')
    ages = list(counts_df.age_group_id.unique())
    locs = list(counts_df.location_id.unique())
    sexes = list(counts_df.sex_id.unique())
    deaths_df = get_deaths_draws(
        cause_id, year_id, ages, locs, sexes, int_cause, end_product,
        get_all_ages=False
    )
    df = convert_count_to_prop(counts_df, deaths_df)
    print_log_message("Creating age standardized rates")
    age_std_df = create_age_standardized_rates(df, DRAW_COLS, DEM_COLS)
    df = pd.concat([age_std_df, df], ignore_index=True, sort=True)
    print_log_message("Saving output")
    McauseResult(
        int_cause=int_cause, end_product=end_product, process='calculate_props',
        year_id=year_id, cause_id=cause_id, description=description, conf=CONF
    ).write_results(df)


if __name__ == '__main__':
    description = str(sys.argv[1])
    int_cause = str(sys.argv[2])
    end_product = str(sys.argv[3])
    assert end_product != 'rdp', 'Not used for redistribution props'
    task_id = os.environ.get('SGE_TASK_ID')
    if task_id:
        print(f'Running as array job, task_id: {task_id}')
        base_out_dir = McauseResult(
            int_cause=int_cause, end_product=end_product,
            process='run_model', description=description, conf=CONF
        ).parent_model_dir
        task_row = pd.read_csv("FILEPATH").iloc[int(task_id) - 1]
        cause_id = int(task_row['cause_id'])
        year_id = int(task_row['year_id'])
    else:
        cause_id = int(sys.argv[4])
        year_id = int(sys.argv[5])
    print_log_message(f'Running year {year_id}, cause {cause_id}')
    main(description, cause_id, year_id, end_product, int_cause)
