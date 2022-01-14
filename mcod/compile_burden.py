import os
import sys
import argparse
import pandas as pd
from get_draws.api import get_draws
from multiprocessing import Pool
from functools import partial
from cod_prep.downloaders import pretty_print, get_current_cause_hierarchy, get_cod_ages
from cod_prep.claude.configurator import Configurator
from cod_prep.utils import print_log_message, report_if_merge_fail
from mcod_prep.utils.mcause_io import McauseResult

CONF = Configurator('standard')
DEM_COLS = ['cause_id', 'location_id', 'sex_id', 'year_id', 'age_group_id']
DRAW_COLS = ["draw_" + str(x) for x in range(0, 1000)]


def convert_int_cols(df):
    df[DEM_COLS] = df[DEM_COLS].apply(pd.to_numeric, downcast='integer')
    return df


def summarize_draws(df, draw_cols=DRAW_COLS, prefix=''):
    df[prefix + 'mean'] = df[draw_cols].mean(axis=1)
    df[prefix + 'upper'] = df[draw_cols].quantile(.975, axis=1)
    df[prefix + 'lower'] = df[draw_cols].quantile(.025, axis=1)
    df = df.drop(draw_cols, axis=1)
    return df


def compile_metrics(cause, year, description, end_product, int_cause):
    print_log_message(f"Reading in cause {cause} count files")
    df = convert_int_cols(
        McauseResult(
            int_cause=int_cause, end_product=end_product,
            process='cause_aggregation', description=description,
            year_id=year, cause_id=cause, conf=CONF
        ).read_results()
    )

    if end_product == 'mortality':
        print_log_message(f"Calculating cause {cause} specific chain fraction")
        denom_df = convert_int_cols(
            McauseResult(
                int_cause=int_cause, end_product=end_product,
                process='cause_aggregation', description=description,
                year_id=year, cause_id=294, conf=CONF
            ).read_results()
        ).drop('cause_id', axis=1)
        cscf_df = df.merge(denom_df, on=[x for x in DEM_COLS if x != 'cause_id'], how='left')
        cscf_df[DRAW_COLS] = cscf_df.filter(regex="draw_.*_x").rename(
            columns=lambda c: c[:-2]) / cscf_df.filter(regex="draw_.*_y").rename(
                columns=lambda c: c[:-2])
        cscf_df = summarize_draws(cscf_df, prefix='cscf_')[
            DEM_COLS + ['cscf_mean', 'cscf_upper', 'cscf_lower']
        ]

    print_log_message(f"Reading in cause {cause} proportion files")
    cf_df = McauseResult(
        int_cause=int_cause, end_product=end_product,
        process='calculate_props', description=description,
        year_id=year, cause_id=cause, conf=CONF
    ).read_results()
    cf_df = convert_int_cols(cf_df)
    prefix = 'cf_'
    if end_product == 'incidence':
        prefix = 'cfr_'
    cf_df = summarize_draws(cf_df, prefix=prefix)[DEM_COLS + ['cf_mean', 'cf_upper', 'cf_lower']]

    print_log_message(f"Merging {cause} metrics together")
    df = summarize_draws(df)[DEM_COLS + ['mean', 'upper', 'lower']]
    merge_df = df.merge(cf_df, on=DEM_COLS, how='outer', validate='one_to_one')\
        .merge(cscf_df, on=DEM_COLS, how='outer', validate='one_to_one')

    report_if_merge_fail(merge_df, 'mean', DEM_COLS)
    report_if_merge_fail(merge_df, 'cf_mean', DEM_COLS)
    report_if_merge_fail(merge_df, 'cscf_mean', DEM_COLS)
    assert merge_df.notnull().values.all()

    print_log_message(f"Job done: {cause}")
    return merge_df


def main(description, end_product, int_cause, custom, year_id, cause_id):
    df = compile_metrics(cause_id, year_id, description, end_product, int_cause)
    McauseResult(
        int_cause=int_cause, end_product=end_product,
        process='compile', description=description,
        year_id=year_id, cause_id=cause_id, conf=CONF
    ).write_results(df)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compile summary files')
    parser.add_argument('description', type=str)
    parser.add_argument('end_product', type=str)
    parser.add_argument('int_cause', type=str)
    parser.add_argument('cause_id', type=int)
    parser.add_argument('year_id', type=int)
    parser.add_argument('--custom', action='store_true')
    args = parser.parse_args()
    main(**vars(args))
