import sys
import os
import pandas as pd
import argparse
from mcod_prep.run_phase_redistributionworker import main as worker_main
from mcod_prep.utils.mcod_cluster_tools import submit_mcod
from mcod_prep.utils.logging import ymd_timestamp
from mcod_prep.mcod_mapping import MCoDMapper
from cod_prep.downloaders import (
    get_current_location_hierarchy,
    get_cause_map,
    get_redistribution_locations,
    add_code_metadata,
    add_age_metadata,
)
from cod_prep.utils import (
    report_if_merge_fail,
    fill_missing_df,
    just_keep_trying,
    print_log_message,
    wait_for_job_ids,
)
from cod_prep.claude.configurator import Configurator
from mcod_prep.utils.causes import get_all_related_syndromes
from mcod_prep.utils.mcause_io import (
    get_phase_output, write_phase_output,
    makedirs_safely
)

CONF = Configurator("standard")
CACHE_DIR = CONF.get_directory("db_cache")
CODE_DIR = CONF.get_directory("mcod_code")
RD_INPUTS_DIR = "FILEPATH"
PACKAGE_DIR = "FILEPATH"
SG_DIR = "FILEPATH"
SG_PARENT = "FILEPATH"


def convert_deaths_to_died(df, value_col, id_cols):
    """
    Convert deaths column to died boolean for datasets with two value columns.
    """
    return (
        df.eval(f"alive = {value_col} - deaths")
        .drop(columns=value_col)
        .melt(id_vars=id_cols, var_name="died", value_name=value_col)
        .assign(died=lambda d: (d.died == "deaths") * 1)
    )


def has_garbage(df):
    """Determine whether or not there are any garbage codes."""
    return (df["cause_id"] == 743).any()


def format_age_groups(df):
    """Convert age groups to simple ages."""
    df = add_age_metadata(
        df, ["simple_age"], force_rerun=False, block_rerun=True, cache_dir=CACHE_DIR
    )
    df.rename(columns={"simple_age": "age"}, inplace=True)
    return df


def add_rd_locations(df, lhh):
    """Merge on location hierarchy specific to redistribution."""
    rd_lhh = get_redistribution_locations(lhh)
    df = pd.merge(df, rd_lhh, on="location_id", how="left")
    report_if_merge_fail(df, "global", "location_id")
    report_if_merge_fail(df, "dev_status", "location_id")
    return df


def add_split_group_id_column(
    df, nid, extract_type_id, data_type_id, int_cause, id_col="split_group"
):
    """Add group IDs to a dataset."""
    group_cols = [
        "country",
        "subnational_level1",
        "nid",
        "extract_type_id",
        "year_id",
        int_cause,
    ]
    if "admissions" in df.columns:
        group_cols += ["died"]
    if int_cause == "infectious_syndrome":
        group_cols += ["sepsis"]
    g = df.groupby(group_cols)[df.columns[-1]].min().reset_index()[group_cols]
    g[id_col] = g.index + 1
    return df.merge(g, on=group_cols)


def get_columns_added():
    return [
        "global",
        "dev_status",
        "super_region",
        "region",
        "country",
        "subnational_level1",
        "subnational_level2",
        "sex",
        "age",
        "split_group",
        "freq",
        "site_id",
    ]


def format_columns_for_rd(df, int_cause, value_col="deaths"):
    """Ensure necessary columns are appropriately named and present."""
    df.rename(columns={"value": "cause", value_col: "freq"}, inplace=True)
    df["sex"] = df["sex_id"].copy()
    df["site_id"] = 2
    add_cols = get_columns_added()
    assert (
        len(set(add_cols) - set(df.columns)) == 0
    ), f"{set(add_cols) - set(df.columns)} are missing"
    return df


def read_append_split_groups(
    sg_list, nid, extract_type_id, cause_map, int_cause, csid, cmvid
):
    """Read and append split groups after redistribution."""
    sg_dfs = []
    for sg in sg_list:
        filepath = "FILEPATH"
        sg = just_keep_trying(
            pd.read_csv,
            args=[filepath],
            kwargs={"dtype": {"cause": object}},
            max_tries=250,
            seconds_between_tries=6,
            verbose=True,
        )
        if int_cause in MCoDMapper.infectious_syndromes:
            sg = reassign_syndrome(sg, int_cause, csid, cmvid)
        sg = merge_acause_and_collapse(sg, cause_map)
        sg_dfs.append(sg)
    df = pd.concat(sg_dfs)
    return df


def reassign_syndrome(df, int_cause, code_system_id, code_map_version_id):
    """
    When the underlying/admission cause is an infectious syndrome,
    then this should be the infectious syndrome associated with that case/death.
    Redistribution resets the underlying cause, so re-assign
    syndrome after this.
    """
    mapper = MCoDMapper(
        int_cause,
        drop_p2=False,
        code_system_id=code_system_id,
        code_map_version_id=code_map_version_id,
    )
    int_cause_map = mapper.prep_int_cause_map()
    df = MCoDMapper.map_cause_codes(df, int_cause_map, int_cause, cols_to_map=["cause"])
    mapped_col = "cause_" + int_cause
    df[mapped_col] = df[mapped_col].fillna(MCoDMapper.remainder_syndromes[0])
    if int_cause == "infectious_syndrome":
        # If the new underlying cause is a real syndrome, this automatically
        # becomes the syndrome
        df.loc[~df[mapped_col].isin(MCoDMapper.remainder_syndromes), int_cause] = df[
            mapped_col
        ]
    else:
        # If the new underlying cause is a real syndrome, override the syndrome
        # indicator where appropriate
        if "_comm" not in int_cause and "_hosp" not in int_cause:
            target_syndromes = get_all_related_syndromes(int_cause)

            df.loc[df[mapped_col].isin(target_syndromes), int_cause] = 1

            df.loc[
                ~df[mapped_col].isin(target_syndromes)
                & ~df[mapped_col].isin(MCoDMapper.remainder_syndromes),
                int_cause,
            ] = 0
        if "_comm" in int_cause:
            target_syndromes = get_all_related_syndromes(int_cause.replace("_comm", ""))

            df.loc[df[mapped_col].isin(target_syndromes), int_cause] = 1

            df.loc[~df[mapped_col].isin(target_syndromes), int_cause] = 0

        if "_hosp" in int_cause:
            target_syndromes = get_all_related_syndromes(int_cause.replace("_hosp", ""))

            df.loc[df[mapped_col].isin(target_syndromes), int_cause] = 0

            df.loc[
                ~df[mapped_col].isin(MCoDMapper.remainder_syndromes)
                & ~df[mapped_col].isin(target_syndromes),
                int_cause,
            ] = 0

    df = df.drop([mapped_col], axis="columns")
    return df


def revert_variables(df, value_col="deaths"):
    """Change things back to standard columns."""
    add_cols = get_columns_added()
    id_cols = list(set(df.columns) - set(add_cols))
    df.rename(columns={"freq": value_col}, inplace=True)
    if value_col == "admissions":
        id_cols.remove("died")
        df = (
            df.assign(deaths=lambda d: d["died"] * d[value_col])
            .groupby(id_cols, as_index=False)[value_col, "deaths"]
            .sum()
        )
    else:
        df = df[[value_col] + id_cols]
    return df


def submit_split_group(
    nid, extract_type_id, split_groups, csid, int_cause, data_type_id
):
    """Submit jobs by split group."""
    jobname = f"{int_cause}_redistributionworker_{nid}_{extract_type_id}"
    worker = "FILEPATH/run_phase_redistributionworker.py"
    params = [nid, extract_type_id, int_cause, csid, data_type_id]
    array_path = "FILEPATH"
    pd.DataFrame({"split_group": split_groups}).to_csv(
        "FILEPATH", index=False
    )
    jid = submit_mcod(
        jobname,
        "python",
        worker,
        params=params,
        cores=1,
        memory="3G",
        runtime="00:30:00",
        verbose=True,
        logging=True,
        num_tasks=len(split_groups),
    )
    return jid


def write_split_group_input(df, nid, extract_type_id, sg, int_cause):
    """Write split group input."""
    indir = "FILEPATH"
    makedirs_safely(indir)
    df.to_csv("FILEPATH", index=False)


def delete_split_group_output(nid, extract_type_id, sg, int_cause):
    """Delete the existing intermediate split group files."""
    indir = "FILEPATH"
    for_rd_path = "FILEPATH"
    post_rd_path = "FILEPATH"
    for path in [for_rd_path, post_rd_path]:
        if os.path.exists(path):
            os.unlink(path)


def merge_acause_and_collapse(df, cause_map):
    """Add acause column and collapse before appending split groups."""
    cause_map = cause_map[["cause_id", "value"]].copy()
    cause_map = cause_map.rename(columns={"value": "cause"})
    df = df.merge(cause_map, how="left", on="cause")
    df = df.drop(["cause", "split_group"], axis=1)
    df = df.groupby([col for col in df.columns if col != "freq"], as_index=False).sum()
    return df


def run_phase(
    df,
    csvid,
    nid,
    extract_type_id,
    lsvid,
    cmvid,
    csid,
    remove_decimal,
    value_col,
    data_type_id,
    int_cause,
):
    """String together processes for redistribution."""
    read_file_cache_options = {
        "block_rerun": True,
        "cache_dir": CACHE_DIR,
        "force_rerun": False,
        "cache_results": False,
    }
    cause_map = get_cause_map(code_map_version_id=cmvid, **read_file_cache_options)
    lhh = get_current_location_hierarchy(
        location_set_version_id=lsvid, **read_file_cache_options
    )

    print_log_message("Formatting data for redistribution")
    orig_sum = int(df[value_col].sum())
    if value_col == "admissions":
        id_cols = list(set(df.columns) - set(["admissions", "deaths"]))
        df = convert_deaths_to_died(df, value_col, id_cols)
    if remove_decimal:
        print_log_message("Removing decimal from code map")
        cause_map["value"] = cause_map["value"].apply(lambda x: x.replace(".", ""))
    df = add_code_metadata(
        df, ["value", "code_system_id"], code_map=cause_map, **read_file_cache_options
    )
    df = format_age_groups(df)
    df = add_rd_locations(df, lhh)
    df = fill_missing_df(df, verify_all=True)
    df = add_split_group_id_column(df, nid, extract_type_id, data_type_id, int_cause)
    df = format_columns_for_rd(df, int_cause, value_col=value_col)

    split_groups = list(df.split_group.unique())
    print_log_message("Submitting/Running split groups")
    for sg in split_groups:
        delete_split_group_output(nid, extract_type_id, sg, int_cause)
        write_split_group_input(
            df.query(f"split_group == {sg}"), nid, extract_type_id, sg, int_cause
        )
    if len(split_groups) > 1:
        jid = submit_split_group(
            nid, extract_type_id, split_groups, csid, int_cause, data_type_id
        )
        print_log_message("Waiting for splits to complete...")
        wait_for_job_ids([jid], 30)
        print_log_message("Done waiting. Appending them together")
    else:
        worker_main(
            nid, extract_type_id, split_groups[0], csid, int_cause, data_type_id
        )
    df = read_append_split_groups(
        split_groups, nid, extract_type_id, cause_map, int_cause, csid, cmvid
    )
    print_log_message(f"Done appending files - {len(df)} rows assembled")
    df = revert_variables(df, value_col=value_col)

    post_sum = int(df[value_col].sum())
    before_after_text = f"""
        Before GC redistribution: {orig_sum}
        After GC redistribution: {post_sum}
    """
    diff = abs(orig_sum - post_sum)
    diff_threshold = max(0.02 * orig_sum, 5)
    assert diff < diff_threshold, f"{value_col} not close"
    print_log_message(before_after_text)
    return df


def main(
    nid,
    extract_type_id,
    csvid,
    lsvid,
    csid,
    cmvid,
    remove_decimal,
    data_type_id,
    int_cause,
):
    """Download data, run phase, and output result."""
    df = get_phase_output("format_map", nid, extract_type_id, sub_dirs=int_cause)
    if "admissions" in df.columns:
        value_col = "admissions"
    else:
        value_col = "deaths"
    if has_garbage(df):
        print_log_message("Running redistribution")
        df = run_phase(
            df,
            csvid,
            nid,
            extract_type_id,
            lsvid,
            cmvid,
            csid,
            remove_decimal,
            value_col,
            data_type_id,
            int_cause,
        )
    else:
        print_log_message("No redistribution to do.")
        group_cols = list(set(df.columns) - set(["code_id"]))
        df = df.groupby(group_cols, as_index=False)[value_col, "deaths"].sum()

    write_phase_output(
        df, "redistribution", nid, extract_type_id, ymd_timestamp(), sub_dirs=int_cause
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Parameters for running redistribution"
    )
    parser.add_argument("nid", type=int)
    parser.add_argument("extract_type_id", type=int)
    parser.add_argument("csvid", type=int)
    parser.add_argument("lsvid", type=int)
    parser.add_argument("csid", type=int)
    parser.add_argument("cmvid", type=int)
    parser.add_argument("remove_decimal", type=bool)
    parser.add_argument("data_type_id", type=int)
    parser.add_argument("int_cause", type=str)
    args = parser.parse_args()

    main(**vars(args))
