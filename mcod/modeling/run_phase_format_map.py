"""Run formatting and mapping steps for mcod data."""
from builtins import str
import sys
import pandas as pd
from importlib import import_module
from mcod_prep.mcod_mapping import MCoDMapper
from mcod_prep.utils.logging import ymd_timestamp, print_memory_timestamp
from mcod_prep.utils.mcause_io import write_phase_output, makedirs_safely
from cod_prep.claude.unmodeled_demographics_remap import AgeSexRemap
from cod_prep.utils import print_log_message, report_duplicates
from cod_prep.utils.schema import ColType
from cod_prep.downloaders.causes import get_all_related_causes, get_current_cause_hierarchy
from cod_prep.claude.configurator import Configurator

CONF = Configurator()
PROCESS_DIR = CONF.get_directory('mapping_process_data')
DIAG_DIR = PROCESS_DIR + '/{nid}/{extract_type_id}/{int_cause}'
ID_COLS = ['year_id', 'sex_id', 'age_group_id', 'location_id',
           'cause_id', 'code_id', 'nid', 'extract_type_id']


def get_drop_part2(int_cause, source):
    """
    Determine whether or not to drop part II of the death certificate.
    """
    if int_cause in MCoDMapper.infectious_syndromes:
        return True
    else:
        df = pd.read_csv(CONF.get_resource('drop_p2'), dtype={'drop_p2': bool})
        df = df[df['int_cause'] == int_cause]
        report_duplicates(df, 'int_cause')
        try:
            return df['drop_p2'].iloc[0]
        except IndexError:
            print(f"{int_cause} was not found in {CONF.get_resource('drop_p2')}")


def get_formatting_method(source, data_type_id, year, drop_p2):
    """Return the formatting method by source."""
    if data_type_id == 3:
        clean_source = 'clean_hospital_data'
        args = [source, year]
    else:
        clean_source = 'clean_' + source.lower()
        args = [year, drop_p2]
    try:
        formatting_method = getattr(
            import_module(f"mcod_prep.datasets.{clean_source}"), f"{clean_source}"
        )
    except AttributeError:
        print(f"No formatting method found! Check module & main function are named clean_{source}")
    return formatting_method, args


def drop_non_mcause(df, int_cause):
    chain_cols = [x for x in df.columns if ('multiple_cause_' in x) and ('pII' not in x)]
    df['num_chain_causes'] = 0
    for chain_col in chain_cols:
        df.loc[df[chain_col] != '0000', 'num_chain_causes'] += 1
    df['non_missing_chain'] = chain_cols[0]
    for chain_col in chain_cols:
        df.loc[
            (df['num_chain_causes'] == 1) & (df[chain_col] != '0000'), 'non_missing_chain'
        ] = df[chain_col].copy()
    drop_rows = (
        ((df['num_chain_causes'] == 1) & (df['cause'] == df['non_missing_chain']))
        | (df['num_chain_causes'] == 0)
    )
    df['drop_rows'] = 0
    df.loc[drop_rows, 'drop_rows'] = 1
    diag_df = df.groupby(
        ['nid', 'extract_type_id', 'location_id', 'year_id'], as_index=False
    ).agg({'num_chain_causes': ['mean', 'std', 'min', 'max'], 'drop_rows': sum})
    write_phase_output(diag_df, "format_map", nid, extract_type_id,
                       ymd_timestamp(), sub_dirs='diagnostic')
    df = df[~drop_rows]
    df = df.drop(['num_chain_causes', 'non_missing_chain', 'drop_rows'], axis=1)
    return df


def collect_cause_specific_diagnostics(df, acause_list):
    """Save diagnostic dataframe for certain causes."""
    cause_meta_df = get_current_cause_hierarchy(force_rerun=False, block_rerun=True)
    if type(acause_list) != list:
        acause_list = [acause_list]
    df_list = []
    for acause in acause_list:
        diag_df = df.loc[df['cause_id'].isin(get_all_related_causes(acause, cause_meta_df))]
        df_list.append(diag_df)
    diag_df = pd.concat(df_list, ignore_index=True)
    return diag_df


def get_id_value_cols(df, int_cause, data_type_id, inj_diag, drop_p2):
    """Collapse data so it is not longer at individual record level."""
    group_cols = ID_COLS + [int_cause]
    if int_cause == 'infectious_syndrome':
        group_cols += ['sepsis', 'pathogen_from_cause', 'cause_infectious_syndrome']
    if not drop_p2:
        group_cols += ['pII_' + int_cause]
    # make sure nothing was duplicated by accident
    group_cols = list(set(group_cols))
    # set value columns
    if data_type_id == 3:
        value_cols = ['admissions', 'deaths']
    else:
        value_cols = ['deaths']
    return group_cols, value_cols


def run_pipeline(year, source, int_cause, code_system_id, code_map_version_id,
                 nid, extract_type_id, data_type_id,
                 diagnostic_acauses='lri', inj_diag=False):
    """Clean, map, and prep data."""
    drop_p2 = get_drop_part2(int_cause, source)

    print_log_message("Prepping data")
    formatting_method, args = get_formatting_method(source, data_type_id, year, drop_p2)
    df = formatting_method(*args)

    print_log_message("Dropping rows without multiple cause")
    df = drop_non_mcause(df, int_cause)
    assert len(df) > 0, "No multiple cause data here!"

    # mapping from ICD code to code_id
    print_log_message("Mapping data")
    Mapper = MCoDMapper(
        int_cause, drop_p2, code_system_id=code_system_id,
        code_map_version_id=code_map_version_id)
    df = Mapper.get_computed_dataframe(df)

    if int_cause in MCoDMapper.infectious_syndromes:
        print_log_message("Mapping data again to get sepsis-related observations")
        Mapper = MCoDMapper(
            'sepsis', drop_p2, code_system_id=code_system_id,
            code_map_version_id=code_map_version_id)
        df = Mapper.get_computed_dataframe(df, map_underlying_cause=False)

        if int_cause != 'infectious_syndrome':
            print_log_message("Dropping non-sepsis related observations")
            df = df.loc[
                (df['sepsis'] != 'no_sepsis')
                & (df['infectious_syndrome'] != 'unspecified_and_none')
            ]
            assert len(df) > 0, "No sepsis-related rows!"

    # underlying cause specific explorations
    if diagnostic_acauses is not None:
        print_log_message(f"Saving {diagnostic_acauses} diagnostic")
        diag_df = collect_cause_specific_diagnostics(df, diagnostic_acauses)
        write_phase_output(diag_df, "diagnostic_causes", nid, extract_type_id,
                           ymd_timestamp(), sub_dirs=int_cause)

    group_cols, value_cols = get_id_value_cols(df, int_cause, data_type_id,
                                               inj_diag, drop_p2)
    print_memory_timestamp(df, f"Collapsing {value_cols} across {group_cols}")
    df = df.groupby(group_cols, as_index=False)[value_cols].sum()

    print_memory_timestamp(df, "Filtering cause-age-sex restrictions")
    cause_set_version_id = CONF.get_id('cause_set_version')
    column_metadata = {
        **{group_col: {'col_type': ColType.DEMOGRAPHIC} for group_col in group_cols},
        **{value_col: {'col_type': ColType.VALUE} for value_col in value_cols}
    }
    Remapper = AgeSexRemap(
        code_system_id, cause_set_version_id, collect_diagnostics=False, verbose=True,
        column_metadata=column_metadata
    )
    df = Remapper.get_computed_dataframe(df)

    return df


def main(year, source, int_cause, code_system_id, code_map_version_id,
         nid, extract_type_id, data_type_id):
    """Run pipeline."""
    df = run_pipeline(year, source, int_cause, code_system_id, code_map_version_id,
                      nid, extract_type_id, data_type_id)
    print_log_message(f"Writing nid {nid}, extract_type_id {extract_type_id}")
    write_phase_output(
        df, "format_map", nid, extract_type_id, ymd_timestamp(), sub_dirs=int_cause
    )


if __name__ == '__main__':
    year = int(sys.argv[1])
    source = str(sys.argv[2])
    int_cause = str(sys.argv[3])
    code_system_id = int(sys.argv[4])
    code_map_version_id = int(sys.argv[5])
    nid = int(sys.argv[6])
    extract_type_id = int(sys.argv[7])
    data_type_id = int(sys.argv[8])

    main(year, source, int_cause, code_system_id, code_map_version_id,
         nid, extract_type_id, data_type_id)
