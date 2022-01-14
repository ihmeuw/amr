import os
import argparse
import glob
import pandas as pd
import warnings
from cod_prep.utils import print_log_message, clean_icd_codes
from cod_prep.downloaders import get_map_version
from cod_prep.claude.configurator import Configurator
from mcod_prep.mcod_mapping import MCoDMapper
import numpy as np
from multiprocessing import Pool
from functools import partial
from amr_prep.utils.amr_io import save_amr_data

this_dir = "FILEPATH"
repo_dir = "FILEPATH"

MAP_DIR = "FILEPATH"
L_DIR = "FILEPATH"
CONF = Configurator()
DEFAULT_CAUSE = "none"
DEFAULT_RESISTANCE = "unknown"

pathogen_map = pd.read_csv("FILEPATH".format(MAP_DIR))
antibiotic_map = pd.read_csv("FILEPATH".format(MAP_DIR))
specimen_map = pd.read_csv("FILEPATH".format(MAP_DIR))
specimen_map = specimen_map[['raw_specimen', 'specimen']].drop_duplicates()


def get_formatted_data(source_dir, arg):

    print_log_message('Getting unmapped data from:')
    print(source_dir)

    files = glob.glob(source_dir + "FILEPATH")
    if any('FILENAME' in file for file in files):
        if arg == 'new':
            if any('FILENAME' not in file for file in files):
                df = pd.read_csv(source_dir + "FILEPATH")
        else:
            df = pd.read_csv(source_dir + "FILEPATH")
    else:
        return None
    df = df.rename(columns={'specimen': 'raw_specimen'})
    if 'sample_id' in df:
        assert not ((df.sample_id.dtype == int) or (df.sample_id.dtype == float)),\
            "Need to do some kind of source + id appending here"

    return df


def map_pathogens(df, filepath):
    print_log_message('Mapping pathogens...')
    pathogen_map['raw_pathogen'] = pathogen_map['raw_pathogen'].str.replace('\r', '')
    df = df.merge(pathogen_map, how='left', on='raw_pathogen', validate='many_to_one')
    grant_bugs = [
        'escherichia_coli', 'klebsiella_pneumoniae', 'mycobacterium_tuberculosis',
        'neisseria_gonorrheae', 'non_typhoidal_salmonellae', 'salmonella_typhi_paratyphi',
        'shigella_spp', 'staphylococcus_aureus', 'streptococcus_pneumoniae'
    ]
    if ((df['pathogen'].loc[df['pathogen'].notna()].isin(grant_bugs).all()) &
            (df['pathogen'].isna().any())):
        warnings.warn(
            "Not all data mapped, but seems only have grant pathogens in " +
            str(filepath) + " Contact DATA_PROVIDER"
        )
    if ((df['pathogen'].loc[df['pathogen'].notna()].isin(grant_bugs).all()) &
            (df['pathogen'].notna().any())):
        warnings.warn("Only have grant pathogens in " + str(filepath) + " Contact DATA PROVIDER")
    return df


def map_antibiotics(df):
    print_log_message('Mapping antibiotics...')
    df = df.merge(antibiotic_map, how='left', on='raw_antibiotic', validate='many_to_one')

    return df


def map_specimens(df):
    print_log_message('Mapping specimen...')
    df = df.merge(specimen_map, how='left', on='raw_specimen', validate='many_to_one')

    print_log_message('Mapping specimen to syndrome...')
    specimen_to_syndrome = pd.read_csv("FILEPATH").set_index(
        "specimen")['infectious_syndrome'].to_dict()
    df['specimen_syndrome'] = df['specimen'].map(specimen_to_syndrome)

    return df


def map_cause_info(df, nid_metadata, filepath):
    cause_cols = [c for c in df if 'cause' in c]
    has_cause = len(cause_cols) > 0
    has_ucause = 'cause' in cause_cols
    if has_cause:
        code_system_id = nid_metadata.query(
            f"file_path == '{filepath}'").code_system_id.unique()[0]
        if not np.isnan(code_system_id):
            code_map_version_id = get_map_version(
                code_system_id, gbd_round_id=CONF.get_id("gbd_round"))
            if code_system_id in [1, 6]:
                for col in cause_cols:
                    df[col] = clean_icd_codes(df[col], remove_decimal=True)
            mapper = MCoDMapper(
                'infectious_syndrome', drop_p2=True,
                code_system_id=code_system_id, code_map_version_id=code_map_version_id
            )
        else:
            path_to_custom_map = "FILEPATH"
            mapper = MCoDMapper(
                'infectious_syndrome', drop_p2=True,
                path_to_ucause_map=path_to_custom_map,
                path_to_int_cause_map=path_to_custom_map
            )
        df = mapper.get_computed_dataframe(df, map_underlying_cause=has_ucause)
        df = df.drop(cause_cols, axis='columns')

    if has_cause or 'specimen_syndrome' in df:
        df['final_syndrome'] = np.NaN
        if has_cause:
            df['final_syndrome'] = df['infectious_syndrome']
        if 'specimen_syndrome' in df:
            df.loc[
                (df.final_syndrome.isnull() |
                 df.final_syndrome.isin(MCoDMapper.remainder_syndromes)),
                'final_syndrome'
            ] = df['specimen_syndrome']
        if 'pathogen' in df:
            pathogen_to_syndrome = {
                'salmonella_typhi': 'typhoid_fever',
                'salmonella_paratyphi': 'paratyphoid_fever',
                'salmonella_typhi_paratyphi': 'typhoid_paratyphoid_ints',
                'non_typhoidal_salmonellae': 'ints',
                'mycobacterium_tuberculosis': 'tb',
                'neisseria_meningitidis': 'cns_infectious',
                'neisseria_gonorrheae': 'chlamydia_and_gonorrheae',
            }
            for pathogen, syndrome in pathogen_to_syndrome.items():
                df.loc[
                    (df.pathogen == pathogen) &
                    (df.final_syndrome.isnull() |
                     df.final_syndrome.isin(MCoDMapper.remainder_syndromes)),
                    'final_syndrome'
                ] = syndrome
        df['infectious_syndrome'] = df['final_syndrome']
        assert df['infectious_syndrome'].notnull().all()
    if 'cause_id' in df:
        assert df['cause_id'].notnull().all()
    df = df.drop(
        [
            c for c in df if c in
            ['specimen_syndrome', 'code_id', 'cause_infectious_syndrome',
             'final_syndrome']
        ],
        axis='columns'
    )
    return df


def map_cause_info_in_parallel(df, num_workers, nid_metadata, filepath):
    print_log_message("Working on creating data chunks")
    samples = df['sample_id'].unique().tolist()
    groups = list(range(0, num_workers)) * (len(samples) // num_workers + 1)
    groups = groups[0:len(samples)]
    sample_to_group = dict(zip(samples, groups))
    df['group_id'] = df['sample_id'].map(sample_to_group)

    print_log_message("Checking for data with no cause info")
    cause_cols = [c for c in df if 'cause' in c]
    assert 'cause' in cause_cols
    assert len(cause_cols) > 0
    df['no_info'] = (df[cause_cols] == DEFAULT_CAUSE).all(axis=1)
    df['no_info'] = df.groupby('sample_id').no_info.transform(np.all)
    print(
        f"{df.no_info.sum()} out of {len(df)} rows have no cause info,"
        f" mapping separately..."
    )
    save_df = df.loc[df.no_info].copy().drop(cause_cols, axis='columns')
    save_df = map_cause_info(
        save_df, nid_metadata=nid_metadata,
        filepath=filepath
    )
    assert 'infectious_syndrome' in save_df
    save_df['cause_id'] = 919
    df = df.loc[~df.no_info]

    data_chunks = [chunk for group, chunk in df.groupby('group_id')]
    _map_cause_info = partial(
        map_cause_info, nid_metadata=nid_metadata, filepath=filepath
    )

    print_log_message(f"Working on running cause mapping")
    pool = Pool(num_workers)
    df_list = pool.map(_map_cause_info, data_chunks)
    pool.close()
    pool.join()
    print_log_message("Done, putting chunks back together")
    df = pd.concat(df_list, axis='index', ignore_index=True, sort=True)
    df = df.append(save_df)
    df = df.drop(['group_id', 'no_info'], axis='columns')
    return df


def map_all(arg, write=False):
    if write:
        print_log_message(
            "You passed write = True, this should only be done if your cleaning script"
            " has been approved in a PR review"
        )
    print_log_message("Searching for " + arg + " formatted files to map")
    nid_metadata = pd.read_csv("FILEPATH")

    if (arg != 'all') & (arg != 'new'):
        assert arg in nid_metadata['source'].unique()
        filepaths = nid_metadata.query(f"source == '{arg}'")['file_path'].unique().tolist()
    else:
        filepaths = nid_metadata['file_path'].unique().tolist()

    for filepath in filepaths:
        source = nid_metadata.query(f"file_path == '{filepath}'")['source'].unique()[0]
        df = get_formatted_data(filepath, arg)
        if df is None:
            print_log_message('No unmapped output found for: ' + filepath
                              + ' Has this been formatted?:')
            print_log_message(filepath)
            continue

        has_spec = False
        has_drug = False

        assert df[[c for c in df if c not in ['deaths']]].notnull().values.all(),\
            f"The following columns have nulls that need to be resolve before mapping: {df.isnull().any()}"
        assert 'nid' in df, "You must have a real NID"
        assert df.nid.dtype == 'int64', "You must have a real NID"
        assert (df.nid != -1).all(), "You must have a real NID"

        if 'raw_pathogen' in df:
            df = map_pathogens(df, filepath)
        if 'raw_antibiotic' in df:
            df = map_antibiotics(df)
            has_drug = True
        if 'raw_specimen' in df:
            df = map_specimens(df)
            has_spec = True

        print_log_message(
            'Proceding to read out full mapped file, or unmapped specimen, drug, bugs')

        mapped_cols = {
            'specimen': 'raw_specimen',
            'abx_class': 'raw_antibiotic',
            'pathogen': 'raw_pathogen'
        }

        if not has_spec:
            del mapped_cols['specimen']
        if not has_drug:
            del mapped_cols['abx_class']

        any_unmapped = False

        for key, value in mapped_cols.items():
            unmapped = df[key].isna()
            if unmapped.values.any():
                print_log_message('Unmapped {} in the data, outputting a list of unmapped {}. '
                                  'Please update the map whenever possible'.format(value, value))
                unmapped_values = df.loc[unmapped, [value]].drop_duplicates()
                any_unmapped = True
                if not os.path.exists(L_DIR + "FILEPATH"):
                    unmapped_values.to_csv(L_DIR + "FILEPATH")
                else:
                    existing_unmap = pd.read_csv(L_DIR + "FILEPATH")
                    unmapped_values = existing_unmap.append(unmapped_values)
                    unmapped_values = unmapped_values.drop_duplicates()
                    unmapped_values.to_csv(L_DIR + "FILEPATH")

        if source == 'lri_lit':
            df['cause'] = 'lower respiratory tract infections'

        if source == 'SOURCE_NAME':
            warnings.warn("This is SOURCE_NAME, make sure you're running with 50 threads")
            df = map_cause_info_in_parallel(df, 50, nid_metadata, filepath)
        else:
            df = map_cause_info(df, nid_metadata, filepath)

        if "hospital_type" in df.columns.tolist() or "hospital_name" in df.columns.tolist():
            assert "hospital_name" in df.columns.tolist()
            assert df[["hospital_type", "hospital_name"]].notnull().values.all(), "Not every hospital has been classified"

            df.drop("hospital_name", axis=1, inplace=True)
            df.loc[df['hospital_type'].isin(["unknown", "mixed"]), "hospital_type"] = "m/u"

            assert df.hospital_type.isin(["tertiary", "non-tertiary", "m/u"]).values.all(), "there are more hospital classifications in the data than allowed"

        if not any_unmapped:
            print_log_message('All specimen, raw_pathogen, and raw_antibiotic mapped!')
            df = df.drop(columns=list(mapped_cols.values()))
            if write:
                print_log_message('Saving standardized_mapped file')
                save_amr_data(df, phase='mapped', ldir=filepath)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="remap all or only new data")
    parser.add_argument(
        "type", help="'override' to remap everything, 'new' for any datasets not mapped yet, "
        "or an individual source name"
    )
    write_parser = parser.add_mutually_exclusive_group(required=True)
    write_parser.add_argument(
        '--write', dest='write', action='store_true',
        help='Write the output file from mapping '
        '(only use once your cleaning script has been approved)'
    )
    write_parser.add_argument(
        '--no-write', dest='write', action='store_false',
        help='Do not write the output file from mapping (used for testing purposes)'
    )
    arg = parser.parse_args()

    map_all(arg.type, write=arg.write)
