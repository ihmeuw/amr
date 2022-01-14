import pandas as pd
from cod_prep.claude.configurator import Configurator
from cod_prep.utils.cod_db_tools import (
    get_function_results,
    add_tbl_metadata
)
from cod_prep.utils.misc import report_if_merge_fail
from cod_prep.downloaders.locations import add_location_metadata
from db_tools import ezfuncs, query_tools
CONF = Configurator()


def get_datasets(nid=None, extract_type_id=None, source=None,
                 nid_extract_records=None, location_id=None, year_id=None,
                 data_type_id=None, code_system_id=None, iso3=None, location_set_id=None,
                 location_set_version_id=None, region_id=None, parent_nid=None, is_active=None,
                 verbose=False, force_rerun=False, block_rerun=True, cache_results=False):

    if is_active is not None:
        assert isinstance(is_active, bool), \
            "Pass either True or False to is_active"
        is_active = 1 * is_active

    dataset_filters = {
        'nid': nid,
        'extract_type_id': extract_type_id,
        'nid_extract_records': nid_extract_records,
        'source': source,
        'location_id': location_id,
        'year_id': year_id,
        'data_type_id': data_type_id,
        'code_system_id': code_system_id,
        'iso3': iso3,
        'region_id': region_id,
        'parent_nid': parent_nid,
        'is_active': is_active
    }
    dataset_filters = {
        k: dataset_filters[k] for k in dataset_filters
        if dataset_filters[k] is not None and k is not 'nid_extract_records'
    }
    use_nid_extract_records = False
    if nid_extract_records is not None:
        use_nid_extract_records = True

    cache_options = {
        'verbose': verbose,
        'force_rerun': force_rerun,
        'block_rerun': block_rerun,
        'cache_results': cache_results,
        'cache_dir': 'standard'
    }

    datasets = get_nidlocyear_map(**cache_options)
    add_cols = ['source', 'data_type_id', 'code_system_id', 'parent_nid', 'is_active']
    datasets = add_nid_metadata(datasets, add_cols, **cache_options)
    for need_col in ['source', 'data_type_id', 'code_system_id', 'is_active']:
        report_if_merge_fail(datasets, need_col, ['nid', 'extract_type_id'])

    if location_set_version_id is None:
        location_set_version_id = CONF.get_id('location_set_version')
    datasets = add_location_metadata(
        datasets, ['ihme_loc_id', 'region_id'],
        location_set_version_id=location_set_version_id,
        **cache_options
    )

    if location_set_id is not None:
        datasets = datasets.loc[datasets['ihme_loc_id'].notnull()]
    datasets['iso3'] = datasets['ihme_loc_id'].str.slice(0, 3)

    if use_nid_extract_records:
        nid_extract_df = pd.DataFrame.from_records(
            nid_extract_records, columns=['nid', 'extract_type_id']
        )
        datasets = datasets.merge(nid_extract_df)

    for var in list(dataset_filters.keys()):
        vals = dataset_filters[var]
        if not isinstance(vals, list):
            vals = [vals]
        datasets = datasets.loc[(datasets[var].isin(vals))]

    if len(datasets) == 0:
        raise AssertionError(
            "Given dataset filters produced no "
            "datasets: \n{}".format(dataset_filters)
        )

    return datasets


def get_nidlocyear_map(**cache_kwargs):
    nid_query = """SELECT nid, extract_type_id, location_id, year_id, representative_id
                    FROM cod.mcause_nid_location_year"""
    function = ezfuncs.query
    args = [nid_query]
    kwargs = {'conn_def': CONF.get_database_setup('db')}
    cache_name = "nid_locyears"

    df = get_function_results(
        function,
        args,
        kwargs,
        cache_name,
        **cache_kwargs
    )

    return df


def get_nid_metadata(**cache_kwargs):
    nid_query = """
        SELECT * FROM cod.mcause_nid_metadata
    """
    function = ezfuncs.query
    args = [nid_query]
    kwargs = {'conn_def': CONF.get_database_setup('db')}
    cache_name = "nid_metadata"

    df = get_function_results(
        function,
        args,
        kwargs,
        cache_name,
        **cache_kwargs
    )
    return df


def add_nid_metadata(df, add_cols, merge_col=['nid', 'extract_type_id'],
                     nid_meta_df=None, **cache_kwargs):

    if nid_meta_df is not None:
        tbl_or_function = nid_meta_df
    else:
        tbl_or_function = get_nid_metadata
    return add_tbl_metadata(
        tbl_or_function, df, add_cols, merge_col, **cache_kwargs
    )


def get_value_from_nid(nid, extract_type_id, value, force_rerun=False,
                       block_rerun=True, location_set_version_id=None, nid_meta_df=None):

    nid = int(nid)
    extract_type_id = int(extract_type_id)

    if nid_meta_df is None:
        nid_meta_df = get_datasets(
            nid=nid, extract_type_id=extract_type_id,
            force_rerun=force_rerun, block_rerun=block_rerun,
            location_set_version_id=location_set_version_id
        )
    else:
        nid_meta_df = nid_meta_df.query(f'nid == {nid} & extract_type_id == {extract_type_id}')
    assert value in nid_meta_df, f"{value} not available from get_datasets"
    nid_df = nid_meta_df[['nid', 'extract_type_id', value]].drop_duplicates().set_index('nid')
    assert len(nid_df) == 1, f"Expected one row for nid {nid}, found {len(nid_df)}"
    value = str(value)
    return nid_df.loc[nid, value]


def update_nid_metadata_status(nid, extract_type_id, is_active, conn_def='prodcod'):
    assert type(is_active) == bool
    is_active = int(is_active)

    myconn = ezfuncs.get_connection(conn_def)

    query_tools.exec_query("""
        UPDATE cod.mcause_nid_metadata
        SET is_active = {}
        WHERE nid = {}
        AND extract_type_id = {}
        """.format(is_active, nid, extract_type_id), cxn=myconn, close=True)
