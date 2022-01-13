from cod_prep.claude.configurator import Configurator
from cod_prep.claude.formatting import (induce_extraction_type, pull_extract_type_id,
                                        pull_or_create_extract_type)
from cod_prep.utils import report_if_merge_fail, cod_timestamp, clean_icd_codes
from mcod_prep.utils.nids import get_nid_metadata, get_nidlocyear_map
from mcod_prep.utils.mcause_io import write_to_db_nid_table


def validate_cause_column(df, column):
    assert df[column].notnull().all()

    assert (df[column] != "").all()

    assert df[column].str.len().isin(range(3, 7)).all()

    if int(df['code_system_id'].unique()) == 1:
        assert df[column].str.match("^([A-Z][0-9]{2,}|0000)$").all(), \
            'not all columns are icd10-like'
    else:
        assert df[column].str.match("^(([0-9]|E[89])[0-9]{2,}|0000)$").all(), \
            'not all columns are icd9-like'


class MultipleCauseFormatter():

    def __init__(self, df, source, data_type_id, drop_p2, value_cols=['deaths'],
                 machine_learning=False):
        self.df = df
        self.source = source
        self.value_cols = value_cols
        self.data_type_id = int(data_type_id)
        self.icd_regex = "^[A-Z]?[0-9]+$"
        self.dem_cols = ['year_id', 'sex_id', 'age_group_id', 'nid', 'extract_type_id',
                         'location_id', 'code_system_id']
        if self.source in ["SOURCE"]:
            self.drop_p2 = True
        else:
            self.drop_p2 = drop_p2
        self.split_nids = [107077, 394317]
        self.conf = Configurator()
        self.machine_learning = machine_learning

    def check_formatting(self, df):
        if not isinstance(self.value_cols, list):
            self.value_cols = [self.value_cols]
        final_cols = self.dem_cols + self.value_cols
        for col in final_cols:
            df[col] = df[col].astype(int)
        cause_cols = ['cause'] + [x for x in df.columns if 'multiple_cause_' in x]
        if self.machine_learning:
            cause_cols = cause_cols + ['cause_code_original']
        assert len(cause_cols) > 1, "MCOD should have more than 1 cause"
        assert 'cause' in cause_cols, "need underlying cause of death"
        for col in cause_cols:
            if col != 'cause':
                df[col] = df[col].fillna('0000')
            df[col] = clean_icd_codes(df[col], remove_decimal=True)
        p2_cols = [x for x in df.columns if 'pII' in x]
        if not self.drop_p2:
            final_cols += p2_cols
        else:
            cause_cols = list(set(cause_cols) - set(p2_cols))
        final_cols += cause_cols
        df = df[final_cols]
        assert df.notnull().values.all(), f'these columns are null {df.loc[:, df.isna().any()]}'
        if self.source not in ["SOURCES"]:
            assert (df[self.value_cols].isin([1, 0])).all().all(), 'need individual records'
        return df

    def get_formatted_dataframe(self):
        df = self.df.copy()
        assert len(df[['nid', 'year_id']].drop_duplicates()) == 1,\
            "Found more than one year, this will overwrite data in format_map phase."
        self.map_extract_type_id(df)
        df = self.check_formatting(df)
        return df

    def get_diagnostic_dataframe(self):
        raise NotImplementedError

    def get_year_specific_extract_type(self, df):
        year = df['year_id'].iloc[0]
        nid = df['nid'].iloc[0]

        extract_type = f"{nid}: {year} data"
        extract_type_id = pull_or_create_extract_type(extract_type)

        df['extract_type_id'] = extract_type_id

        return df

    def map_extract_type_id(self, df, extract_type=None, conn_def='prodcod'):
        nid_extract_types = {}
        start_extract_type = extract_type
        for nid in df.nid.unique():
            if start_extract_type is None:
                extract_type = None
                extract_type = induce_extraction_type(df.loc[df['nid'] == nid], conn_def=conn_def)
            else:
                extract_type = extract_type

            extract_type = str(extract_type).strip()
            extract_type_id = pull_extract_type_id(extract_type, conn_def=conn_def)
            assert extract_type_id is not None
            nid_extract_types[nid] = extract_type_id

        df['extract_type_id'] = df['nid'].map(nid_extract_types)
        report_if_merge_fail(df, 'extract_type_id', 'nid')

        if (df['nid'].isin(self.split_nids).any()):
            df = self.get_year_specific_extract_type(df)

    def write_metadata(self, df, replace=True):
        df = df.copy()
        format_timestamp = cod_timestamp()

        for table in ['mcause_nid_metadata', 'mcause_nid_location_year']:
            if table == "mcause_nid_metadata":
                keep_cols = ['nid', 'extract_type_id', 'source', 'data_type_id', 'code_system_id',
                             'parent_nid', 'is_active']
                df['parent_nid'] = None
                df['is_active'] = 1
                df['source'] = self.source
                df['data_type_id'] = self.data_type_id
            elif table == "mcause_nid_location_year":
                keep_cols = ['nid', 'extract_type_id', 'location_id',
                             'year_id', 'representative_id']
                df['representative_id'] = 1

            idf = df[keep_cols].drop_duplicates()
            nid_extracts = df[
                ['nid', 'extract_type_id']
            ].drop_duplicates().to_records(index=False)
            for nid, extract_type_id in nid_extracts:
                nid = int(nid)
                extract_type_id = int(extract_type_id)
                idf = idf.loc[
                    (idf['nid'] == nid) &
                    (idf['extract_type_id'] == extract_type_id)
                ].copy()

                int_cols = [x for x in idf.columns if 'id' in x and 'parent_nid' not in x]
                idf[int_cols] = idf[int_cols].astype(int)

                idf["last_formatted_timestamp"] = format_timestamp

                write_to_db_nid_table(idf, table, replace=True)

        self.cache_mcause_nid_files()

    def cache_mcause_nid_files(self):
        print("\nRefreshing mcause nid metadata cache files")
        force_cache_options = {
            'force_rerun': True,
            'block_rerun': False,
            'cache_dir': "standard",
            'cache_results': True,
            'verbose': True
        }
        get_nid_metadata(**force_cache_options)
        get_nidlocyear_map(**force_cache_options)
