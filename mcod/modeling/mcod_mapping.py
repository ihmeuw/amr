from builtins import zip
import pandas as pd
import re
from cod_prep.utils import (
    print_log_message, report_duplicates, clean_icd_codes, report_if_merge_fail
)
from cod_prep.downloaders import (
    get_cause_map, add_code_metadata, add_age_metadata,
    add_cause_metadata, get_cause_package_map
)
from cod_prep.claude.configurator import Configurator
from mcod_prep.utils.causes import (
    get_infsyn_hierarchy, get_child_to_available_parent_syndrome,
    get_all_related_syndromes
)


class MCoDMapper():
    cache_options = {'force_rerun': False, 'block_rerun': True}
    conf = Configurator()
    int_cause_name_dict = {
        'sepsis': '',
    }
    infectious_syndromes = [
        # Special option that preps data for all syndromes rather than for a
        # binomial regression on one syndrome
        'infectious_syndrome',
        # Level 1 infectious syndromes
        'blood_stream_infectious',
        'skin_infectious',
        'cns_infectious',
        'diarrhea',
        'peritoneal_and_intra_abdomen_infectious',
        'bone_joint_infection',
        'cardiac_infectious',
        'tb',
        'typhoid_paratyphoid_ints',
        'chlamydia_and_gonorrheae',
        'others_and_non_bacterial_infectious',
        'uti_plus_hosp',
        'uti_plus_comm',
        'respiratory_infectious_hosp',
        'respiratory_infectious_comm',
    ]
    possible_int_causes = list(int_cause_name_dict.keys()) + infectious_syndromes
    # Remainder syndrome
    remainder_syndromes = [
        'unspecified_and_none', 'non_infectious_cause', 'inj_medical', 'other_injuries'
    ]
    # Placeholder code for nulls
    null_code = '0000'

    def __init__(self, int_cause, drop_p2, code_system_id=None, code_map_version_id=None,
                 path_to_ucause_map=None, path_to_int_cause_map=None):
        self.int_cause = int_cause
        self.drop_p2 = drop_p2

        self.code_system_id = code_system_id
        self.code_map_version_id = code_map_version_id
        self.path_to_ucause_map = path_to_ucause_map
        self.path_to_int_cause_map = path_to_int_cause_map

        if self.code_system_id and self.code_map_version_id:
            self.is_real_code_system = True
        elif self.path_to_ucause_map or self.path_to_int_cause_map:
            assert self.int_cause == 'infectious_syndrome',\
                "You specified paths to custom maps but you are not trying"\
                "to capture infectious syndromes - this option is not implemented"
            self.is_real_code_system = False
        else:
            raise AssertionError(
                "You must pass either code system and code map version "
                "or paths to ucause and int_cause maps"
            )

        assert self.int_cause in self.possible_int_causes, \
            f"{self.int_cause} is not a valid intermediate cause"
        try:
            self.full_cause_name = self.int_cause_name_dict[self.int_cause]
        except KeyError:
            self.full_cause_name = self.int_cause
        if self.full_cause_name == '':
            self.full_cause_name = self.int_cause
        if type(self.full_cause_name) != list:
            self.full_cause_name = [self.full_cause_name]

    @staticmethod
    def get_code_columns(df):
        """Get a list of raw cause columns with ICD codes as values."""
        col_names = list(df.columns)
        code_cols = [x for x in col_names if "multiple_cause" in x and "pII" not in x]
        if 'cause' in df:
            code_cols += ["cause"]
        return code_cols

    @staticmethod
    def _get_cause_num(mcod_col):
        """Get sort order for cause columns."""
        if mcod_col.startswith('cause'):
            return '0'
        else:
            assert re.match(r"^multiple_cause_[a-z]*[0-9]*", mcod_col), \
                f"column {mcod_col} does not match expected format: multiple_cause_x"
            return mcod_col.split('_')[2]

    @staticmethod
    def prep_raw_mapped_cause_dictionary(raw_cols, mapped_cols):
        """
        Create dictionary of raw cause columns to mapped cause columns.
        """
        raw_cols = sorted(raw_cols, key=MCoDMapper._get_cause_num)
        mapped_cols = sorted(mapped_cols, key=MCoDMapper._get_cause_num)
        return dict(list(zip(raw_cols, mapped_cols)))

    @staticmethod
    def fix_icd_codes(df, codes, code_system_id):
        """Adjustment to icd9/10 cause codes."""
        if code_system_id == 6:
            df.loc[df['cause'].str.contains('^[89]'), 'cause'] = 'E' + df['cause']
        return df

    @staticmethod
    def prep_cause_package_map(cause_package_map):
        check_map = cause_package_map[['map_id', 'map_type']].drop_duplicates()
        report_duplicates(check_map, 'map_id')
        return cause_package_map.set_index('value')['map_id'].to_dict()

    @staticmethod
    def prep_cause_map(cause_map):
        """Clean up cause map."""
        cause_map['value'] = clean_icd_codes(cause_map['value'], remove_decimal=True)
        cause_map = cause_map.drop_duplicates(['code_system_id', 'value'])
        cause_map['code_id'] = cause_map['code_id'].astype(int)
        return cause_map.set_index('value')['code_id'].to_dict()

    @staticmethod
    def map_cause_codes(df, coi_map, coi, cols_to_map=None):
        """Map cause codes to any given value (e.g. acause, category, etc.).

        Inputs
        df (pd dataframe): incoming, unmapped data with ICD codes
        coi_map (pd dataframe): special map designed just for one cause of interest
        coi (string): cause of interest
        Returns
        df (pd dataframe): mapped dataframe with additional columns for each cause
        """
        df = df.copy()
        if not cols_to_map:
            cols_to_map = MCoDMapper.get_code_columns(df)
        # map chain causes using cause of interest map
        for col in cols_to_map:
            df[col] = df[col].fillna(MCoDMapper.null_code)
            df[col] = df[col].astype(object)
            df[col + '_' + coi] = df[col].map(coi_map)
        return df

    @staticmethod
    def trim_and_remap(df, code_dict, cause_map):
        """Trim ICD codes and map again."""
        df = df.copy()
        # before trimming, map "null" chain causes to MCoDMapper.null_code
        for code, mapped_code in list(code_dict.items()):
            df.loc[df[code] == MCoDMapper.null_code, mapped_code] = MCoDMapper.null_code

        # trim and re map null mappings
        for n in reversed(range(3, 6)):
            for code, mapped_code in list(code_dict.items()):
                temp_code = 'temp_' + code
                df[temp_code] = df[code].copy()
                try:
                    df.loc[df[mapped_code].isnull(), temp_code] = df[temp_code].apply(
                        lambda x: x[0:n])
                except TypeError:
                    if mapped_code != 'cause_mapped':
                        df[mapped_code] = MCoDMapper.null_code
                    else:
                        print("problem code here..." + df[code])
                df.loc[df[mapped_code].isnull(), mapped_code] = df[temp_code].map(cause_map)
                df = df.drop(temp_code, axis=1)
        return df

    def prep_int_cause_map(self):
        map_dir = self.conf.get_directory('process_inputs')
        code_system_name = {1: 'icd10', 6: 'icd9'}.get(self.code_system_id, None)
        if self.int_cause == 'sepsis':
            df = pd.read_csv(
                "FILEPATH",
                dtype={'icd_code': object})
            df['icd_code'] = clean_icd_codes(df['icd_code'], remove_decimal=True)
            df = df[['icd_code', 'sepsis']].drop_duplicates()
            report_duplicates(df, 'icd_code')
            mcod_map = dict(list(zip(df['icd_code'], df['sepsis'])))
            mcod_map.update({
                'acause_inj_trans_road_4wheel': 'no_sepsis',
                'acause_inj_fires': 'no_sepsis',
                'acause_inj_poisoning_other': 'no_sepsis'})
            assert set(mcod_map.values()).issubset({
                'ucod_no need of-implicit',
                'explicit',
                'ucod_need of-implicit',
                'no_sepsis',
                'implicit-of-code'})

        elif self.int_cause in self.infectious_syndromes:
            if self.is_real_code_system:
                df = pd.read_csv(
                    "FILEPATH",
                    dtype={'icd_code': object})
                if self.code_system_id in [1, 6]:
                    df['icd_code'] = clean_icd_codes(df['icd_code'], remove_decimal=True)
                df = df[['icd_code', 'infectious_syndrome']].drop_duplicates()
                report_duplicates(df, 'icd_code')
                mcod_map = dict(list(zip(df['icd_code'], df['infectious_syndrome'])))
            else:
                df = pd.read_csv(self.path_to_int_cause_map)[['cause', 'infectious_syndrome']]
                report_duplicates(df, 'cause')
                mcod_map = dict(list(zip(df['cause'], df['infectious_syndrome'])))

        return mcod_map

    def capture_int_cause(self, df, int_cause_cols):
        """Flag deaths related to the intermediate cause."""
        if self.int_cause == 'sepsis':
            df[int_cause_cols] = df[int_cause_cols].fillna("no_sepsis")
            # First set explicit sepsis
            df.loc[(df[int_cause_cols] == 'explicit').any(axis=1), 'sepsis'] = 'explicit'
            # Next set implicit sepsis
            not_explicit = df.sepsis != 'explicit'
            ucod_no_need_of = df.cause_sepsis == 'ucod_no need of-implicit'
            ucod_need_of = df.cause_sepsis == 'ucod_need of-implicit'
            has_of = (df[int_cause_cols] == 'implicit-of-code').any(axis=1)
            df.loc[not_explicit & (ucod_no_need_of | (ucod_need_of & has_of)),
                   'sepsis'] = 'implicit'
            df['sepsis'] = df['sepsis'].fillna('no_sepsis')
        elif self.int_cause in self.infectious_syndromes:
            df = self.apply_syndrome_severity(df, int_cause_cols, self.int_cause)
        assert df[self.int_cause].notnull().values.all()
        return df

    def apply_syndrome_severity(self, df, int_cause_cols: list, syndrome: str):
        assert {
            "_row_id", "pathogen_from_cause",
            'infectious_syndrome', 'l1_syndrome'
        }.isdisjoint(set(df))
        infsyn = get_infsyn_hierarchy()
        severity_map = infsyn.set_index("infectious_syndrome")['severity'].to_dict()
        has_ucause = 'cause' in df
        if has_ucause:
            cause_col_dict = {
                MCoDMapper._get_cause_num(cause_col): cause_col for cause_col in int_cause_cols
            }
            ucod_col = cause_col_dict['0']  # Underlying cause column mapped to syndrome
        true_syndromes = list(set(severity_map.keys()) - set(MCoDMapper.remainder_syndromes))

        # First, determine the syndrome responsible for each case
        # based on severity
        reverse_severity = {v: k for k, v in severity_map.items()}
        df = df.assign(
            infectious_syndrome=lambda x: x[int_cause_cols].fillna('unspecified_and_none').apply(
                lambda x: x.map(severity_map)
            ).min(axis=1).map(reverse_severity)
        )
        assert df['infectious_syndrome'].notnull().all()
        # In neonates, BSI should override everything
        bsi = get_all_related_syndromes('blood_stream_infectious', infsyn=infsyn)
        neonatal_bsi = get_all_related_syndromes('bsi_neonat_sp_bactria', infsyn=infsyn)
        df = add_age_metadata(df, ['age_group_days_end'], **self.cache_options)
        report_if_merge_fail(df, 'age_group_days_end', 'age_group_id')
        neonatal = df.age_group_days_end <= 27
        has_bsi = df[int_cause_cols].isin(bsi).any(axis=1)
        has_neonatal_bsi = df[int_cause_cols].isin(neonatal_bsi).any(axis=1)
        df.loc[(neonatal & has_bsi) | has_neonatal_bsi, 'infectious_syndrome'] =\
            df[int_cause_cols]\
            .apply(lambda x: list(set(x).intersection(bsi)), axis=1)\
            .apply(
                lambda x: reverse_severity[min([severity_map[y] for y in x])]
                if len(x) > 0 else ""
            )
        df = df.drop('age_group_days_end', axis='columns')

        if has_ucause:
            # Where the ucod is an infectious syndrome, override syndrome severity
            df.loc[df[ucod_col].isin(true_syndromes), 'infectious_syndrome'] = df[ucod_col]

        # Gather all the syndromes within the level 1 syndrome and
        # extract any pathogens
        syndrome_to_level_1_map = get_child_to_available_parent_syndrome(
            infsyn.query("level == 1").infectious_syndrome.unique(),
            infsyn=infsyn
        )
        level_1_to_children = {
            l1: get_all_related_syndromes(l1, infsyn=infsyn)
            for l1 in infsyn.query("level == 1").infectious_syndrome.unique()
        }
        df['l1_syndrome'] = df['infectious_syndrome'].map(syndrome_to_level_1_map)
        assert df['l1_syndrome'].notnull().all()
        df['_row_id'] = list(range(0, df.shape[0]))
        syn_df = df.loc[:, ['_row_id', 'l1_syndrome'] + int_cause_cols]\
            .fillna('unspecified_and_none')\
            .melt(
                id_vars=['_row_id', 'l1_syndrome'],
                value_vars=int_cause_cols, value_name='detailed_infsyn',
                var_name='chain_col')\
            .loc[:, ['_row_id', 'l1_syndrome', 'detailed_infsyn']]\
            .drop_duplicates()
        syn_df = pd.concat([
            group.loc[group.detailed_infsyn.isin(level_1_to_children[key])]
            for key, group in syn_df.groupby(['l1_syndrome'])])
        syndrome_to_pathogen = infsyn.set_index('infectious_syndrome')['pathogen'].to_dict()
        syn_df['pathogen_from_cause'] = syn_df['detailed_infsyn'].map(
            syndrome_to_pathogen
        )

        sample_to_pathogens = syn_df.loc[
            syn_df.pathogen_from_cause.notnull()
        ].groupby('_row_id')['pathogen_from_cause'].aggregate(
            lambda x: ','.join(x.drop_duplicates().sort_values())
        ).to_dict()
        df['pathogen_from_cause'] = df['_row_id'].map(sample_to_pathogens).fillna("none")

        # Finally, create a flag column indicating whether the case
        # is associated with the syndrome of interest, if applicable
        if self.int_cause != 'infectious_syndrome':
            # Get a list containing the syndrome of interest plus
            # any of its children
            if "_hosp" in self.int_cause:
                target_syndromes = get_all_related_syndromes(
                    self.int_cause.replace("_hosp", ""), infsyn)

            if "_comm" in self.int_cause:
                target_syndromes = get_all_related_syndromes(
                    self.int_cause.replace("_comm", ""), infsyn)

            if "_comm" not in self.int_cause and "_hosp" not in self.int_cause:
                target_syndromes = get_all_related_syndromes(self.int_cause, infsyn)

            # Find all samples that have these syndromes
            samples_to_flag = syn_df.loc[
                syn_df.detailed_infsyn.isin(target_syndromes), '_row_id'
            ].unique()
            df[self.int_cause] = df._row_id.isin(samples_to_flag).astype(int)

            # Adjusting flag for HAI vs CAI
            if "_comm" in self.int_cause:
                df.loc[(df['infectious_syndrome'] != df[ucod_col]) &
                       (df[self.int_cause] == 1), self.int_cause] = 0

            if "_hosp" in self.int_cause:
                df.loc[(df['infectious_syndrome'] == df[ucod_col]) &
                       (df[self.int_cause] == 1), self.int_cause] = 0

        df = df.drop(['_row_id', 'l1_syndrome'], axis='columns')

        return df

    def set_part2_flag(self, df):
        """Mark whether or not the int cause is from part 2 of the death certificate."""
        p2_cols = [x for x in df.columns if 'pII' in x]
        int_cause_chains = [x for x in df.columns if (self.int_cause in x) and ('multiple' in x)]
        p2_chain_dict = dict(list(zip(p2_cols, int_cause_chains)))
        df['pII_' + self.int_cause] = 0
        for p2_col, chain in sorted(p2_chain_dict.items()):
            df.loc[
                (df[chain].isin(self.full_cause_name)) &
                (df[p2_col] == 1), 'pII_' + self.int_cause
            ] = 1
        return df

    def get_computed_dataframe(self, df, map_underlying_cause=True):
        """
        Return mapped dataframe.
        """
        raw_cause_cols = MCoDMapper.get_code_columns(df)
        df = MCoDMapper.fix_icd_codes(df, raw_cause_cols, self.code_system_id)

        if map_underlying_cause:
            print_log_message("Mapping underlying cause/primary diagnosis")
            if self.is_real_code_system:
                cause_map = get_cause_map(
                    code_map_version_id=self.code_map_version_id, **self.cache_options
                )
                code_map = MCoDMapper.prep_cause_map(cause_map)
                df['cause_mapped'] = df['cause'].map(code_map)

                if self.code_system_id in [1, 6]:
                    print_log_message(
                        "Trimming ICD codes and remapping underlying cause/primary diagnosis")
                    df = MCoDMapper.trim_and_remap(
                        df, {'cause': 'cause_mapped'}, code_map)

                # merge on the cause_id for the underlying cause
                df = df.rename(columns={'cause_mapped': 'code_id'})
                df = add_code_metadata(df, 'cause_id', code_map_version_id=self.code_map_version_id,
                                       **self.cache_options)
                df.loc[df.cause == 'none', 'cause_id'] = 919
                report_if_merge_fail(df, 'cause_id', 'cause')
            else:
                cause_map = pd.read_csv(self.path_to_ucause_map)[['cause', 'cause_id']]
                df = df.merge(cause_map, how='left', on='cause', validate='many_to_one')
                report_if_merge_fail(df, 'cause_id', 'cause')

        print_log_message("Mapping chain causes")
        int_cause_map = self.prep_int_cause_map()
        df = MCoDMapper.map_cause_codes(df, int_cause_map, self.int_cause)
        # These are the "multiple?_cause_X_{int_cause}" columns from mapping
        # columns should include every cause listed on death certificate/hospital record
        int_cause_cols = [x for x in df.columns if self.int_cause in x]

        if self.code_system_id in [1, 6]:
            print_log_message("Trimming ICD codes and remapping chain causes")
            int_cause_col_dict = MCoDMapper.prep_raw_mapped_cause_dictionary(
                raw_cause_cols, int_cause_cols)
            df = MCoDMapper.trim_and_remap(df, int_cause_col_dict, int_cause_map)

        if self.int_cause == 'sepsis':
            report_if_merge_fail(df, 'cause_' + self.int_cause, 'cause')

        print_log_message("Identifying rows with intermediate cause of interest")
        df = self.capture_int_cause(df, int_cause_cols)
        if not self.drop_p2:
            df = self.set_part2_flag(df)

        # drop all of the mapped multiple causes except the underlying
        drop_cols = ['_'.join([x, self.int_cause]) for x in raw_cause_cols if 'multiple' in x]
        drop_cols = [c for c in drop_cols if c in df]
        print_log_message(f"Dropping columns: {drop_cols}")
        df = df.drop(drop_cols, axis=1)

        return df
