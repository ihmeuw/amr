import pandas as pd
import numpy as np
from warnings import warn
from cod_prep.claude.claude_io import Configurator
from mcod_prep.utils.causes import get_infsyn_hierarchy, get_all_related_syndromes
from mcod_prep.utils.mcause_io import get_mcause_data
from amr_prep.utils.amr_io import get_amr_data
from cod_prep.utils import (
    print_log_message, report_duplicates, create_square_df,
    report_if_merge_fail
)
from mcod_prep.utils.nids import add_nid_metadata, get_datasets
from cod_prep.claude.relative_rate_split import relative_rate_split
from cod_prep.downloaders import (
    add_age_metadata, get_current_location_hierarchy,
    add_location_metadata, get_cod_ages, pretty_print,
    get_pop, add_population, getcache_age_aggregate_to_detail_map,
    get_country_level_location_id, get_ages,
    prep_age_aggregate_to_detail_map
)


class PathogenFormatter():
    blank_pathogens = ['none', 'unknown']
    unique_cols = [
        'source', 'nid', 'data_type', 'location_id', 'year_id', 'age_group_id', 'sex_id',
        'hosp', 'infectious_syndrome', 'main_diagnosis', 'sample_id', 'pathogen'
    ]

    def __init__(self, model_type, infectious_syndrome, keep_pathogens,
                 age_weights_use=None, cache_kwargs=None):
        self.model_type = model_type
        self.infectious_syndrome = infectious_syndrome
        self.keep_pathogens = keep_pathogens
        self.age_weights_use = age_weights_use or {}
        if not cache_kwargs:
            self.cache_kwargs = {
                'force_rerun': False, 'block_rerun': True, 'cache_results': False
            }
        else:
            self.cache_kwargs = cache_kwargs

        self.conf = Configurator()
        self.set_sources()
        self.validate_inputs()

    def set_sources(self):
        self.amr_metadata = pd.read_csv("FILEPATH")
        self.mcause_metadata = get_datasets(
            is_active=True, data_type_id=[3, 9, 13],
            **self.cache_kwargs
        )

        if self.model_type == 'network':
            amr_source_query = "active_3a == 1"
            mcause_data_type_query = "data_type_id in (3, 9, 13)"
        elif self.model_type == 'cfr':
            amr_source_query = "active_3b == 1"
            mcause_data_type_query = "data_type_id == 3"
        self.pull_sources = set(self.amr_metadata.query(
            amr_source_query + " | inform_age_weights == 1"
        ).source)
        self.final_sources = set(
            self.amr_metadata.query(amr_source_query).source
        ).union(set(
            self.mcause_metadata.query(mcause_data_type_query).source
        ))

        self.age_weight_sources = set(self.amr_metadata.query(
            "inform_age_weights == 1"
        ).source).union(self.mcause_metadata.source)

        self.other_node_sources = set(self.amr_metadata.query(
            "inform_other_node == 1"
        ).source)

        self.deaths_only_sources = set(
            self.amr_metadata.query("deaths_only == 1").source
        ).union(
            self.mcause_metadata.query("data_type_id in (9, 13)").source
        )

        self.icu_only_sources = set(
            self.amr_metadata.query("icu_only == 1").source
        )

    def validate_inputs(self):
        assert self.conf.config_type == 'amr'
        assert self.model_type in ['network', 'cfr']

        self.infsyn = get_infsyn_hierarchy()
        assert self.infectious_syndrome in self.infsyn.infectious_syndrome.unique()
        self.target_syndromes = get_all_related_syndromes(self.infectious_syndrome, self.infsyn)

        pathogens = pd.read_csv("FILEPATH")
        assert set(self.keep_pathogens) <= set(pathogens.pathogen)
        assert set(self.age_weights_use) <= set(pathogens.pathogen).union(
            {'all', 'other'})
        assert set(self.age_weights_use.values()) <= set(pathogens.pathogen).union(
            {'all', 'other'})

    def fix_ages(self, df):
        df.loc[df.age_group_id.isin([388, 389]), 'age_group_id'] = 4
        df.loc[df.age_group_id.isin([34, 238]), 'age_group_id'] = 5
        df.loc[df.age_group_id == 325, 'age_group_id'] = 5
        df.loc[df.age_group_id.isin([247, 361, 382]), 'age_group_id'] = 183
        df.loc[df.age_group_id.isin([164, 429]), 'age_group_id'] = 2
        return df

    def get_data(self):
        print_log_message("Reading MCoD, hospital, and linkage data...")
        df_mcause = get_mcause_data(
            'format_map', is_active=True, data_type_id=[3, 9, 13],
            sub_dirs='infectious_syndrome',
            **self.cache_kwargs
        )
        df_mcause = add_nid_metadata(
            df_mcause, ['source', 'data_type_id'], **self.cache_kwargs)

        df_mcause['infectious_syndrome'] = df_mcause['infectious_syndrome'].astype(str)
        df_mcause = df_mcause.loc[df_mcause.infectious_syndrome.isin(self.target_syndromes)]

        df_mcause['pathogen_from_cause'] = df_mcause.pathogen_from_cause.str.split(',')
        df_mcause['sample_id'] = ['mcause_' + str(x) for x in range(0, len(df_mcause))]
        df_mcause = df_mcause.explode('pathogen_from_cause')
        df_mcause = df_mcause.rename(columns={
            'admissions': 'cases',
            'pathogen_from_cause': 'pathogen',
            'cause_infectious_syndrome': 'main_diagnosis'
        })
        df_mcause['cases'] = df_mcause['cases'].fillna(df_mcause['deaths'])

        print_log_message("Reading AMR data...")
        df_amr = get_amr_data(
            sources=list(self.pull_sources),
            infectious_syndrome=self.target_syndromes
        )
        df_amr = df_amr.loc[
            ~((df_amr.source == 'maternal_neonatal_sepsis_lit')
              & (df_amr.cause_id == 368))
        ]
        if self.infectious_syndrome == 'peritoneal_and_intra_abdomen_infectious':
            drop_nids = [471994, 471992, 471851, 471998, 471996, 471853, 471876]
            df_amr = df_amr.loc[~df_amr.nid.isin(drop_nids)]

        df_amr['main_diagnosis'] = 'microbiology'

        df = pd.concat([df_mcause, df_amr], sort=False)
        df['hosp'] = df['hosp'].fillna('unknown')

        nid_meta = pd.read_csv("FILEPATH")
        nid_unmapped_meta = pd.read_csv("FILEPATH")

        data_types = nid_meta.append(nid_unmapped_meta).loc[:, [
            'source', 'data_type']].drop_duplicates().reset_index()
        df = df.merge(data_types, on='source', how='left')
        df['data_type'] = df['data_type'].fillna(
            df['data_type_id'].map({9: 'mcod', 3: 'hospital', 13: 'linkage'})
        )
        no_dt_sources = list(df.loc[df['data_type'].isna(), 'source'].unique())
        assert len(no_dt_sources) == 0,\
            f"{no_dt_sources} don't have a data_type in their NID metadata, fix imminently!"
        df = df.drop('data_type_id', axis='columns')

        df = df.reset_index(drop=True)
        tabulated_data = df.sample_id.isnull()
        df.loc[tabulated_data, 'sample_id'] = df['source'] + '_'
        df.loc[tabulated_data, 'sample_id'] += [str(x) for x in range(0, tabulated_data.sum())]
        assert df[self.unique_cols].notnull().values.all()
        inconsistent_hosp = df[['source', 'sample_id', 'hosp']].drop_duplicates()\
            .duplicated(subset=['source', 'sample_id']).sum()
        warn(
            f"{inconsistent_hosp} out of {df.sample_id.nunique()} sample_ids"
            f" show inconsistent hosp"
        )
        inconsistent_counts = df[self.unique_cols + ['deaths', 'cases']].drop_duplicates()\
            .duplicated(subset=self.unique_cols, keep=False).sum()
        df = df.sort_values(by='cases', ascending=False).drop_duplicates(
            subset=self.unique_cols, keep='first'
        )
        warn(
            f"{inconsistent_counts} out of {len(df)} records show inconsistent counts"
        )
        df = df[self.unique_cols + ['deaths', 'cases']]

        df = self.fix_ages(df)
        df.loc[df.year_id == 2020, 'year_id'] = 2019
        df = df.loc[df.year_id >= 1975]
        df.loc[df.year_id < 1980, 'year_id'] = 1980

        return df

    def select_pathogens(self, df):
        assert set(PathogenFormatter.blank_pathogens).isdisjoint(self.keep_pathogens)
        df = df.loc[~df.pathogen.isin(PathogenFormatter.blank_pathogens)]

        df = df.loc[df['pathogen'] != 'contaminant']

        parents_to_children = {
            "virus": r"virus",
            "escherichia_coli": r"(enteropathogenic|enterotoxigenic)_escherichia_coli"
        }

        for parent, children in parents_to_children.items():
            df.loc[df.pathogen.str.contains(children), "pathogen"] = parent
        df = df.drop_duplicates(subset=self.unique_cols)

        if self.infectious_syndrome == 'uti_plus':
            df = df.loc[~df.pathogen.isin(['chlamydia_spp', 'neisseria_gonorrheae'])]
        if self.infectious_syndrome == 'respiratory_infectious':
            df = df.loc[df.pathogen != 'malaria']
        if self.infectious_syndrome == 'skin_infectious':
            df = df.loc[df.pathogen != 'fungus']

        if self.infectious_syndrome in ['respiratory_infectious', 'blood_stream_infectious', 'uti_plus', 'peritoneal_and_intra_abdomen_infectious', 'skin_infectious']:
            df = df.loc[(df['pathogen'] != 'coagulase_negative_staph')
                        | (df['age_group_id'].isin([2, 3])), :]

        sample_cols = [c for c in self.unique_cols if c != 'pathogen']
        df['can_add_to_poly'] = ~df.pathogen.str.contains("virus|fungus|other_unsp_bacteria")
        df = df.loc[~(
            df.duplicated(subset=sample_cols, keep=False) & ~df.can_add_to_poly &
            df.groupby(sample_cols).can_add_to_poly.transform(np.any)
        )]
        priority = {'other_unsp_bacteria': 1, 'virus': 2, 'fungus': 3}
        rows_to_prioritize = (df.duplicated(subset=sample_cols, keep=False) &
                              ~df.groupby(sample_cols).can_add_to_poly.transform(np.any))
        sep_df = df.loc[rows_to_prioritize].copy()
        df = df.loc[~rows_to_prioritize]
        sep_df['priority'] = sep_df['pathogen'].map(priority)
        assert sep_df['priority'].notnull().all()
        sep_df = sep_df.sort_values(by='priority').drop_duplicates(
            subset=sample_cols, keep='first'
        )
        sep_df = sep_df.drop('priority', axis='columns')
        df = df.append(sep_df, sort=False)

        df.loc[df.duplicated(subset=sample_cols, keep=False), 'pathogen'] = 'polymicrobial'
        df = df.drop_duplicates(subset=self.unique_cols)
        report_duplicates(df, sample_cols)

        if self.model_type == 'network':
            others = df.loc[(df.source.isin(self.other_node_sources) | df.nid.isin([468273, 468279])) & (
                ~df.pathogen.isin(self.keep_pathogens)), :].groupby('pathogen')['cases'].agg('sum').reset_index()
            others['syndrome'] = self.infectious_syndrome
            others.to_csv("FILEPATH")
            df.loc[
                (df.source.isin(self.other_node_sources.union({'gbs_meningitis_lit'})) | (
                    df.nid.isin([468273, 468279]))
                 ) & (~df.pathogen.isin(self.keep_pathogens)),
                'pathogen'
            ] = 'true_other'
            df = df.loc[df.pathogen.isin(self.keep_pathogens + ['true_other'])]
            df.loc[df.pathogen == 'true_other', 'pathogen'] = 'other'
            df = df.loc[
                (df.pathogen != 'polymicrobial') | df.source.isin(self.other_node_sources)
            ]
        elif self.model_type == 'cfr':
            df.loc[~df.pathogen.isin(self.keep_pathogens), 'pathogen'] = 'other'
        return df

    def prep_age_sex_weights(self, ref_df, pop_df, value_col):
        ref_df = ref_df.loc[ref_df[value_col].notnull()]
        ages = get_cod_ages(gbd_round_id=self.conf.get_id("gbd_round"))
        ref_df = ref_df.loc[
            ref_df.age_group_id.isin(ages.age_group_id.unique().tolist()) &
            ref_df.sex_id.isin([1, 2])
        ]
        if value_col == 'cases':
            age_weight_sources = self.age_weight_sources - self.deaths_only_sources
        else:
            age_weight_sources = self.age_weight_sources
        ref_df = ref_df.loc[ref_df.source.isin(age_weight_sources)]
        missing_ages = set(ages.age_group_id) - set(ref_df.age_group_id)
        if len(missing_ages) > 0:
            warn(
                f"The reference data for age/sex weights is missing the following"
                f" age_group_ids: {missing_ages}"
            )

        ref_df = ref_df.groupby(
            ['nid', 'location_id', 'year_id',
             'age_group_id', 'sex_id', 'pathogen'], as_index=False
        )[value_col].sum()

        ref_df = add_population(ref_df, pop_df=pop_df)
        report_if_merge_fail(
            ref_df, 'population', ['location_id', 'age_group_id', 'sex_id', 'year_id'])
        weight_cols = ['age_group_id', 'sex_id', 'pathogen']
        ref_df = ref_df.groupby(weight_cols, as_index=False)[
            [value_col, 'population']].sum()
        ref_df = ref_df.append(
            ref_df.groupby(['age_group_id', 'sex_id'], as_index=False)
            [[value_col, 'population']].sum().assign(pathogen='all')
        )
        ref_df['weight'] = ref_df[value_col] / ref_df['population']
        ref_df_sq = create_square_df(
            ref_df.append(pd.DataFrame(missing_ages, columns=['age_group_id'])),
            weight_cols
        )
        ref_df = ref_df_sq.merge(
            ref_df, how='left', on=weight_cols, validate='one_to_one').fillna(0)
        ref_df = ref_df[weight_cols + ['weight']]
        return ref_df

    def add_country_location_id(self, df):
        country_locs = get_country_level_location_id(
            df.location_id.unique().tolist(),
            get_current_location_hierarchy(
                location_set_version_id=self.conf.get_id("location_set_version"),
                **self.cache_kwargs
            )
        )
        df = df.merge(country_locs, how='left', on='location_id', validate='many_to_one')
        report_if_merge_fail(df, 'country_location_id', 'location_id')
        return df

    def age_sex_split(self, df, pop_df, weights, value_col):
        pop_id_cols = ['location_id', 'age_group_id', 'sex_id', 'year_id']

        print_log_message("    Prepping age map")
        age_detail_map = getcache_age_aggregate_to_detail_map(
            gbd_round_id=self.conf.get_id("gbd_round"), **self.cache_kwargs
        )
        sex_detail_map = pd.DataFrame(
            columns=['agg_sex_id', 'sex_id'],
            data=[
                [3, 1],
                [3, 2],
                [9, 1],
                [9, 2],
                [1, 1],
                [2, 2]
            ]
        )
        detail_maps = {
            'age_group_id': age_detail_map,
            'sex_id': sex_detail_map
        }

        print_log_message("    Prep which distributions should be used")
        pathogen_to_weight_pathogen_map = df[['pathogen']].drop_duplicates()\
            .assign(dist_pathogen=lambda d: d['pathogen'])

        pathogen_to_weight_pathogen_map.loc[
            ~pathogen_to_weight_pathogen_map.pathogen.isin(weights.pathogen),
            'dist_pathogen'
        ] = 'all'
        pathogen_to_weight_pathogen_map['dist_pathogen'].update(
            pathogen_to_weight_pathogen_map['pathogen'].map(
                self.age_weights_use
            )
        )
        self.pathogen_to_weight_pathogen_map = pathogen_to_weight_pathogen_map
        val_to_dist_maps = {
            'pathogen': pathogen_to_weight_pathogen_map
        }
        split_cols = ['age_group_id', 'sex_id']
        split_inform_cols = ['pathogen']
        value_cols = [value_col]
        start_val = df[value_col].sum()
        start_cols = df.columns.tolist()

        df = self.add_country_location_id(df)
        df['orig_location_id'] = df['location_id']
        df['location_id'] = df['country_location_id']

        print_log_message("    Running split")
        df = relative_rate_split(
            df,
            pop_df,
            weights,
            detail_maps,
            split_cols,
            split_inform_cols,
            pop_id_cols,
            value_cols,
            pop_val_name='population',
            val_to_dist_map_dict=val_to_dist_maps,
            verbose=False
        )
        assert np.isclose(start_val, df[value_col].sum())
        return df[start_cols]

    def split_one_metric(self, df, pop_df, weights, value_col):
        keep_cols = self.unique_cols + [value_col]
        df = df[keep_cols].copy()
        df = df.loc[df[value_col].notnull()]
        if value_col == 'cases':
            df = df.loc[~df.source.isin(self.deaths_only_sources)]
        df = self.age_sex_split(
            df, pop_df, weights, value_col
        )
        self.dem_cols = [c for c in self.unique_cols if c != 'sample_id']
        df = df.groupby(self.dem_cols, as_index=False)[value_col].sum()
        return df

    def format_data(self):
        print_log_message(f"Getting data for model type {self.model_type}")
        df = self.get_data()

        print_log_message(f"Selecting specified pathogens")
        df = self.select_pathogens(df)

        print_log_message(f"Generating age/sex splitting weights")
        pop_df = get_pop(
            pop_run_id=self.conf.get_id("pop_run"), **self.cache_kwargs
        )
        self.death_weights = self.prep_age_sex_weights(
            df.copy(), pop_df, value_col='deaths'
        )
        self.case_weights = self.prep_age_sex_weights(
            df.copy(), pop_df, value_col='cases'
        )

        print_log_message("Splitting deaths and cases")
        if self.model_type == 'cfr':
            df = df.loc[df.deaths.notnull() & df.cases.notnull()]
        if self.infectious_syndrome == 'cns_infectious':
            gbs_neonates = df.loc[
                (df.pathogen == 'group_b_strep') & (df.age_group_id.isin([2, 3, 42])), 'nid'
            ].unique().tolist()
        deaths_df = self.split_one_metric(
            df, pop_df, self.death_weights, value_col='deaths'
        )
        cases_df = self.split_one_metric(
            df, pop_df, self.case_weights, value_col='cases'
        )

        print_log_message("Recombining data frames")
        df = pd.merge(
            deaths_df, cases_df, how='outer',
            on=self.dem_cols,
            validate='one_to_one'
        )

        print_log_message(
            f"Subsetting to final data needed for {self.model_type}"
        )
        df = df.loc[df.source.isin(self.final_sources)]
        if self.model_type == 'cfr':
            cfrs = df.groupby('source').apply(lambda x: pd.Series(
                {'cfr': sum(x['deaths']) / sum(x['cases'])})).reset_index()
            sources_to_drop = cfrs.loc[cfrs['cfr'].isin([0, 1]), 'source'].tolist()
            df = df.loc[~df['source'].isin(sources_to_drop), :]
            df.loc[df.source.isin(self.icu_only_sources), 'ICU'] = 'ICU_only'
            df.loc[~df.source.isin(self.icu_only_sources), 'ICU'] = 'mixed'
        elif self.model_type == 'network':
            if self.infectious_syndrome == 'cns_infectious':
                df = df.loc[~(
                    ~df.nid.isin(gbs_neonates)
                    & (df.pathogen == 'group_b_strep')
                    & (df.age_group_id.isin([2, 3]))
                )]
            df.loc[df.source.isin(self.icu_only_sources), 'cases'] = np.NaN
            assert df.main_diagnosis.notnull().all()
            assert df.infectious_syndrome.notnull().all()
            community = df.infectious_syndrome == df.main_diagnosis
            microbio = df.main_diagnosis == "microbiology"
            df.loc[community & ~microbio, 'hosp'] = 'community'
            df.loc[~community & ~microbio, 'hosp'] = 'hospital'
        return df
