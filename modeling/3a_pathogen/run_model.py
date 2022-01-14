"""
Run network meta-analysis to model the distribution of pathogens
causing a given infectious syndrome.
"""
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import sys
import getpass
import yaml
import warnings
from scipy.stats import t
import pickle
from cod_prep.utils import (
    report_duplicates, report_if_merge_fail, create_square_df,
    wait_for_job_ids
)
import itertools
from cod_prep.claude.claude_io import Configurator, makedirs_safely
from mcod_prep.utils.covariates import merge_covariate
from cod_prep.downloaders import (
    add_age_metadata, get_current_location_hierarchy,
    add_location_metadata, get_cod_ages, pretty_print,
    get_pop, add_population, getcache_age_aggregate_to_detail_map,
    get_country_level_location_id, get_ages,
    prep_age_aggregate_to_detail_map
)
from mcod_prep.utils.causes import get_infsyn_hierarchy, get_all_related_syndromes
from cod_prep.utils import print_log_message, warn_if_merge_fail
from amr_prep.utils.pathogen_formatting import PathogenFormatter
from mcod_prep.utils.mcod_cluster_tools import submit_mcod
from multiprocessing import Pool
from functools import partial
sys.path.append('FILEPATH')
from netprop.data import Data
from netprop.dorm_model import DormModel
from netprop.model import Model


CONF = Configurator()


class PathogenNetwork(object):
    """Model a pathogen network"""
    blank_pathogens = ['none', 'unknown']
    out_dir = Path(CONF.get_directory("process_data").format(model_step='03a_pathogen'))
    id_cols = ['location_id', 'year_id', 'age_group_id', 'sex_id', 'hosp']

    def __init__(self, model_version, infectious_syndrome, keep_pathogens,
                 covariates, agg_age_group_ids, cfr_use, age_weights_use,
                 year_ids, cfr_ready=None, factor_covs=None, ref_pathogen=None,
                 aggregate_cols=None, unknown=0.5, age_start=None,
                 age_end=None, study_weights=None, gprior_sd=None):
        self.model_version = model_version
        self.infectious_syndrome = infectious_syndrome
        self.keep_pathogens = keep_pathogens
        self.ref_pathogen = ref_pathogen
        self.covariates = covariates
        self.agg_age_group_ids = agg_age_group_ids
        self.cfr_use = cfr_use or {}
        self.age_weights_use = age_weights_use
        self.year_ids = year_ids
        self.factor_covs = factor_covs or []
        self.cfr_ready = cfr_ready
        self.aggregate_cols = aggregate_cols or []
        self.unknown = unknown
        self.age_start = age_start
        self.age_end = age_end
        self.study_weights = study_weights
        self.gprior_sd = gprior_sd

        self.conf = Configurator()
        self.cache_kwargs = {'force_rerun': True, 'block_rerun': False, 'cache_results': True}
        self.validate_inputs()
        self.model_dir = "FILEPATH"
        makedirs_safely(str(self.model_dir))

    def validate_inputs(self):
        """Validate inputs"""
        assert self.conf.config_type == 'amr'

        self.infsyn = get_infsyn_hierarchy()
        assert self.infectious_syndrome in self.infsyn.infectious_syndrome.unique()
        self.target_syndromes = get_all_related_syndromes(self.infectious_syndrome, self.infsyn)

        pathogens = pd.read_csv("FILEPATH")
        assert set(self.keep_pathogens) <= set(pathogens.pathogen)
        if self.ref_pathogen is not None:
            assert self.ref_pathogen in pathogens.pathogen.unique()
            assert self.ref_pathogen in self.keep_pathogens

        ages = get_ages(**self.cache_kwargs)
        assert set(self.agg_age_group_ids) <= set(ages.age_group_id)

    def drop_data(self, df):
        # DATA DROPS/OUTLIERS
        if self.infectious_syndrome == 'respiratory_infectious':
            df = self.add_covariates(
                df, covariates=['microbiology'],
                add_model_age=True)
            df = df.loc[~(
                (df.pathogen == 'streptococcus_pneumoniae') & (
                    (df.microbiology == 1) | (df.source.str.endswith("_lit"))
                )
            )]
            df['total_cases'] = df.groupby(
                ['pathogen', 'agg_age_group_id', 'hosp'])['cases'].transform(sum)
            df = df.loc[df['total_cases'] >= 10]
            df = df.drop(['microbiology', 'total_cases', 'agg_age_group_id'], axis='columns')
        elif self.infectious_syndrome == 'cns_infectious':
            df = df.loc[~(
                (df.source == 'SOURCE') & (df.pathogen == 'other')
            )]
            if 'for_priors' in self.model_version:
                pass
            else:
                df = df.loc[df.source != 'SOURCE']
            if 'non_neonatal' in self.model_version:
                df = df.loc[df.source != 'SOURCE']
        elif self.infectious_syndrome == 'skin_infectious':
            df = self.add_covariates(df, covariates=['microbiology'], add_model_age=False)
            df = df.loc[
                (df['microbiology'] != 1) & (~df['nid'].isin(["SOURCE"]))]
            df = df.drop('microbiology', axis='columns')
        elif self.infectious_syndrome == 'uti_plus':
            df = df.loc[~(
                (df.source == 'SOURCE') & (df.pathogen == 'other')
            )]
        elif self.infectious_syndrome == 'peritoneal_and_intra_abdomen_infectious':
            df = df.loc[~df.nid.isin(["SOURCE"])]
            df = df.loc[df.pathogen != 'anaerobe']
            self.keep_pathogens.remove('anaerobe')
        return df

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

    def get_cfrs(self):
        """Get CFR models for pathogen/syndromes"""
        self.cfr_cols = self.id_cols + ['pathogen']
        value_col = 'predict'
        if self.cfr_ready:
            input_dir = "FILEPATH"
            df = pd.concat(
                [pd.read_csv(file) for file in input_dir.iterdir()],
                sort=False
            )
            assert df.notnull().values.all()
            if 'hosp' not in df:
                self.cfr_cols.remove('hosp')
            if (df.sex_id == 3).all():
                self.cfr_cols.remove('sex_id')

            if 'hosp' in df and 'unknown' not in df['hosp'].unique():
                df = df.append(
                    df.assign(hosp='unknown')
                    .groupby(self.cfr_cols, as_index=False)[value_col].mean()
                )
            if 'all' not in df.pathogen.unique():
                df = df.append(
                    df.query("pathogen == 'other'").assign(pathogen='all')
                )

            for pathogen, use_pathogen in self.cfr_use.items():
                df = df.loc[df.pathogen != pathogen]
                df = df.append(
                    df.query(f"pathogen == '{use_pathogen}'")
                    .assign(pathogen=pathogen)
                )

            if self.infectious_syndrome == 'urogenital_infectious':
                assert len(df.loc[(df.age_group_id == 26) & (df.hosp == 'community')]) == 0
                df = df.append(
                    df.query("age_group_id == 25 & hosp == 'community'")
                    .assign(age_group_id=26)
                )

            df = df.drop_duplicates()
            report_duplicates(df, self.cfr_cols)
            df = df.rename(columns={value_col: 'cfr'})
            df = df[self.cfr_cols + ['cfr']]
            self.cfr = df.copy()
        else:
            self.cfr = pd.DataFrame(columns=self.cfr_cols + ['cfr'])

    def subset_ages(self, df):
        df = add_age_metadata(df, ['simple_age'], **self.cache_kwargs)
        report_if_merge_fail(df, 'simple_age', 'age_group_id')
        if self.age_start is not None:
            df = df.loc[df.simple_age >= self.age_start]
        if self.age_end is not None:
            df = df.loc[df.simple_age <= self.age_end]
        df = df.drop('simple_age', axis='columns')
        return df

    def add_model_ages(self, df, allow_fail=False):
        age_group_table = get_ages()
        good_age_group_ids = df.age_group_id.unique().tolist()
        age_detail_map = prep_age_aggregate_to_detail_map(age_group_table, good_age_group_ids)\
            .query("agg_age_group_id in @self.agg_age_group_ids")
        df = df.merge(
            age_detail_map[['age_group_id', 'agg_age_group_id']], how='left',
            on='age_group_id', validate='many_to_one'
        )
        if allow_fail:
            warn_if_merge_fail(df, 'agg_age_group_id', 'age_group_id')
        else:
            report_if_merge_fail(df, 'agg_age_group_id', 'age_group_id')
        return df

    def apply_cfrs(self, df):
        # Apply CFRs whereever admissions are null
        cfr = self.cfr.rename(columns={
            'location_id': 'country_location_id',
            'age_group_id': 'agg_age_group_id'
        })
        df = self.add_model_ages(df, allow_fail=False)
        df = self.add_country_location_id(df)
        merge_cols = [
            c for c in self.cfr_cols if c not in ['location_id', 'age_group_id']
        ] + ['country_location_id', 'agg_age_group_id']
        df = df.merge(
            cfr[merge_cols + ['cfr']], how='left', validate='many_to_one',
            on=merge_cols
        )
        print_log_message(f"Filling {df.loc[df.cfr.isnull(), merge_cols].drop_duplicates()}")
        df = df.merge(
            cfr.loc[cfr.pathogen == 'all'].drop('pathogen', axis='columns')[
                [c for c in merge_cols if c != 'pathogen'] + ['cfr']
            ], how='left', validate='many_to_one',
            on=[c for c in merge_cols if c != 'pathogen']
        )
        df['cfr'] = df['cfr_x'].fillna(df['cfr_y'])
        if self.cfr_ready:
            report_if_merge_fail(df, 'cfr', merge_cols)
        else:
            df['cfr'] = df['cfr'].fillna(1)
        df['cases'] = df['cases'].fillna(df['deaths'] / df['cfr'])
        assert df['cases'].notnull().all()
        df = df.drop(
            ['agg_age_group_id', 'country_location_id', 'cfr',
             'cfr_x', 'cfr_y'], axis='columns'
        )

        df = df.groupby(
            self.id_cols + ['source', 'nid', 'pathogen'],
            as_index=False
        )['cases'].sum()
        return df

    def aggregate_data(self, df, cols):
        can_aggregate = {'age_group_id', 'sex_id', 'hosp'}
        assert set(cols) <= can_aggregate,\
            f"Don't know how to aggregate {set(cols) - can_aggregate}"
        # Aggregate age
        if 'age_group_id' in cols:
            df = self.add_model_ages(df, allow_fail=False)
            df['age_group_id'] = df['agg_age_group_id']
        if 'sex_id' in cols:
            df['sex_id'] = 3
        if 'hosp' in cols:
            df['hosp'] = 'all'
        df = df.groupby(
            self.id_cols + ['source', 'nid', 'pathogen'],
            as_index=False
        )['cases'].sum()
        return df

    def calculate_se(self, df):
        df = df.groupby(
            self.id_cols + ['source', 'nid', 'pathogen'],
            as_index=False
        )['cases'].sum()
        df['total_all_paths'] = df.groupby(
            [c for c in self.id_cols if c != 'pathogen'] + ['nid']
        )['cases'].transform(sum)
        df = df.loc[df.total_all_paths != 0]
        df = df.loc[df.cases != 0]
        df['prop'] = df['cases'] / df['total_all_paths']
        df['se'] = np.sqrt(df['prop'] * (1 - df['prop']) * df['total_all_paths'])
        assert df.notnull().values.all()
        df = df.drop(['total_all_paths', 'prop'], axis='columns')
        return df

    def set_composites(self, df):
        """
        Set composite pathogens for specific studies.
        """
        if self.infectious_syndrome == 'skin_infectious':
            df.loc[(df['nid'] == "SOURCE") & (df['pathogen'] == 'other'), 'pathogen']\
                = 'other-acinetobacter_baumanii-enterococcus_faecalis-enterococcus_spp'
        elif self.infectious_syndrome == 'cns_infectious':
            df.loc[
                (df.source == 'SOURCE') & (df.pathogen == 'other'),
                'pathogen'
            ] = '-'.join(set(self.keep_pathogens).union({'other'}) - {'group_b_strep', 'virus'})
        return df

    def get_matches(self, df, ref, case_def, match_cols):
        """
        Match observations for the crosswalk
        """
        assert ref in df[case_def].unique()
        alts = [x for x in df[case_def].unique() if x != ref]
        assert df[match_cols].notnull().values.all()
        # First match ref to alts
        matches = []
        for alt in alts:
            match = df.query(f"{case_def} == '{ref}'").merge(
                df.query(f"{case_def} == '{alt}'"),
                on=match_cols, validate='one_to_one'
            )
            matches.append(match)

        # Now match alts to alts
        for combo in list(itertools.combinations(alts, 2)):
            match = df.query(f"{case_def} == '{combo[0]}'").merge(
                df.query(f"{case_def} == '{combo[1]}'"),
                on=match_cols, validate='one_to_one'
            )
            matches.append(match)
        df_matched = pd.concat(matches, sort=True)

        # Calculate log ratios and standard errors using delta method
        df_matched['log_ratio'] = np.log(df_matched['cases_y'] / df_matched['cases_x'])
        df_matched['log_ratio_se'] = np.sqrt(
            (df_matched['se_x'] / df_matched['cases_x'])**2
            + (df_matched['se_y'] / df_matched['cases_y'])**2
        )
        return df_matched

    def add_covariates(self, df, covariates=None, add_model_age=True):
        if not covariates:
            covariates = self.covariates
        if add_model_age:
            df = self.add_model_ages(df)
        for cov in covariates:
            if cov not in df:
                if cov in ['europe', 'south_se_asia']:
                    lh = get_current_location_hierarchy(
                        location_set_version_id=self.conf.get_id('location_set_version'),
                        **self.cache_kwargs
                    )
                    region_map = lh.set_index('location_id')['region_name'].to_dict()
                    cov_to_string_match = {
                        'europe': ["Central Europe", "Eastern Europe", "Western Europe"],
                        'south_se_asia': ["South Asia", "Southeast Asia"]
                    }
                    df[cov] = df['location_id'].map(region_map).isin(
                        cov_to_string_match[cov]).astype(int)
                elif cov == 'hosp_continuous':
                    print_log_message(f"Using unknown = {self.unknown}")
                    df['hosp_continuous'] = df['hosp'].map({
                        'community': 0,
                        'unknown': self.unknown,
                        'hospital': 1
                    })
                elif cov == 'microbiology':
                    nid_metadata = pd.read_csv("FILEPATH")
                    if 'source' in df.columns:
                        df.loc[df['source'].isin(nid_metadata['source']), 'microbiology'] = 1
                        df.loc[~df['source'].isin(nid_metadata['source']), 'microbiology'] = 0
                        df.loc[df.source.str.endswith("_lit"), 'microbiology'] = 0
                    else:
                        df['microbiology'] = 0
                else:
                    df = merge_covariate(
                        df, cov, decomp_step=self.conf.get_id("decomp_step"),
                        gbd_round_id=self.conf.get_id("gbd_round")
                    )
                assert df[cov].notnull().all()
        return df

    def create_predictions_template(self, pathogens, ages, years):
        """
        Create a predictions template. Should always have all of the id_cols
        """
        lh = get_current_location_hierarchy(
            location_set_version_id=self.conf.get_id('location_set_version'))
        locs = lh.query("level == 3").location_id.unique().tolist()
        if 'sex_id' in self.covariates:
            sexes = [1, 2]
        else:
            sexes = [3]
        if 'hosp' in self.covariates or 'hosp_continuous' in self.covariates:
            hosp = ['community', 'hospital', 'unknown']
        else:
            hosp = ['all']
        index = pd.MultiIndex.from_product(
            [locs, ages, pathogens, years, sexes, hosp],
            names=[
                'location_id', 'age_group_id', 'pathogen',
                'year_id', 'sex_id', 'hosp'
            ]
        )
        square_df = pd.DataFrame(index=index).reset_index()
        square_df = self.add_covariates(square_df, add_model_age=True)
        square_df = square_df.drop('age_group_id', axis='columns')

        # Add CFR
        merge_cols = [c for c in self.cfr_cols if c != "age_group_id"] + ['agg_age_group_id']
        cfr = self.cfr.rename(columns={'age_group_id': 'agg_age_group_id'})
        square_df = square_df.merge(
            cfr, how='left', validate='many_to_one', on=merge_cols
        )
        print_log_message(
            f"Filling {square_df.loc[square_df.cfr.isnull(), merge_cols].drop_duplicates()}")
        square_df = square_df.merge(
            cfr.loc[cfr.pathogen == 'all'].drop('pathogen', axis='columns')[
                [c for c in merge_cols if c != 'pathogen'] + ['cfr']
            ], how='left', validate='many_to_one',
            on=[c for c in merge_cols if c != 'pathogen']
        )
        square_df['cfr'] = square_df['cfr_x'].fillna(square_df['cfr_y'])
        square_df = square_df.drop(['cfr_x', 'cfr_y'], axis='columns')
        if self.cfr_ready:
            report_if_merge_fail(square_df, 'cfr', merge_cols)
        else:
            square_df['cfr'] = square_df['cfr'].fillna(1)

        # Add metadata
        square_df = add_location_metadata(
            square_df, ['location_name', 'super_region_id'])
        square_df['super_region_name'] = square_df['super_region_id'].map(
            lh.set_index('location_id')['location_name_short'].to_dict()
        )
        return square_df

    def one_hot_encode(self, df):
        # Drop all of the factor covs from the list of covariates
        self.covariates = [cov for cov in self.covariates if cov not in self.factor_covs]
        new_factor_covs = []
        for col in self.factor_covs:
            if df[col].dtype == 'float64':
                df[col] = df[col].astype('int')
            add_cols = pd.get_dummies(
                df[col], prefix=col, prefix_sep='', drop_first=True)
            self.covariates = list(set(self.covariates + add_cols.columns.tolist()))
            new_factor_covs = list(set(new_factor_covs + add_cols.columns.tolist()))
            df = pd.concat([df, add_cols], axis=1)
        self.new_factor_covs = new_factor_covs
        return df

    def set_gaussian_priors(self):
        # Set 0 Gaussian priors on all covariates for every bug
        if self.gprior_sd is not None:
            self.gpriors = {
                pathogen: {cov: [0, self.gprior_sd] for cov in self.covariates}
                for pathogen in self.keep_pathogens + ['other']
            }
        else:
            self.gpriors = {
                pathogen: None for pathogen in self.keep_pathogens + ['other']
            }
        # Exceptions
        if self.infectious_syndrome == 'uti_plus':
            self.gpriors['acinetobacter_baumanii']['haqi'][1] = 0.02
        elif (self.infectious_syndrome == 'cns_infectious') and (
            'for_priors' not in self.model_version
        ):
            if 'non_neonatal' in self.model_version:
                prior_mv = "2021_08_25_for_priors_non_neonatal_retry"
            elif 'neonatal' in self.model_version:
                prior_mv = "2021_08_25_for_priors_neonatal_retry"
            priorbetas = pd.read_csv("FILEPATH")
            covs_for_prior = ['cv_menafrivac']
            for pathogen in self.gpriors.keys():
                for cov in covs_for_prior:
                    self.gpriors[pathogen][cov][0] =\
                        priorbetas.loc[(priorbetas['cov_names'] == cov) & (
                            priorbetas['pathogen'] == pathogen), 'beta'].values[0]

    def run_models(self, df_matched, read_model_cache):
        """Run models with leave one country out (LOCO) cross validation"""
        df_matched = df_matched.reset_index(drop=True)
        # Add iso3 to distinguish countries
        df_matched = add_location_metadata(df_matched, 'iso3', **self.cache_kwargs)
        assert df_matched.iso3.notnull().all()
        iso3s = sorted(df_matched.iso3.unique().tolist())
        print(f"Running out-of-sample validation with {len(iso3s)} holdouts")
        data_splits = [
            df_matched.assign(
               train=lambda d: d['iso3'] != iso3,
               test=lambda d: d['iso3'] == iso3
            ) for iso3 in iso3s
        ] + [df_matched.assign(train=True, test=True)]
        if not read_model_cache:
            dorm_models = [
                DormModel(
                    name=pathogen, covs=self.covariates + ['intercept'],
                    gprior=gprior)
                for pathogen, gprior in self.gpriors.items()
            ]
            with open("FILEPATH", 'wb') as file:
                pickle.dump(dorm_models, file)
            df_matched.to_csv("FILEPATH", index=False)
            worker = "FILEPATH/model_worker.py"
            Path("FILEPATH").mkdir(exist_ok=True)
            print_log_message("Launching workers for modelling...")
            jobs = []
            for holdout in iso3s + ['no_holdout']:
                jobname = f"modelworker_{self.model_version}_"\
                    f"{self.infectious_syndrome}_{holdout}"
                params = [
                    self.model_version, self.infectious_syndrome,
                    holdout, self.ref_pathogen
                ]
                jid = submit_mcod(
                    jobname, language='python', worker=worker,
                    cores=10, memory="10G", params=params,
                    runtime="10:00:00", logging=True, queue="long.q",
                    log_base_dir=self.model_dir
                )
                jobs.append(jid)
            print_log_message("Waiting...")
            wait_for_job_ids(jobs)
            print_log_message("Jobs complete!")

        print_log_message("Reading cached models...")
        models = {}
        for holdout in iso3s + ['no_holdout']:
            out_file = "FILEPATH"
            with open(out_file, 'rb') as file:
                models[holdout] = pickle.load(file)
        data_splits = dict(zip(iso3s + ['no_holdout'], data_splits))
        return data_splits, models

    def create_beta_df(self, model):
        buglist = []
        cov_name = []
        betas = []
        for bug in model.dorms:
            betas.extend(model.beta[model.dorm_model_index[bug]].tolist())
            cov_name.extend(model.dorm_models[bug].covs)
            buglist.extend([bug] * len(model.dorm_models[bug].covs))

        betas = pd.DataFrame({'pathogen': buglist, 'cov_names': cov_name, 'beta': betas})

        betas['beta_sd'] = np.sqrt(np.diagonal(model.beta_vcov))

        betas['beta_pval'] = 2 * t.pdf(
            -1 * np.abs(betas['beta'] / betas['beta_sd']),
            model.data.shape[0] - model.beta.shape[0]
        )

        betas.loc[betas['beta_pval'] <= 0.05, 'signif'] = True
        betas.loc[betas['beta_pval'] > 0.05, 'signif'] = False
        return betas

    def get_residuals(self, model, newdata):
        # Load new data into a new data object
        newdata = newdata.query("test").reset_index().assign(intercept=1)
        newdata_obj = Data.load(
            newdata,
            obs="log_ratio",
            obs_se="log_ratio_se",
            ref_dorm="pathogen_x",
            alt_dorm="pathogen_y",
            dorm_separator="-"
        )
        # Step 1 - construct dorm_model_mats
        # Model matrices (X) for each pathogen
        dorm_model_mats = {
            name: model.dorm_models[name].get_mat(newdata)
            for name in model.dorms
        }
        # X * beta for each pathogen
        dorm_values = model.get_dorm_values(
            beta=model.beta, dorm_model_mats=dorm_model_mats)
        # Weights w for each pathogen
        ref_dorm_weights = model.get_dorm_weights(newdata_obj.ref_dorm)
        alt_dorm_weights = model.get_dorm_weights(newdata_obj.alt_dorm)
        # Calculate residuals
        ref_pred = np.log(np.sum(ref_dorm_weights * dorm_values, axis=1))
        alt_pred = np.log(np.sum(alt_dorm_weights * dorm_values, axis=1))
        newdata['resid'] = newdata_obj.obs.values - (alt_pred - ref_pred)
        return newdata

    def get_prop_residuals(self, df, models):
        df = df.loc[df.iso3.isin(models.keys())]
        df = df.assign(
            pathogen_x=self.ref_pathogen,
            pathogen_y=lambda d: d['pathogen'],
            # Log ratio and log ratio SE are both arbitrary here
            log_ratio=lambda d: np.log(d['cases']),
            log_ratio_se=1
        )
        prop_resids = pd.concat([
            self.get_residuals(
                models[holdout],
                df.assign(
                    test=lambda d: d['iso3'] == holdout if holdout != 'no_holdout' else True
                )
            ).assign(holdout=holdout)
            for holdout in models.keys()
        ])
        prop_resids['prop'] = prop_resids['cases'] / prop_resids.groupby(
            self.id_cols + ['nid', 'holdout']
        ).cases.transform(sum)
        prop_resids['cases_pred'] = np.exp(
            prop_resids['log_ratio'] - prop_resids['resid']
        )
        prop_resids['prop_pred'] = prop_resids['cases_pred'] / prop_resids.groupby(
            self.id_cols + ['nid', 'holdout']
        ).cases_pred.transform(sum)
        return prop_resids

    def generate_point_predictions(self, model, preds):
        print_log_message('Generating point predictions...')
        preds = self.one_hot_encode(preds)
        preds['intercept'] = 1

        # generate case point predictions
        preds = preds.reset_index(drop=True)
        # Model.predict returns columns for each pathogen, sorted
        # alphabetically
        pathogens = sorted(self.keep_pathogens + ['other'])
        preds = pd.concat([
            preds, pd.DataFrame(model.predict(preds), columns=pathogens)
        ], axis=1)
        assert preds.notnull().values.all()
        # Now select the correct column for each row
        preds['prop_cases'] = preds.apply(lambda x: x[x['pathogen']], axis=1)
        preds = preds.drop(pathogens, axis='columns')

        # Calculate prop_deaths
        dem_cols = ['location_id', 'agg_age_group_id', 'sex_id', 'year_id', 'hosp']
        preds['prop_deaths'] = preds['prop_cases'] * preds['cfr']
        preds['prop_deaths'] = preds['prop_deaths'] / preds.groupby(
            dem_cols)['prop_deaths'].transform(sum)

        # Apply any post-processing
        if self.infectious_syndrome == 'respiratory_infectious':
            # Zero out polymicrobial in neonates, community and re-normalize
            zero_out = (preds.agg_age_group_id == 42) & (preds.hosp == 'community') & (
                preds.pathogen == 'polymicrobial'
            )
            preds.loc[zero_out, 'prop_cases'] = 0
            preds.loc[zero_out, 'prop_deaths'] = 0
            preds['prop_cases'] = preds['prop_cases'] / preds.groupby(
                dem_cols)['prop_cases'].transform(sum)
            preds['prop_deaths'] = preds['prop_deaths'] / preds.groupby(
                dem_cols)['prop_deaths'].transform(sum)
        return preds

    def make_plots(self):
        worker = "FILEPATH/make_plots.R"
        params = [
            self.model_dir, self.infectious_syndrome,
            '--cov_cols'
        ] + self.covariates
        if len(self.factor_covs) > 0:
            params += ['--encode_cols'] + self.new_factor_covs

        # Launch a job
        print_log_message("Launching job for making plots...")
        jobname = f"make_plots_{self.model_version}_{self.infectious_syndrome}"
        jid = submit_mcod(
            jobname, language='r', worker=worker,
            cores=1, memory="10G", params=params,
            runtime="00:30:00", logging=True,
            log_base_dir=self.model_dir
        )
        print_log_message("Waiting...")
        wait_for_job_ids([jid])
        print_log_message("Job complete")

    def run(self, read_data_cache=False, read_model_cache=False):
        print_log_message("Getting data...")
        if not read_data_cache:
            formatter = PathogenFormatter(
                model_type='network',
                infectious_syndrome=self.infectious_syndrome,
                keep_pathogens=self.keep_pathogens,
                age_weights_use=self.age_weights_use,
                cache_kwargs=self.cache_kwargs
            )
            df = formatter.format_data()

            print_log_message("Saving age/sex weights for vetting")
            pretty_print(formatter.death_weights).to_csv(
                "FILEPATH",
                index=False)
            pretty_print(formatter.case_weights).to_csv(
                "FILEPATH",
                index=False
            )
            df.to_csv("FILEPATH", index=False)
        else:
            print_log_message("Reading from cache")
            df = pd.read_csv("FILEPATH")
        df['nid'] = df['nid'].apply(lambda x: int(x) if type(x) in [int, float] else x)

        print_log_message("Subsetting ages...")
        df = self.subset_ages(df)

        print_log_message("Getting CFRs...")
        self.get_cfrs()
        print_log_message("Applying CFRs...")
        df = self.apply_cfrs(df)
        print_log_message("Aggregating data...")
        df = self.aggregate_data(df, cols=self.aggregate_cols)
        print_log_message("Calculating standard errors...")
        df = self.calculate_se(df)

        print_log_message("Applying data drops...")
        df = self.drop_data(df)

        print_log_message("Getting matches...")
        if self.ref_pathogen is None:
            self.ref_pathogen = df.groupby(['pathogen'])['cases'].sum().idxmax()
            print_log_message(
                f"No reference pathogen specified, setting to {self.ref_pathogen}"
            )
        df = self.set_composites(df)
        df_matched = self.get_matches(
            df, self.ref_pathogen, 'pathogen',
            self.id_cols + ['nid', 'source']
        )

        print_log_message("Adding covariates...")
        df = self.add_covariates(df)
        df_matched = self.add_covariates(df_matched)

        if self.infectious_syndrome == 'respiratory_infectious':
            # Add Strep pneumo PAFs from vaccine efficacy analysis
            vedata = pd.read_csv("FILEPATH")
            vedata['source'] = vedata['study_type']
            vedata['log_ratio_se'] = vedata['mod_se']
            df_matched = df_matched.append(vedata, sort=False)

        print_log_message("Applying study weights...")
        if self.study_weights is not None:
            df_matched['log_ratio_se'] /= df_matched['source'].map(
                self.study_weights).fillna(1)

        print_log_message("Saving input data...")
        df = pretty_print(df, exclude=['nid'])
        df['iso3'] = df['ihme_loc_id'].str[0:3]
        df.to_csv("FILEPATH", index=False)

        print_log_message("Creating predictions template...")
        preds = self.create_predictions_template(
            self.keep_pathogens + ['other'],
            df.agg_age_group_id.unique().tolist(),
            self.year_ids
        )

        print_log_message("One-hot encoding factor variables...")
        df_matched = self.one_hot_encode(df_matched)
        df = self.one_hot_encode(df)

        print_log_message("Setting Gaussian priors...")
        self.set_gaussian_priors()

        print_log_message("Running model in-sample and with LOCO CV")
        data_splits, models = self.run_models(df_matched, read_model_cache)

        print_log_message("Getting betas and residuals...")
        betas = pd.concat([
            self.create_beta_df(model).assign(holdout=holdout)
            for holdout, model in models.items()
        ])
        # Ratio-space residuals
        resids = pd.concat([
            self.get_residuals(models[holdout], data_splits[holdout]).assign(
                holdout=holdout
            ) for holdout in models.keys()
        ])
        # Proportion-space residuals
        prop_resids = self.get_prop_residuals(df, models)

        print_log_message("Making predictions...")
        preds = self.generate_point_predictions(models['no_holdout'], preds)

        print_log_message("Saving results")
        resids.to_csv("FILEPATH", index=False)
        prop_resids.to_csv("FILEPATH", index=False)
        betas.to_csv("FILEPATH", index=False)
        preds.to_csv("FILEPATH", index=False)

        self.make_plots()
        print_log_message("Done!")


def parse_config(model_version, infectious_syndrome):
    # Read config and parse args
    config_file = pd.read_excel(
        "FILEPATH",
        sheet_name='run'
    )
    config = config_file.query(
        f"model_version == '{model_version}' & "
        f"infectious_syndrome == '{infectious_syndrome}'"
    )
    assert len(config) == 1
    config = config.iloc[0].to_dict()
    # Parse lists and dictionaries
    for param in {
        'covariates', 'agg_age_group_ids', 'keep_pathogens',
        'cfr_use', 'age_weights_use', 'year_ids', 'factor_covs',
        'aggregate_cols', 'study_weights'
    }.intersection(set(config.keys())):
        if not config[param] == 'None':
            config[param] = str(config[param]).split(',')
            if param in ['agg_age_group_ids', 'year_ids']:
                config[param] = [int(x) for x in config[param]]
            if param in ['cfr_use', 'age_weights_use', 'study_weights']:
                param_to_value_type = {
                    'cfr_use': str, 'age_weights_use': str,
                    'study_weights': float
                }
                value_type = param_to_value_type[param]
                config[param] = {
                    x.split(':')[0]: value_type(x.split(':')[1]) for x in config[param]
                }
        else:
            config[param] = None
    for param in ['ref_pathogen', 'age_start', 'age_end', 'unknown', 'gprior_sd']:
        if config[param] == 'None':
            config[param] = None
    for param in ['cfr_ready']:
        config[param] = bool(config[param])

    # Deduce keep_pathogens
    if 'keep_pathogens' not in config.keys() or config['keep_pathogens'] is None:
        syn_path = pd.read_csv("FILEPATH")
        syn_path = syn_path.query(
            f"infectious_syndrome == '{config['infectious_syndrome']}'"
        )
        config['keep_pathogens'] = syn_path['pathogens'].iloc[0].split(', ')
    return config


def save_config(out_dir, config):
    with open("FILEPATH", 'w') as outfile:
        yaml.dump(config, outfile)


if __name__ == '__main__':
    # Input args
    parser = argparse.ArgumentParser(description='Launch pathogen network model')
    parser.add_argument('model_version', type=str)
    parser.add_argument('infectious_syndrome', type=str)
    parser.add_argument('--read_data_cache', action='store_true')
    parser.add_argument('--read_model_cache', action='store_true')
    args = parser.parse_args()
    config = parse_config(args.model_version, args.infectious_syndrome)
    print_log_message(
        f"You submitted the following config: {config}"
    )
    network = PathogenNetwork(**config)
    save_config(network.model_dir, config)
    network.run(read_data_cache=args.read_data_cache, read_model_cache=args.read_model_cache)
