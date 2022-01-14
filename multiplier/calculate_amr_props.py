import pandas as pd
import numpy as np
import sys
import os
from cod_prep.claude.configurator import Configurator
from amr_prep.utils.amr_io import AmrResult
from cod_prep.utils import (
    print_log_message, get_norway_subnational_mapping,
    report_if_merge_fail)
from pathlib import Path
from cod_prep.downloaders import (
    get_current_location_hierarchy, add_location_metadata,
    add_population, get_current_cause_hierarchy)
import warnings
from get_draws.api import get_draws
from amr_prep.utils.misc import get_prepped_csbg_universe
from cod_prep.utils.misc import report_duplicates


CONF = Configurator()


class AmrPropCalculator(object):

    def __init__(self, burden_type, year_id, pathogen, aggregate_fraction_of_resistance=False):
        self.burden_type = burden_type
        self.year_id = year_id
        self.pathogen = pathogen
        self.aggregate_fraction_of_resistance = aggregate_fraction_of_resistance
        self.conf = Configurator()
        self.cache_kwargs = {'force_rerun': False, 'block_rerun': True, 'cache_results': False}
        self.dem_cols = ['location_id', 'age_group_id', 'sex_id', 'year_id']
        self.draw_cols = ["draw_" + str(x) for x in range(0, 1000)]

    def subset_or_exceptions(self, df, step, abx_class):
        df = df.copy()
        if self.pathogen in [
            'enteropathogenic_escherichia_coli',
            'enterotoxigenic_escherichia_coli'
        ]:
            return df.query(
                f"pathogen == 'escherichia_coli' & abx_class == '{abx_class}'"
            ).assign(pathogen=self.pathogen)
        elif (self.pathogen == 'salmonella_paratyphi') & (abx_class != 'mdr') & (step == 'rr'):
            return df.query(
                f"pathogen == 'salmonella_typhi' & abx_class == '{abx_class}'"
            ).assign(pathogen=self.pathogen)
        elif (
            (self.pathogen in ['salmonella_paratyphi', 'salmonella_typhi']) & (
                abx_class == 'mdr') & (step == 'rr')
        ):
            return df.query(
                f"pathogen == 'salmonella_typhi' & abx_class == 'sulfa'"
            ).assign(pathogen=self.pathogen, abx_class=abx_class)
        elif (
            (self.pathogen == 'shigella_spp') & (abx_class == 'fluoroquinolone') & (
                step == 'rr') & (self.burden_type == 'fatal')
        ):
            return df.query(
                f"pathogen == 'salmonella_typhi' & abx_class == '{abx_class}'"
            ).assign(pathogen=self.pathogen)
        elif (
            (self.pathogen in ['mycobacterium_tuberculosis']) & (
                abx_class == 'isoniazid')
        ):
            return df.query(
                f"pathogen == 'mycobacterium_tuberculosis' \
                    & (abx_class == 'isoniazid_new' | abx_class == 'isoniazid')"
            ).assign(pathogen=self.pathogen, abx_class=abx_class)
        elif (
            (self.pathogen in ['mycobacterium_tuberculosis']) & (
                abx_class == 'rifampicin')
        ):
            return df.query(
                f"pathogen == 'mycobacterium_tuberculosis' \
                    & (abx_class == 'rifampicin_new' | abx_class == 'rifampicin')"
            ).assign(pathogen=self.pathogen, abx_class=abx_class)
        else:
            return df.query(
                f"pathogen == '{self.pathogen}' & abx_class == '{abx_class}'"
            )

    def get_frac_resist_run_id(self, abx_class):
        if not (
            (self.pathogen == 'mycobacterium_tuberculosis') & (
                abx_class in ['mdr', 'xdr'])
        ):
            bug_drug_universe = pd.read_csv("FILEPATH")
            bug_drug = self.subset_or_exceptions(bug_drug_universe, 'pr', abx_class)
            assert len(bug_drug) == 1,\
                f"{self.pathogen}/{abx_class} is not a valid bug/drug"
            return int(bug_drug['resistance_run_id'].unique()[0])

    def get_stgpr_run_params(self, run_id):
        stgpr_run_params = pd.read_csv("FILEPATH")
        assert len(stgpr_run_params) == 1
        return stgpr_run_params.loc[0].to_dict()

    def get_stgpr_fraction_of_resistance(self, abx_class, run_id):
        if (
            (self.pathogen == 'mycobacterium_tuberculosis') & (
                abx_class in ['mdr', 'xdr', 'isoniazid', 'rifampicin'])
        ):
            tb_path = 'DIRECTORY'
            if abx_class == 'mdr':
                df = pd.read_csv(tb_path + "FILEPATH")
            elif abx_class == 'xdr':
                df = pd.read_csv(tb_path + "FILEPATH")
            elif abx_class in ['isoniazid', 'rifampicin']:
                df = self.get_tb_weighted(abx_class, run_id)
        else:
            run_params = self.get_stgpr_run_params(run_id)
            st_gpr_draws = Path("DIRECTORY")
            lh = get_current_location_hierarchy(
                location_set_version_id=self.conf.get_id("location_set_version"),
                **self.cache_kwargs
            ).query("level == 3").assign(bad_loc=0)
            if run_params['gbd_round_id'] == 7:
                nor = get_norway_subnational_mapping().set_index(
                    'location_id_old')['location_id_new'].to_dict()
                lh['location_id'] = lh['location_id'].replace(nor)
            df = pd.DataFrame()
            for location_id in lh.location_id.unique().tolist():
                loc_file = (st_gpr_draws + "FILEPATH")
                if not loc_file.exists():
                    lh.loc[lh.location_id == location_id, 'bad_loc'] = 1
                else:
                    loc_df = pd.read_csv(loc_file)
                    if self.year_id not in loc_df.year_id.unique():
                        assert self.year_id == 2019 and 2018 in loc_df.year_id.unique()
                        loc_df = loc_df.query(f"year_id == 2018").assign(
                            year_id=2019
                        )
                    else:
                        loc_df = loc_df.query(f"year_id == {self.year_id}")
                    df = df.append(loc_df)
            if lh.bad_loc.sum() > 0:
                raise AssertionError(
                    f"The following locations are missing:\n "
                    f"{lh.loc[lh.bad_loc == 1, 'location_name'].unique().tolist()}"
                )
            if run_params['gbd_round_id'] == 7:
                nor = get_norway_subnational_mapping().rename(
                    columns={'location_id_new': 'location_id'}
                )
                df = df.merge(nor, how='left', on=['location_id'])
                df.loc[df.location_id_old.notnull(), 'location_id'] = df['location_id_old']
            df = df[self.dem_cols + [c for c in df if 'draw' in c]]

            assert run_params['draws'] == 1000
            assert run_params['data_transform'] == 'logit'
            assert (df[self.draw_cols] >= 0).values.all()
            assert (df[self.draw_cols] <= 1).values.all()

        if self.aggregate_fraction_of_resistance:
            print_log_message("Aggregating to super region")
            df = add_location_metadata(
                df, 'super_region_id',
                location_set_version_id=self.conf.get_id('location_set_version'),
                **self.cache_kwargs
            )
            report_if_merge_fail(df, 'super_region_id', 'location_id')
            df = add_population(
                df, pop_run_id=self.conf.get_id('pop_run'), **self.cache_kwargs
            )
            report_if_merge_fail(df, 'population', 'location_id')
            group_cols = ['super_region_id', 'age_group_id', 'sex_id', 'year_id']
            df['weight'] = df['population'] / df.groupby(group_cols).population.transform(sum)
            df[self.draw_cols] = pd.DataFrame(
                df[self.draw_cols].to_numpy() * df[['weight']].to_numpy(),
                index=df.index
            )
            df[self.draw_cols] = df.groupby(group_cols)[self.draw_cols].transform(sum)
            assert (df[self.draw_cols] <= 1).values.all()
            df = df.drop(['super_region_id', 'population', 'weight'], axis='columns')
        df = df.assign(pathogen=self.pathogen, abx_class=abx_class)
        return df

    def get_tb_weighted(self, abx_class, run_id):
        dem_cols = ['location_id', 'age_group_id', 'sex_id', 'year_id']
        draw_cols = ["draw_" + str(x) for x in range(0, 1000)]

        run_id_ret = self.get_frac_resist_run_id(abx_class + '_retreated')
        run_id_prop_ret = self.get_frac_resist_run_id('prop_retreated')

        st_gpr_weights_ret_draws = Path("FILEPATH")
        st_gpr_draws_new = Path("FILEPATH")
        st_gpr_draws_ret = Path("FILEPATH")

        lh = get_current_location_hierarchy(
            location_set_version_id=self.conf.get_id('location_set_version'),
            **self.cache_kwargs
        ).query("level == 3")

        nor = get_norway_subnational_mapping().set_index(
            'location_id_old')['location_id_new'].to_dict()
        lh['location_id'] = lh['location_id'].replace(nor)

        df = pd.DataFrame()
        for location_id in lh.location_id.unique().tolist():
            loc_file_new = st_gpr_draws_new + "FILEPATH"
            loc_file_ret = st_gpr_draws_ret + "FILEPATH"
            loc_file_weight_ret = st_gpr_weights_ret_draws + "FILEPATH"

            loc_df_weight_ret = pd.read_csv(loc_file_weight_ret)
            loc_df_new = pd.read_csv(loc_file_new)
            loc_df_ret = pd.read_csv(loc_file_ret)

            loc_df_weight_ret = loc_df_weight_ret.query(f"year_id == 2018").assign(year_id=2019)
            loc_df_weight_ret = loc_df_weight_ret.loc[
                loc_df_weight_ret["year_id"] == self.year_id]
            loc_df_new = loc_df_new.query(f"year_id == 2018").assign(year_id=2019)
            loc_df_new = loc_df_new.loc[loc_df_new["year_id"] == self.year_id]
            loc_df_ret = loc_df_ret.query(f"year_id == 2018").assign(year_id=2019)
            loc_df_ret = loc_df_ret.loc[loc_df_ret["year_id"] == self.year_id]

            assert len(loc_df_weight_ret) == 1
            assert len(loc_df_new) == 1
            assert len(loc_df_ret) == 1

            loc_df_weight_new = 1 - loc_df_weight_ret[draw_cols]

            loc_df_ret[draw_cols] = pd.DataFrame(
                loc_df_ret[draw_cols].to_numpy() * loc_df_weight_ret[draw_cols].to_numpy(),
                index=loc_df_ret.index
            )
            loc_df_new[draw_cols] = pd.DataFrame(
                loc_df_new[draw_cols].to_numpy() * loc_df_weight_new[draw_cols].to_numpy(),
                index=loc_df_new.index
            )

            df_loc = pd.DataFrame()
            df_loc = df_loc.append(loc_df_ret)
            df_loc = df_loc.append(loc_df_new)
            df_loc[draw_cols] = df_loc.groupby(dem_cols)[draw_cols].transform(sum)
            df_loc = df_loc.drop_duplicates().reset_index(drop=True)
            assert len(df_loc) == 1

            df = df.append(df_loc)

        assert (df[self.draw_cols] >= 0).values.all()
        assert (df[self.draw_cols] <= 1).values.all()

        nor = get_norway_subnational_mapping().rename(
            columns={'location_id_new': 'location_id'}
        )
        df = df.merge(nor, how='left', on=['location_id'])
        df.loc[df.location_id_old.notnull(), 'location_id'] = df['location_id_old']
        df = df[self.dem_cols + [c for c in df if 'draw' in c]]

        df = self.rescale_tb_mono(df)
        df = df.assign(pathogen=self.pathogen, abx_class=abx_class)
        return df

    def rescale_tb_mono(self, df):
        lh = get_current_location_hierarchy(
            location_set_version_id=self.conf.get_id("location_set_version"),
            **self.cache_kwargs
        ).query("level == 3").assign(bad_loc=0)
        ch = get_current_cause_hierarchy(
            cause_set_id=CONF.get_id('computation_cause_set'),
            cause_set_version_id=CONF.get_id('computation_cause_set_version'),
            **self.cache_kwargs
        )
        tbs = []
        for tb in ['tb', 'tb_other']:
            tb = ch.loc[ch['acause'] == tb, 'cause_id'].unique().tolist()
            tb_incd = get_draws(
                "cause_id", tb,
                source="como",
                age_group_id=22,
                sex_id=[3],
                location_id=lh['location_id'].unique().tolist(),
                year_id=self.year_id,
                measure_id=6,
                metric_id=3,
                version_id=CONF.get_id("como_version"),
                release_id=CONF.get_id("release"),
                gbd_round_id=CONF.get_id("gbd_round"),
                num_workers=5
            )

            tb_incd = add_population(
                tb_incd, pop_run_id=CONF.get_id('pop_run'), **self.cache_kwargs)
            report_if_merge_fail(tb_incd, 'population', self.dem_cols)

            tb_incd.loc[:, self.draw_cols] = tb_incd.loc[:, self.draw_cols].multiply(
                tb_incd.loc[:, 'population'], axis="index")
            tb_incd = tb_incd[self.dem_cols + ['cause_id'] + self.draw_cols]
            tbs.append(tb_incd)

        tb_parent = tbs[0].drop(columns=['cause_id'])
        tb_other = tbs[1].drop(columns=['cause_id'])

        parent_over_other = tb_parent.set_index(self.dem_cols).divide(
            tb_other.set_index(self.dem_cols))
        df = df.sort_values(by='location_id').set_index(self.dem_cols).mul(parent_over_other)
        df = df.reset_index()
        assert (df[self.draw_cols] >= 0).values.all()
        assert (df[self.draw_cols] <= 1).values.all()
        return df

    def get_resistance_profiles(self):
        lh = get_current_location_hierarchy(
            location_set_version_id=self.conf.get_id("location_set_version"),
            **self.cache_kwargs
        ).query("level == 3").assign(bad_loc=0)
        profile_path = self.conf.get_resource("resistance_profile")
        assert self.year_id in [2018, 2019], "There are no resistance profiles for before 2018"
        loc_dfs = []
        for location_id in lh.location_id.unique().tolist():
            if self.pathogen in [
                'enteropathogenic_escherichia_coli',
                'enterotoxigenic_escherichia_coli'
            ]:
                pathogen_profile = 'escherichia_coli'
            else:
                pathogen_profile = self.pathogen
            loc_file = Path(profile_path.format(
                pathogen=pathogen_profile, location_id=location_id
            ))
            if not loc_file.exists():
                lh.loc[lh.location_id == location_id, 'bad_loc'] = 1
            else:
                loc_dfs.append(pd.read_csv(loc_file))
        if lh.bad_loc.sum() > 0:
            raise AssertionError(
                f"The following locations are missing:\n "
                f"{lh.loc[lh.bad_loc == 1, 'location_name'].unique().tolist()}"
            )
        df = pd.concat(loc_dfs, sort=False)
        assert {'age_group_id', 'sex_id', 'year_id'}.isdisjoint(df)
        df = df.assign(age_group_id=22, sex_id=3, year_id=self.year_id)
        return df

    def get_relative_risk_of_death(self, abx_classes):
        df = pd.read_csv("FILEPATH")
        df = pd.concat([
            self.subset_or_exceptions(df, 'rr', abx_class)
            for abx_class in abx_classes
        ], sort=False)

        df = df.rename(
            lambda x: 'draw_' + str(int(x.replace("V", "")) - 1) if 'V' in x else x,
            axis='columns'
        )
        if self.pathogen == 'mycobacterium_tuberculosis':
            base_dir = "DIRECTORY"
            if 'mdr' in abx_classes:
                df = df.append(pd.read_csv(base_dir + "FILEPATH"))
            if 'xdr' in abx_classes:
                df = df.append(pd.read_csv(base_dir + "FILEPATH"))
            if 'isoniazid' in abx_classes:
                df = df.append(pd.read_csv(base_dir + "FILEPATH"))
            if 'rifampicin' in abx_classes:
                df = df.append(pd.read_csv(base_dir + "FILEPATH"))
        assert len(df) == len(abx_classes)
        assert (df[self.draw_cols] >= 0).values.all()
        if (df[self.draw_cols] < 1).values.any():
            warnings.warn("Some draws show a protective relative risk")
        return df

    def get_relative_risk_of_los(self, abx_classes):
        all_los = pd.read_csv("FILEPATH")
        all_los_prior = pd.read_csv("FILEPATH")
        all_los_prior = all_los_prior.rename(columns={'Unnamed: 0': 'abx_class'})
        dfs = []
        for abx_class in abx_classes:
            if (self.pathogen == 'mycobacterium_tuberculosis') & (
                abx_class in ['mdr', 'xdr']
            ):
                if abx_class == 'mdr':
                    exposed = 946
                    unexposed = 934
                elif abx_class == 'xdr':
                    exposed = 947
                    unexposed = 934
                get_draws_kwargs = {
                   'gbd_id_type': 'cause_id',
                   'gbd_id': [exposed, unexposed],
                   'source': 'como',
                   'measure_id': [3, 6],
                   'location_id': 1,
                   'year_id': self.year_id,
                   'age_group_id': 22,
                   'sex_id': 3,
                   'metric_id': 3,
                   'version_id': self.conf.get_id("como_version"),
                   'gbd_round_id': self.conf.get_id("gbd_round"),
                   'release_id': self.conf.get_id("release")
                }
                print_log_message("Pulling COMO results for TB LOS RR")
                df = get_draws(**get_draws_kwargs)
                yld_exposed = df.query(
                    f"measure_id == 3 & cause_id == {exposed}").reset_index()[self.draw_cols]
                cases_exposed = df.query(
                    f"measure_id == 6 & cause_id == {exposed}").reset_index()[self.draw_cols]
                yld_unexposed = df.query(
                    f"measure_id == 3 & cause_id == {unexposed}").reset_index()[self.draw_cols]
                cases_unexposed = df.query(
                    f"measure_id == 6 & cause_id == {unexposed}").reset_index()[self.draw_cols]
                df = yld_exposed.div(cases_exposed).div(yld_unexposed.div(cases_unexposed))
                df = df.assign(pathogen=self.pathogen, abx_class=abx_class)
            else:
                df = self.subset_or_exceptions(all_los, 'rr', abx_class)
                if len(df) == 0:
                    df = all_los_prior.query(f"abx_class == '{abx_class}'")
                    df = df.assign(pathogen=self.pathogen)
                assert len(df) == 1
                df = df.rename(
                    lambda x: 'draw_' + str(int(x.replace("V", "")) - 1) if 'V' in x else x,
                    axis='columns'
                )
            dfs.append(df)

        df = pd.concat(dfs, sort=False)
        assert len(df) == len(abx_classes)

        assert (df[self.draw_cols] >= 0).values.all()
        if (df[self.draw_cols] < 1).values.any():
            warnings.warn("Some draws show a protective relative risk")
        return df

    def process_inputs(self, resist_cases, rr):
        rr = pd.merge(
            resist_cases[
                ['pathogen', 'abx_set', 'combinatoric_id', 'abx_class']
            ].drop_duplicates(), rr, how='outer',
            on=['pathogen', 'abx_class'], validate='many_to_one',
            indicator=True
        )
        assert (rr._merge == 'both').all()
        rr = rr.drop('_merge', axis='columns')
        rr_prof = rr.groupby(
            ['pathogen', 'abx_set', 'combinatoric_id'], as_index=False
        )[self.draw_cols].max()

        resist_cases = resist_cases.drop('abx_class', axis=1)\
            .drop_duplicates().reset_index(drop=True)
        resist_unique_cols = [
            'pathogen', 'abx_set', 'combinatoric_id',
            'location_id', 'year_id', 'age_group_id', 'sex_id'
        ]
        rr_unique_cols = [
            'pathogen', 'abx_set', 'combinatoric_id'
        ]

        report_duplicates(resist_cases, resist_unique_cols)
        report_duplicates(rr_prof, rr_unique_cols)
        resist_cases = resist_cases.set_index(resist_unique_cols)
        rr_prof = rr_prof.set_index(rr_unique_cols)

        reall_props = rr.copy().reset_index(drop=True)
        beta_lactams = [
            'penicillin', 'aminopenicillin', 'beta_lactamase_inhibitor',
            'third_gen_ceph', 'fourth_gen_ceph',
            'anti_pseudomonal_penicillin', 'carbapenem'
        ]
        beta_lactams_dict = dict(zip(beta_lactams, list(range(1, 8))))
        reall_props['beta_lactams_rank'] = reall_props['abx_class'].map(beta_lactams_dict)
        reall_props['combo_max_bl'] = reall_props.groupby(
            ['combinatoric_id'], as_index=False).beta_lactams_rank.transform('max')
        reall_props = reall_props.loc[
            (reall_props['beta_lactams_rank'].isna()) | (
                reall_props['beta_lactams_rank'] == reall_props['combo_max_bl']),
        ]
        reall_props.drop(columns=['beta_lactams_rank', 'combo_max_bl'], inplace=True)
        reall_props[self.draw_cols] = reall_props[self.draw_cols] - 1
        all_protective = reall_props.groupby(
            ['pathogen', 'abx_set', 'combinatoric_id']
        )[self.draw_cols].transform(lambda x: (x < 0).all())
        mask = (~all_protective & (reall_props[self.draw_cols] < 0)).to_numpy()
        props = reall_props[self.draw_cols].to_numpy()
        props[mask] = 0
        mask2 = all_protective.to_numpy()
        props[mask2] = 0.001
        reall_props[self.draw_cols] = pd.DataFrame(props, index=reall_props.index)
        reall_props[self.draw_cols] = reall_props[self.draw_cols] / reall_props.groupby(
            ['pathogen', 'abx_set', 'combinatoric_id']
        )[self.draw_cols].transform(sum)

        assert (reall_props[self.draw_cols] >= 0).values.all()
        assert (reall_props[self.draw_cols] <= 1).values.all()
        assert (reall_props != np.inf).values.all()
        reall_props = reall_props.set_index(
            ['pathogen', 'abx_set', 'combinatoric_id', 'abx_class']
        )

        self.resist_unique_cols = resist_unique_cols
        self.resist_cases = resist_cases
        self.rr_prof = rr_prof
        self.reall_props = reall_props

    def calculate_props(self, counterfactual, measure_id):
        print_log_message(
            f"Calculating props for {counterfactual}, measure {measure_id}"
        )
        assert counterfactual in ['no_infection', 'susceptible_infection']
        assert measure_id in [1, 3, 4, 6]  # Deaths, YLDs, YLLS, incidence

        if measure_id in [1, 3, 4]:
            def reord(df):
                return df.reorder_levels(
                    [c for c in self.resist_unique_cols if c in df.index.names]
                )

            group_cols = [
                'pathogen', 'abx_set', 'location_id', 'year_id',
                'age_group_id', 'sex_id'
            ]

            if counterfactual == 'no_infection':
                total_resist = reord(
                    self.resist_cases.groupby(level=group_cols)[self.draw_cols].sum()
                )
                total_relative = reord(self.resist_cases.mul(self.rr_prof).groupby(
                    level=group_cols)[self.draw_cols].sum())
                df = reord(
                    reord(self.resist_cases.mul(self.rr_prof)).div(
                        reord(total_relative.add(1 - total_resist))
                    )
                )
            elif counterfactual == 'susceptible_infection':
                excess_risk = reord(self.resist_cases.mul(self.rr_prof - 1))
                total_excess_risk = reord(
                    excess_risk.groupby(level=group_cols)[self.draw_cols].sum())
                df = reord(excess_risk.div(1 + total_excess_risk))
        elif measure_id == 6:
            if counterfactual == 'no_infection':
                df = self.resist_cases.copy()
            elif counterfactual == 'susceptible_infection':
                df = self.resist_cases.copy()
                df[self.draw_cols] = 0

        assert (df[self.draw_cols] <= 1).values.all()
        if counterfactual == 'no_infection':
            assert (df[self.draw_cols] >= 0).values.all()
        elif counterfactual == 'susceptible_infection':
            assert (df[self.draw_cols] >= -1).values.all()
        assert df.notnull().values.all()
        assert (df != np.inf).values.all()

        print_log_message("Reallocating props back to drugs")
        df = df.mul(self.reall_props)
        df = df.groupby(
            level=['pathogen', 'abx_set', 'abx_class'] + self.dem_cols
        )[self.draw_cols].sum()

        df = df.reset_index(drop=False)
        assert df.notnull().values.all()
        df = df[self.dem_cols + ['pathogen', 'abx_set', 'abx_class'] + self.draw_cols]
        df = df.assign(counterfactual=counterfactual, measure_id=measure_id)
        return df

    def reconcile_counterfactuals(self, df):
        susc = df.loc[df.counterfactual == 'susceptible_infection']
        df = df.loc[df.counterfactual != 'susceptible_infection']
        assert (susc.abx_set == 'all').all()
        assert not susc.abx_class.isin(['all_susceptible', 'all_resistant']).any()

        merge_cols = self.dem_cols + ['pathogen', 'abx_class', 'measure_id']
        no_inf = pd.merge(
            df[merge_cols + self.draw_cols],
            susc[merge_cols], how='right', on=merge_cols,
            validate='one_to_one'
        )
        susc = susc.drop(['abx_set', 'counterfactual'], axis='columns')
        susc = susc.set_index(merge_cols)
        no_inf = no_inf.set_index(merge_cols)
        susc = susc.mask(no_inf.sub(susc) < 0, other=no_inf)
        susc = susc.reset_index().assign(
            abx_set='all', counterfactual='susceptible_infection'
        )[df.columns.tolist()]
        df = df.append(susc)
        assert df.notnull().values.all()
        return df

    def run_calculator(self):
        print_log_message("Loading inputs")
        print_log_message("Getting fraction of resistance in cases")
        csbg = get_prepped_csbg_universe(self.burden_type, add_parents=False)
        abxs = csbg.loc[
            csbg.abx_class.notnull() & (csbg.pathogen == self.pathogen), 'abx_class'
        ].unique().tolist()

        resist_cases = pd.concat([
            self.get_stgpr_fraction_of_resistance(
                abx_class, self.get_frac_resist_run_id(abx_class)
            ) for abx_class in abxs
        ], sort=False).assign(
            abx_set=lambda d: d['abx_class'] + '_only',
            combinatoric_id=lambda d: d['pathogen'] + '-' + d['abx_class'],
        )
        if self.pathogen == 'mycobacterium_tuberculosis':
            resist_cases.loc[
                resist_cases.abx_class.isin(['isoniazid', 'rifampicin']), 'abx_set'
            ] = 'all'
            assert (
                resist_cases.groupby(self.dem_cols + ['pathogen', 'abx_set'])[
                    self.draw_cols].sum() <= 1).values.all()

        if self.pathogen != 'mycobacterium_tuberculosis':
            if len(abxs) > 1:
                print_log_message("Getting fractions for resistance profiles")
                resist_prof = self.get_resistance_profiles()
                assert set(resist_prof.abx_class) == set(abxs)
                assert len(set(resist_prof.combinatoric_id)) == 2 ** len(abxs)
                assert np.allclose(
                    resist_prof.groupby(
                        ['pathogen', 'location_id', 'age_group_id', 'sex_id',
                         'year_id', 'abx_class']
                    )[self.draw_cols].sum(), 1)
                resist_prof = resist_prof.query("resistant == 1").drop(
                    'resistant', axis='columns'
                )
                if self.pathogen in [
                    'enteropathogenic_escherichia_coli',
                    'enterotoxigenic_escherichia_coli'
                ]:
                    resist_prof['pathogen'] = self.pathogen
                    resist_prof['combinatoric_id'] = resist_prof['combinatoric_id'].str.replace(
                        'escherichia_coli', self.pathogen
                    )
                resist_cases = resist_cases.append(resist_prof.assign(abx_set='all'))
            else:
                resist_cases = resist_cases.append(resist_cases.assign(abx_set='all'))

        if self.burden_type == 'fatal':
            print_log_message("Getting relative risk of death")
            rr = self.get_relative_risk_of_death(abxs)
            measures = [1, 4]
        elif self.burden_type == 'nonfatal':
            print_log_message("Getting relative risk of LOS")
            rr = self.get_relative_risk_of_los(abxs)
            measures = [3, 6]

        print_log_message("Prepping inputs")
        self.process_inputs(resist_cases, rr)

        print_log_message("Calculating AMR props for 2 counterfactuals")
        df = pd.concat([
            self.calculate_props(counterfactual, measure_id)
            for counterfactual in ['no_infection', 'susceptible_infection']
            for measure_id in measures
        ], sort=False)
        non_draw_cols = list(set(df.columns.values) - set(self.draw_cols))
        non_draw_cols.remove('abx_class')
        if (self.pathogen == 'mycobacterium_tuberculosis'):
            df['abx_set'] = 'all'
            mono_with = df.loc[
                (df['counterfactual'] == 'no_infection') & (
                    df['abx_class'].isin(['rifampicin', 'isoniazid'])), ]
            mono_with = mono_with.groupby(non_draw_cols, as_index=False)[self.draw_cols].sum()
            mono_with.set_index(non_draw_cols, inplace=True)
            tb_susc = (1 - mono_with)
            tb_susc.reset_index(inplace=True)
            tb_susc['abx_class'] = 'all_susceptible'
            tb_susc['abx_set'] = 'all'
            df = pd.concat([df, tb_susc])
        else:
            no_ddp_with_resis = df.loc[
                (df['counterfactual'] == 'no_infection') & (
                    df['abx_set'] != 'all'), ]
            ddp_with_resis = df.loc[
                (df['counterfactual'] == 'no_infection') & (
                    df['abx_set'] == 'all'), ]
            df = df.loc[(df['counterfactual'] != 'no_infection') & (df['abx_set'] == 'all'), ]

            ddp_with_resis = ddp_with_resis.groupby(
                non_draw_cols, as_index=False)[self.draw_cols].sum()
            ddp_with_resis.set_index(non_draw_cols, inplace=True)
            ddp_susc = (1 - ddp_with_resis)
            ddp_susc.reset_index(inplace=True)
            ddp_with_resis.reset_index(inplace=True)
            ddp_susc['abx_class'] = 'all_susceptible'
            ddp_with_resis['abx_class'] = 'all_resistant'
            df = pd.concat([df, no_ddp_with_resis, ddp_with_resis, ddp_susc])

        df = self.reconcile_counterfactuals(df)

        print_log_message("Saving results")
        AmrResult(
            process='calculate_amr_props',
            burden_type=self.burden_type,
            year_id=self.year_id,
            pathogen=self.pathogen,
        ).write_results(df)


if __name__ == '__main__':
    burden_type = str(sys.argv[1])
    task_id = os.environ.get('SGE_TASK_ID')
    if task_id:
        print(f'Running as array job, task_id: {task_id}')
        task_row = pd.read_csv("FILEPATH").iloc[int(task_id) - 1]
        year_id = int(task_row['year_id'])
        pathogen = str(task_row['pathogen'])
    else:
        year_id = int(sys.argv[2])
        pathogen = str(sys.argv[3])
    calc = AmrPropCalculator(
        burden_type=burden_type,
        year_id=year_id,
        pathogen=pathogen
    )
    calc.run_calculator()
