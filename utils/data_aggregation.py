import pandas as pd
import numpy as np
from db_queries import get_age_metadata
from get_draws.api import get_draws
from amr_prep.utils.pathogen_formatting import PathogenFormatter
from cod_prep.downloaders.ages import add_age_metadata
from cod_prep.utils.misc import report_if_merge_fail, print_log_message
from cod_prep.claude.configurator import Configurator

CONF = Configurator('standard')


def summarize_draws(df, draw_cols):
    df['mean'] = df[draw_cols].mean(axis=1)
    df['upper'] = df[draw_cols].quantile(.975, axis=1)
    df['lower'] = df[draw_cols].quantile(.025, axis=1)

    return df


def create_measure_specific_df(df, measure_name, value_col):
    df = df.copy()
    df["measure"] = measure_name
    df.rename(columns={f"{value_col}": "value"}, inplace=True)
    no_values = list(df.groupby("source", as_index=False)["value"].sum(
    ).query("value==0").source.unique())
    df = df.loc[~(df.source.isin(no_values))]
    return df


def delineate_source_cols(df):
    df = df.copy()
    df.loc[df["source"].str.contains("full_denom"), "source"] = "microbiology_lit_full_denom"
    return df


def assign_measures(df, value_cols):
    df = df.copy()
    df["cases"] = df["admissions"]
    df.loc[(df.admissions > 0) & (df.deaths > 0), "cases"] = df.admissions - df.deaths
    df.loc[df.admissions == df.deaths, "cases"] = 0
    morbidity = create_measure_specific_df(df, "morbidity", value_col="cases")
    mortality = create_measure_specific_df(df, "mortality", value_col="deaths")
    df = pd.concat([morbidity, mortality], sort=True, ignore_index=True)
    df.drop(columns=value_cols, inplace=True)
    df["orig_source"] = df["source"]
    df = delineate_source_cols(df)
    df["source"] = df[["source", "measure"]].apply(
        lambda x: x.source + "_" + x.measure, axis=1)
    return df


def _fix_SOURCE_denom(df, gbd_df, dem_cols, bd_prop_col):
    df = df.merge(gbd_df[dem_cols + ["gbd_prop"]],
                  on=dem_cols)
    df["adjusted_denom"] = df["total_meningitis"] + (df[f"{bd_prop_col}"] * df["gbd_prop"])
    df.drop(["gbd_prop", "total_meningitis"], axis=1, inplace=True)
    df.rename(columns={"adjusted_denom": "total_meningitis"}, inplace=True)
    return df


def fix_bd_denom(df, gbd_df, dem_cols, pathogen):

    df["total_meningitis"] = df.groupby(dem_cols + ["source"])[
        "value"].transform("sum")

    if pathogen == "Streptococcus pneumoniae":
        df = _fix_bd_denom(df.query(f"pathogen=='{pathogen}'"),
                           gbd_df, dem_cols, bd_prop_col="value")
    else:
        strep_df = df.query("pathogen=='Streptococcus pneumoniae'")
        strep_df.rename(columns={"value": "strep_value"}, inplace=True)
        p_df = df.query(f"pathogen=='{pathogen}'")
        df = strep_df.merge(p_df, on=dem_cols + ["measure", "source", "total_meningitis"])
        df = _fix_bd_denom(df.query(f"pathogen_x=='Streptococcus pneumoniae'"),
                           gbd_df, dem_cols, bd_prop_col="strep_value")

    return df


def adjust_bd(df, pathogen, measure):

    dem_cols = ["location_id", "age_group_id", "sex_id", "year_id"]
    draw_cols = ["draw_" + str(x) for x in range(0, 1000)]
    ages = get_age_metadata(gbd_round_id=5, age_group_set_id=12)
    ages = list(ages.age_group_id.unique()) + [28]

    kwargs = {
        "gbd_id_type": "cause_id",
        "source": "como",
        "location_id": 102,
        "gbd_round_id": 5,
        "age_group_id": ages,
        "year_id": range(1990, 2018, 1),
        "sex_id": [1, 2, 3],
        "measure_id": 6}

    morbidity_dict = {
        "measure_id": 6,
        "source": "como"
    }
    mortality_dict = {
        "measure_id": 1,
        "source": "codcorrect"
    }

    if measure == "morbidity":
        kwargs.update(morbidity_dict)
        hib = get_draws(**kwargs, gbd_id=334)
        meningo = get_draws(**kwargs, gbd_id=335)
        strep = get_draws(**kwargs, gbd_id=333)
    elif measure == "mortality":
        kwargs.update(mortality_dict)
        hib = get_draws(**kwargs, gbd_id=334)
        meningo = get_draws(**kwargs, gbd_id=335)
        strep = get_draws(**kwargs, gbd_id=333)

    hm_df = hib.merge(meningo, on=dem_cols)
    for draw in draw_cols:
        hm_df[draw] = hm_df[f"{draw}_x"] + hm_df[f"{draw}_y"]
    hm_df = hm_df[dem_cols + draw_cols]

    gbd_df = hm_df.merge(strep, on=dem_cols)
    for draw in draw_cols:
        gbd_df[draw] = gbd_df[f"{draw}_x"] / gbd_df[f"{draw}_y"]

    gbd_df = summarize_draws(gbd_df, draw_cols)
    gbd_df.rename(columns={"mean": "gbd_prop"}, inplace=True)
    gbd_df = gbd_df[dem_cols + ["gbd_prop"]]

    df = fix_bd_denom(df, gbd_df, dem_cols, pathogen)
    return df


def prep_non_lit_data(pathogen, syndrome, group_cols, value_cols, measure):

    Mapper = PathogenFormatter(infectious_syndrome=syndrome, syndrome_level=1.5)
    df = Mapper.get_computed_dataframe()[group_cols + value_cols]

    df[group_cols] = df[group_cols].fillna("none")
    df[value_cols] = df[value_cols].fillna(0)
    df = df.groupby(group_cols, as_index=False)[value_cols].sum()

    df = assign_measures(df, value_cols)

    print_log_message("adjusting BD denom for meningitis")
    if measure == "morbidity":
        string = ['microbiology_data_partial_denom_morbidity']
        bd_df = adjust_bd(
            df.query(f"location_id==102 & orig_source in @string").reset_index(drop=True),
            pathogen, measure)
    elif measure == "mortality":
        string = ['microbiology_data_partial_denom_mortality']
        bd_df = adjust_bd(
            df.query("location_id==102 & orig_source in @string").reset_index(drop=True),
            pathogen, measure)

    df = df.query("~(location_id==102 & orig_source in @string)").reset_index(drop=True)
    df["total_meningitis"] = df.groupby(["age_group_id", "sex_id", "year_id",
                                         "location_id", "source"])["value"].transform("sum")
    df = df.query(f"pathogen=='{pathogen}'")
    df = pd.concat([df, bd_df], sort=True, ignore_index=True)

    df[f"prop_{pathogen.replace(' ', '_').lower()}"] = df["value"] / df["total_meningitis"]

    return df


def prep_lit_data(group_cols, pathogen_int, syndrome, pathogen):
    df = pd.read_csv("FILEPATH")

    df = add_age_metadata(df, add_cols="age_group_id", merge_col='age_group_name')
    report_if_merge_fail(df, "age_group_name", "age_group_id")

    df["year_id"] = round(((df.year_start + df.year_end) / 2))

    df["source"] = np.where(df.denom_just_certain_bact == 1,
                            "lit_partial_denom_morbidity", "lit_full_denom_morbidity")
    assert len(df.loc[df.source.isnull()]) == 0, "there was an issue in assigning sources"

    df["pathogen"] = pathogen
    df["sex_id"] = df.sex_label.map({'male': 1, 'female': 2, 'unknown': 9})
    df["measure"] = "morbidity"

    df.rename(columns={"prop": f"prop_{pathogen.replace(' ', '_').lower()}",
                       "detailed_sample_size": f"total_{syndrome}", "pathogen_cases": "value"},
              inplace=True)
    df = df[group_cols +
            [f"prop_{pathogen.replace(' ', '_').lower()}", "value", f"total_{syndrome}", "measure"]]

    return df
