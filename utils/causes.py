import pandas as pd
import numpy as np
from cod_prep.claude.configurator import Configurator

CONF = Configurator('standard')


def manual_replacement(syndrome_col_series):
    replace_dict = {"diarrhea-e coli": "diarrhea-sp-bacter",
                    "diarrhea-paracolon bacilli": "diarrhea-sp-bacter",
                    "diarrhea-aerobacter aerogenes": "diarrhea-sp-bacter",
                    "diarrhea-aerobacter proteus (mirabilis) (morganii)": "diarrhea-sp-bacter",
                    "diarrhea-staphylococcus": "diarrhea-sp-bacter",
                    "diarrhea-pseudomonas": "diarrhea-sp-bacter",
                    "sepsis/bacteremia-strep-c": "sepsis/bacteremia-sp",
                    "sepsis/bacteremia-strep-g": "sepsis/bacteremia-sp",
                    "sepsis/bacteremia-pneumo": "sepsis/bacteremia-sp",
                    "diarrhea-e coli": "diarrhea-enterotoxigenic e coli",
                    "sepsis-neonat-streptococcus,b": "neonatal sepsis-streptococcus,b",
                    "sepsis-neonat-streptococcus": "neonatal sepsis-streptococcus",
                    "sepsis-neonat-staphylococcus a": "neonatal sepsis-staphylococcus a",
                    "sepsis-neonat-ecoli": "neonatal sepsis-ecoli",
                    "sepsis-neonat-listeriosis": "neonatal spsis-listeriosis"}
    return syndrome_col_series.replace(replace_dict)


def mark_estimated_pathogens(df, cols):
    df["pathogen"] = np.where(df[cols].isnull().all(axis=1), "none", "other")
    df["pathogen"] = np.where(df["severe_pathogen"].notnull(),
                              df["severe_pathogen"], df["pathogen"])
    df["pathogen"] = np.where((df[cols].astype(
        str).applymap(lambda x: "unsp" in x).any(
        1)) & (df.severe_pathogen.isnull()), "unspecified", df["pathogen"])
    df["pathogen"] = np.where((df[cols].astype(str).applymap(
        lambda x: "meningitis-viral" in x).any(
        1)) & (df.severe_pathogen.isnull()),
        "meningitis-viral", df["pathogen"])

    return df


def map_to_estimated_pathogens(df):

    unspecified_dict = add_on_unspecified()
    pathogen_map = pd.read_csv("FILEPATH")
    pathogen_name_dict = pathogen_map.set_index("infectious syndrome level 3")[
        "estimation_pathogen"].to_dict()
    pathogen_number_dict = pathogen_map.set_index("estimation_pathogen")[
        "pathogen_number"].to_dict()
    pathogen_name_dict.update(unspecified_dict)
    pathogen_name_dict.update({"meningitis-viral": "meningitis-viral"})

    name_df = pd.read_excel("FILEPATH")
    name_dict = name_df.set_index("pathogen_number")["pathogen_name"].to_dict()

    pathogen_severity = pd.read_excel("FILEPATH")
    severe_pathogen_dict = pathogen_severity.set_index("pathogen_name")[
        "rank"].to_dict()

    cols = [x for x in list(df) if "infectious_syndrome" in x]
    pathogen_cols = ["named_pathogen_" + str(x) for x in range(1, len(cols) + 1)]
    df[cols] = df[cols].apply(lambda x: manual_replacement(x))

    df[cols] = df[cols].apply(lambda x: x.map(pathogen_name_dict))

    df[pathogen_cols] = df[cols]

    df[cols] = df[cols].apply(lambda x: x.map(pathogen_number_dict))

    df[cols] = df[cols].apply(lambda x: x.map(name_dict))

    df[cols] = df[cols].apply(lambda x: x.map(
        severe_pathogen_dict)).fillna(100).astype(int)
    df["severe_pathogen_int"] = df[cols].apply(lambda x: x.min(), axis=1)

    df.loc[df["severe_pathogen_int"] == 100, "severe_pathogen_int"] = np.NaN

    reverse_dict = dict([value, key]
                        for key, value in severe_pathogen_dict.items())
    df["severe_pathogen"] = df["severe_pathogen_int"].map(reverse_dict)

    df = mark_estimated_pathogens(df, pathogen_cols)
    df.drop(columns=[x for x in list(df)
                     if "infectious_syndrome" in x], inplace=True)

    return df


def add_on_unspecified():

    unspecified = []

    ten_map = pd.read_excel("FILEPATH")
    nine_map = pd.read_excel("FILEPATH")

    unspecified = unspecified + (list(ten_map.loc[ten_map["infectious syndrome level 2"].str.contains(
        "unsp", na=False), "infectious syndrome level 2"].unique()))
    unspecified = unspecified + (list(nine_map.loc[nine_map["infectious syndrome level 2"].str.contains(
        "unsp", na=False), "infectious syndrome level 2"].unique()))
    unspecified = list(set(unspecified))
    unspecified_dict = (dict(zip(unspecified, unspecified)))
    return unspecified_dict
