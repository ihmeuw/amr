from builtins import range
import pandas as pd
from mcod_prep.datasets.mcod_formatting import MultipleCauseFormatter
from cod_prep.utils import report_if_merge_fail, print_log_message
from cod_prep.downloaders.ages import get_ages

DATA_COL = ['deaths']

SOURCE = "SOURCE"

YEARS = list(range(1998, 2018))


def get_file_path(year):
    if year in range(1998, 2012):
        date = "Y2019M02D15"
    elif year in range(2012, 2016):
        date = "Y2019M02D13"
    elif year in range(2016, 2018):
        date = "STATA_V_14_Y2019M02D15"

    file_path = "FILEPATH"

    return file_path


def get_columns(year):
    columns = ["gru_ed1", "sexo", "sit_defun", "c_dir1", "c_ant1", "c_ant2", "c_ant3",
               "c_pat1", "c_bas1"]

    if year in range(2008, 2018):
        columns += ['c_dir12', 'c_ant12', 'c_ant22', 'c_ant32', 'c_pat2']

    if year in list(range(2002, 2004)) + list(range(2008, 2014)) + list(range(2015, 2018)):
        columns = [col.upper() if col != "sexo" else col for col in columns]

    return columns


def read_data(year, test=False):
    file_path = get_file_path(year)
    columns = get_columns(year)

    if test:
        print_log_message("THIS IS A TEST - READING ONLY 1000 ROWS FOR {}".format(year))
        reader = pd.read_stata(
            file_path, iterator=True, convert_categoricals=False, columns=columns)
        df = reader.get_chunk(1000)
    else:
        print_log_message("Reading data for {}".format(year))
        df = pd.read_stata(file_path, convert_categoricals=False, columns=columns)

    df.columns = df.columns.str.lower()

    return df


def get_age_map(year):

    if year in range(1998, 2008):
        age_map = {
            "01": "Early Neonatal",    # Under 1 day
            "02": "Early Neonatal",    # 1-6 days
            "03": "Late Neonatal",     # 7-27 days
            "04": "Post Neonatal",     # 28-29 days
            "05": "Post Neonatal",     # 1-5 months
            "06": "Post Neonatal",     # 6-11 months
            "07": "1 to 4",            # 1 year
            "08": "1 to 4",            # 2-4 years
            "09": "5 to 9",            # 5-9 years
            "10": "10 to 14",          # 10-14 years
            "11": "15 to 19",          # 15-19 years
            "12": "20 to 24",          # 20-24 years
            "13": "25 to 29",          # 25-29 years
            "14": "30 to 34",          # 30-34 years
            "15": "35 to 39",          # Listed as 35-99 years
            "16": "40 to 44",          # 40-44 years
            "17": "45 to 49",          # 45-49 years
            "18": "50 to 54",          # 50-54 years
            "19": "55 to 59",          # 55-59 years
            "20": "60 to 64",          # 60-64 years
            "21": "65 to 69",          # 65-69 years
            "22": "70 to 74",          # 70-74 years
            "23": "75 to 79",          # 75-79 years
            "24": "80 to 84",          # 80-84 years
            "25": "85 plus",           # 85 and more
            "26": "Unknown"            # Age unknown
        }
    elif year in range(2008, 2018):
        age_map = {
            "00": "Early Neonatal",    # Less than 1 hour
            "01": "Early Neonatal",    # Less than 1 day
            "02": "Early Neonatal",    # 1-6 days
            "03": "Late Neonatal",     # 7-27 days
            "04": "Post Neonatal",     # 28-29 days
            "05": "Post Neonatal",     # 1-5 months
            "06": "Post Neonatal",     # 6-11 months
            "07": "1 to 4",            # 1 year
            "08": "1 to 4",            # 2-4 years
            "09": "5 to 9",            # 5-9 years
            "10": "10 to 14",          # 10-14 years
            "11": "15 to 19",          # 15-19 years
            "12": "20 to 24",          # 20-24 years
            "13": "25 to 29",          # 25-29 years
            "14": "30 to 34",          # 30-34 years
            "15": "35 to 39",          # 35-39 years
            "16": "40 to 44",          # 40-44 years
            "17": "45 to 49",          # 45-49 years
            "18": "50 to 54",          # 50-54 years
            "19": "55 to 59",          # 55-59 years
            "20": "60 to 64",          # 60-64 years
            "21": "65 to 69",          # 65-69 years
            "22": "70 to 74",          # 70-74 years
            "23": "75 to 79",          # 75-79 years
            "24": "80 to 84",          # 80-84 years
            "25": "85 to 89",          # 85-89 years
            "26": "90 to 94",          # 90-94 years
            "27": "95 plus",           # 95-99 years
            "28": "95 plus",           # 100 years and more
            "29": "Unknown"            # Age unknown
        }

    return age_map


def create_age_group_id(df, year):
    age_map = get_age_map(year)

    gbd_ages = get_ages()
    gbd_ages = gbd_ages[['age_group_name', 'age_group_id']]
    assert ~gbd_ages.duplicated().any()

    df['age_group_name'] = df.gru_ed1.map(age_map)
    assert df.age_group_name.notnull().all()

    df = df.merge(gbd_ages, how='left', on='age_group_name')
    report_if_merge_fail(df, 'age_group_id', 'age_group_name')
    return df


def create_sex_id(df):

    sex_map = {
        1: 1,   # Male
        2: 2,   # Female
        3: 9,   # Unknown
    }

    df['sex_id'] = df.sexo.astype(int).map(sex_map)
    assert df.sex_id.notnull().all()
    return df


def create_nid(df, year):

    year_to_nid = {
        1998: 397403,
        1999: 397405,
        2000: 397407,
        2001: 397409,
        2002: 397411,
        2003: 397413,
        2004: 397415,
        2005: 397417,
        2006: 397419,
        2007: 397421,
        2008: 65267,
        2009: 65199,
        2010: 57982,
        2011: 265177,
        2012: 265178,
        2013: 265179,
        2014: 265181,
        2015: 265219,
        2016: 265220,
        2017: 396730
    }

    df['nid'] = year_to_nid[year]

    return df


def mark_hospital_deaths(df):
    df.sit_defun = df.sit_defun.astype(int)
    df['hospital_death'] = 0
    df.loc[df.sit_defun == 1, 'hospital_death'] = 1
    return df


def validate_cause_column(df, column):
    assert df[column].notnull().all()

    assert (df[column] != "").all()

    assert df[column].str.len().isin([3, 4]).all()

    assert df[column].str.match("(^[A-Z][0-9]{2,4}$)|(^0000$)").all()


def clean_cause_column(df, column):
    not_like_ICD = ~df[column].str.match("(^[A-Z][0-9]{2,4}$)|(^0000$)")
    df.loc[not_like_ICD, column] = df[column].apply(lambda x: x[:-1])

    not_like_ICD = ~df[column].str.match("(^[A-Z][0-9]{2,4}$)|(^0000$)")
    weird_stuff = df.loc[not_like_ICD, column].tolist()
    indices_changed = df.loc[not_like_ICD].index.tolist()
    df.loc[not_like_ICD, column] = '0000'

    return df, weird_stuff, indices_changed


def create_underlying_cause(df, verbose=False):
    df['cause'] = df.c_bas1.str.rstrip("X").str.upper()

    (df, ws, ic) = clean_cause_column(df, 'cause')

    if verbose:
        print_log_message(
            "Found the following weird causes in cause column: {}".format(set(ws))
        )
        print_log_message("Changed {} deaths".format(len(set(ic))))

    validate_cause_column(df, 'cause')

    missing_ucod = (df['cause'] == '0000')
    assert len(df[missing_ucod]) < 5, \
        "There are {} rows missing an underlying cause of death!".format(len(df[missing_ucod]))
    df = df[~missing_ucod]

    return df


def create_multiple_causes(df, year, verbose=False):
    if year in range(1998, 2008):
        cols = ['c_dir1', 'c_ant1', 'c_ant2', 'c_ant3', 'c_pat1']
    elif year in range(2008, 2018):
        cols = ["c_dir1", "c_dir12", "c_ant1", "c_ant12", "c_ant2", "c_ant22",
                "c_ant3", "c_ant32", "c_pat1", "c_pat2"]

    unknown_causes = []
    indices_changed = []
    for col in cols:
        suffix = col.split("_")[1]
        mc_col = "multiple_cause_" + suffix
        df[mc_col] = df[col].str.rstrip("X").str.upper()
        df.loc[df[mc_col] == "", mc_col] = '0000'

        (df, uc, ic) = clean_cause_column(df, mc_col)
        unknown_causes += uc
        indices_changed += ic
        validate_cause_column(df, mc_col)

    if verbose:
        print_log_message(
            "Found the following unknown causes in the multiple_cause columns: {}"\
            .format(set(unknown_causes))
        )
        print_log_message("Changed {} deaths".format(len(set(indices_changed))))

    return df


def handle_part_II(df, drop_p2=None):
    assert drop_p2 is not None, "Need to know whether to keep part II"

    p2_cols = [col for col in df.columns if col in ['multiple_cause_pat1', 'multiple_cause_pat2']]

    if drop_p2:
        df.drop(p2_cols, axis='columns', inplace=True)
    else:
        mc_cols = [col for col in df.columns if 'multiple_cause' in col]
        for mc_col in mc_cols:
            suffix = mc_col.split("_")[-1]
            indicator_col = "pII_" + suffix
            if mc_col in p2_cols:
                df[indicator_col] = 1
            else:
                df[indicator_col] = 0
    return df


def clean_SOURCE(year, drop_p2=None, test=False):
    assert year in YEARS, "No data available for {}".format(year)

    df = read_data(year, test=test)

    print_log_message("Cleaning data for {}".format(year))

    df = create_age_group_id(df, year)
    df = create_sex_id(df)
    df['year_id'] = year
    df['location_id'] = 125
    df['code_system_id'] = 1
    df = create_nid(df, year)
    df = mark_hospital_deaths(df)

    df = create_underlying_cause(df)
    df = create_multiple_causes(df, year)

    df = handle_part_II(df, drop_p2=drop_p2)

    df['deaths'] = 1

    mcodformatter = MultipleCauseFormatter(df, source=SOURCE, data_type_id=9, drop_p2=drop_p2)
    df = mcodformatter.get_formatted_dataframe()

    return df


if __name__ == '__main__':
    test = False
    write_metadata = True
    drop_p2 = True
    for year in YEARS:
        df = clean_SOURCE(year, drop_p2=drop_p2, test=test)
        if write_metadata:
            mcodformatter = MultipleCauseFormatter(
                df, source=SOURCE, data_type_id=9, drop_p2=drop_p2)
            mcodformatter.write_metadata(df)
