import pandas as pd
import numpy as np
from pathlib import Path
from db_queries import (
    get_rei_metadata, get_location_metadata,
    get_demographics
)
from get_draws.api import get_draws
from cod_prep.claude.configurator import Configurator
from cod_prep.downloaders import getcache_age_aggregate_to_detail_map

CONF = Configurator()
WORKING_DIR = Path("FILEPATH")


def pull_gbd_deaths_pafs(read_cache=False):
    locs = get_location_metadata(
        location_set_id=CONF.get_id("location_set"),
        gbd_round_id=CONF.get_id("gbd_round")
    )
    country_locs = locs.query("level == 3").location_id.unique().tolist()
    ages = get_demographics("cod", gbd_round_id=CONF.get_id("gbd_round"))['age_group_id']

    if not read_cache:
        df = get_draws(
            ['rei_id'] * 13 + ['cause_id'],
            gbd_id=list(range(173, 186)) + [302],
            source='burdenator',
            location_id=country_locs,
            year_id=2019,
            age_group_id=ages,
            sex_id=[1, 2],
            gbd_round_id=6,
            decomp_step='step5',
            metric_id=2,
            measure_id=1,
            num_workers=10
        )
        df.to_csv("FILEPATH", index=False)
    else:
        df = pd.read_csv("FILEPATH")

    reis = get_rei_metadata(3, gbd_round_id=CONF.get_id("gbd_round"))
    df = df.merge(
        reis[['rei_id', 'rei']], how='left', on='rei_id', validate='many_to_one'
    )
    assert df.rei.notnull().all()
    etiology_to_pathogen = {
      "eti_diarrhea_cholera": "vibrio_cholerae_spp",
      "eti_diarrhea_salmonella": "non_typhoidal_salmonellae",
      "eti_diarrhea_shigellosis": "shigella_spp",
      "eti_diarrhea_epec": "enteropathogenic_escherichia_coli",
      "eti_diarrhea_etec": "enterotoxigenic_escherichia_coli",
      "eti_diarrhea_campylobac": "campylobacter",
      "eti_diarrhea_amoebiasis": "amebiasis",
      "eti_diarrhea_cryptospor": "cryptosporidiosis",
      "eti_diarrhea_rotavirus": "rotavirus",
      "eti_diarrhea_aeromonas": "aeromonas_spp",
      "eti_diarrhea_clostridium": "clostridium_difficile",
      "eti_diarrhea_norovirus": "norovirus",
      "eti_diarrhea_adenovirus": "adenoviruses"
    }
    df['pathogen'] = df['rei'].map(etiology_to_pathogen)
    assert df.pathogen.notnull().all()
    df['infectious_syndrome'] = 'diarrhea'
    keep_cols = [
        'location_id', 'age_group_id', 'sex_id', 'year_id', 'pathogen',
        'infectious_syndrome'
    ] + [c for c in df if 'draw' in c]
    df = df[keep_cols]
    return df


def normalize_deaths(df):
    # normalize any demographic where the sum is greater
    # than 1 to 1
    group_cols = ['infectious_syndrome', 'location_id', 'age_group_id', 'sex_id', 'year_id']
    df = df.set_index(group_cols + ['pathogen'])
    total = df.groupby(level=group_cols).sum()
    df = df.mask((df < np.Inf).mul(total > 1), other=df.div(total))
    assert (df.groupby(level=group_cols).sum() <= 1.0000001).values.all()
    df = df.reset_index(drop=False)
    assert df.notnull().values.all()

    # calculate the "other" remainder for every demographic
    draw_cols = [c for c in df.columns if 'draw' in c]
    other = df.groupby(group_cols, as_index=False)[draw_cols].sum()
    other[draw_cols] = 1 - other[draw_cols]
    assert other[draw_cols].min().min() > -0.0000001
    other[draw_cols] = np.clip(other[draw_cols], 0, None)
    assert (other[draw_cols] >= 0).all().all() and (
        other[draw_cols] <= 1).all().all()
    other['pathogen'] = 'other'
    df = df.append(other, sort=False)
    assert np.allclose(df.groupby(group_cols)[draw_cols].sum(), 1)
    assert (df[draw_cols] >= 0).all().all()
    return df


def get_cfrs():
    input_dir = "FILEPATH"
    df = pd.concat(
        [pd.read_csv(file) for file in input_dir.iterdir()],
        sort=False
    )
    df = df.drop_duplicates()
    return df


def apply_cfrs(df, cfrs):
    df.loc[df.pathogen.str.contains("escherichia_coli"), 'merge_pathogen'] =\
        'escherichia_coli'
    df.loc[df.pathogen.str.contains("virus"), 'merge_pathogen'] = 'virus'
    df['merge_pathogen'] = df['merge_pathogen'].fillna(df['pathogen'])

    age_map = getcache_age_aggregate_to_detail_map(force_rerun=False, block_rerun=True)
    age_map = age_map.loc[age_map.agg_age_group_id.isin(cfrs.age_group_id.unique())]
    df = df.merge(age_map, how='left', on='age_group_id', validate='many_to_one')
    assert df.agg_age_group_id.notnull().all()

    df['merge_sex'] = 3

    cfrs = cfrs.rename(
        columns={
           'pathogen': 'merge_pathogen', 'age_group_id': 'agg_age_group_id',
           'sex_id': 'merge_sex'
        }
    )
    merge_cols = [
        'location_id',
        'year_id',
        'agg_age_group_id',
        'merge_sex',
        'merge_pathogen'
    ]
    df = df.merge(
        cfrs, how='left', on=merge_cols, validate='many_to_one'
    )
    print(f"No cfrs for rows with pathogens {df.loc[df.predict.isnull(), 'pathogen'].unique()}")
    df = df.merge(
        cfrs.query("merge_pathogen == 'other'").drop('merge_pathogen', axis='columns'),
        how='left', on=[c for c in merge_cols if c != 'merge_pathogen'],
        validate='many_to_one'
    )
    df['cfr'] = df['predict_x'].fillna(df['predict_y'])
    assert df.cfr.notnull().all()

    draws = [c for c in df if 'draw' in c]
    df = df.reset_index()
    cases = df.copy()
    cases[draws] = pd.DataFrame(
        cases[draws].to_numpy() / cases[['cfr']].to_numpy(), index=cases.index)
    group_cols = ['infectious_syndrome', 'location_id', 'age_group_id', 'sex_id', 'year_id']
    cases[draws] = pd.DataFrame(
        cases[draws].to_numpy()
        / cases.groupby(group_cols)[draws].transform(sum).to_numpy(),
        index=cases.index
    )

    df = df.assign(measure_id=1).append(cases.assign(measure_id=6))
    df = df[group_cols + ['pathogen', 'measure_id'] + draws]
    print("Saving results")
    df.to_csv("FILEPATH", index=False)
    return df


if __name__ == '__main__':
    print("Pulling death PAFs from GBD")
    df = pull_gbd_deaths_pafs(read_cache=True)
    print("Normalizing deaths PAFs")
    df = normalize_deaths(df)
    print("Applying CFRs")
    cfrs = get_cfrs()
    df = apply_cfrs(df, cfrs)
