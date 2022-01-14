import pandas as pd
import pickle
import sys
import getpass
import numpy as np
import random
user = getpass.getuser()
netproppath = 'FILEPATH'
if netproppath not in sys.path:
    sys.path.append('FILEPATH')

BASE_DIR = "FILEPATH"

N_DRAWS = 1000


def make_predictions(infectious_syndrome, model_version, model):
    # Get the point predictions
    df = pd.read_csv("FILEPATH")
    df = df.drop(['prop_cases', 'prop_deaths'], axis='columns')

    random.seed(42)
    print("Generating samples")
    beta_samples = np.random.multivariate_normal(
        mean=model.beta, cov=model.beta_vcov, size=N_DRAWS
    )

    # make draws
    print("Making draws")
    pathogens = sorted(df.pathogen.unique().tolist())
    draws = pd.DataFrame()
    temp = df.drop(['pathogen', 'cfr'], axis='columns')
    temp = temp.drop_duplicates()
    dem_cols = ['location_id', 'agg_age_group_id', 'sex_id', 'year_id', 'hosp']
    assert not temp.duplicated(subset=dem_cols).any()
    temp = temp.reset_index(drop=True)

    for i in range(0, N_DRAWS):
        if i % 50 == 0:
            print(f"Working on draw {i}")
        draw = pd.DataFrame(
            model.predict(temp, beta=beta_samples[i]), columns=pathogens
        )
        draw = pd.concat([temp, draw], axis=1)
        assert draw.notnull().values.all()
        draw['draw'] = f'draw_{i}'
        draws = draws.append(draw, sort=False)
    # Melt pathogens
    draws = pd.melt(
        draws, id_vars=dem_cols + ['draw'],
        value_vars=pathogens, var_name='pathogen', value_name='prop_cases'
    )
    # Merge back on CFR
    cfr = df[dem_cols + ['pathogen', 'cfr']]
    draws = draws.merge(cfr, how='left', on=dem_cols + ['pathogen'], validate='many_to_one')
    # Calculate prop_deaths
    draws['prop_deaths'] = draws['prop_cases'] * draws['cfr']
    draws['prop_deaths'] = draws['prop_deaths'] / draws.groupby(
        dem_cols + ['draw'])['prop_deaths'].transform(sum)

    # Apply any post=processing
    if infectious_syndrome == 'respiratory_infectious':
        # Zero out polymicrobial in neonates, community and re-normalize
        zero_out = (draws.agg_age_group_id == 42) & (draws.hosp == 'community') & (
            draws.pathogen == 'polymicrobial'
        )
        draws.loc[zero_out, 'prop_cases'] = 0
        draws.loc[zero_out, 'prop_deaths'] = 0
        draws['prop_cases'] = draws['prop_cases'] / draws.groupby(
            dem_cols + ['draw'])['prop_cases'].transform(sum)
        draws['prop_deaths'] = draws['prop_deaths'] / draws.groupby(
            dem_cols + ['draw'])['prop_deaths'].transform(sum)

    # Melt cases/deaths
    draws = pd.melt(
        draws, id_vars=dem_cols + ['pathogen', 'draw'],
        value_vars=['prop_cases', 'prop_deaths'],
        var_name='measure_id', value_name='prop'
    )
    draws['measure_id'] = draws['measure_id'].map({'prop_cases': 6, 'prop_deaths': 1})

    # Pivot draws
    draws = pd.pivot_table(
        draws, values='prop', index=dem_cols + ['pathogen', 'measure_id'],
        columns='draw'
    ).reset_index()

    # Some minor renaming
    draws['infectious_syndrome'] = infectious_syndrome
    draws = draws.rename(columns={'agg_age_group_id': 'age_group_id'})
    return draws


def validate_results(df, infectious_syndrome, model_version, diag_dir):
    """Validate against the existing point estimates"""
    # point predictions
    df = df.copy()
    point = pd.read_csv("FILEPATH")

    draws = [f'draw_{i}' for i in range(0, N_DRAWS)]
    df['mean'] = df[draws].mean(axis=1)
    df['upper'] = df[draws].quantile(q=0.975, axis=1)
    df['lower'] = df[draws].quantile(q=0.025, axis=1)

    id_vars = [
        'location_id', 'year_id', 'age_group_id', 'sex_id', 'hosp',
        'pathogen'
    ]
    point = point.rename(columns={'agg_age_group_id': 'age_group_id'})
    point = pd.melt(
        point, id_vars=id_vars,
        value_vars=['prop_cases', 'prop_deaths'],
        var_name='measure_id', value_name='prop'
    )
    point['measure_id'] = point['measure_id'].map({'prop_cases': 6, 'prop_deaths': 1})
    compare = pd.merge(
        df, point, how='outer', on=id_vars + ['measure_id'], validate='one_to_one',
        indicator=True
    )
    assert (compare._merge == 'both').all()

    # rescale the draws
    if infectious_syndrome == 'skin_infectious':
        compare = compare.reset_index()
        compare[draws] = pd.DataFrame(
            compare[draws].to_numpy() * compare[['prop']].to_numpy()
            / compare[['mean']].to_numpy(),
            index=compare.index
        )
        # Renormalize to sum to 1
        group_cols = ['location_id', 'year_id', 'sex_id', 'age_group_id', 'hosp', 'measure_id']
        compare[draws] = pd.DataFrame(
            compare[draws].to_numpy()
            / compare.groupby(group_cols)[draws].transform(sum).to_numpy(),
            index=compare.index
        )
        compare['mean'] = compare[draws].mean(axis=1)
        compare['upper'] = compare[draws].quantile(q=0.975, axis=1)
        compare['lower'] = compare[draws].quantile(q=0.025, axis=1)

    save_comp = compare.drop(draws, axis='columns')
    save_comp.to_csv("FILEPATH")
    return compare


def main(infectious_syndrome, model_version, diag_dir):
    # Read model object
    file_name = "FILEPATH"
    with open(file_name, 'rb') as file:
        model = pickle.load(file)
    df = make_predictions(infectious_syndrome, model_version, model)
    no_dupes = [
        'location_id', 'age_group_id', 'sex_id', 'year_id',
        'hosp', 'measure_id', 'pathogen'
    ]
    assert not df[no_dupes].duplicated().any()

    compare = validate_results(
        df, infectious_syndrome, model_version, diag_dir)
    if infectious_syndrome == 'skin_infectious':
        compare = compare[list(df)]
        df = compare.copy()

    df.to_csv("FILEPATH")
    print('Finished successfully!')


if __name__ == '__main__':
    infectious_syndrome = str(sys.argv[1])
    model_version = str(sys.argv[2])
    diag_dir = str(sys.argv[3])
    main(infectious_syndrome, model_version, diag_dir)
