import sys
import getpass
import pandas as pd
import argparse
import pickle
from cod_prep.utils import print_log_message
from modelling.pathogen.run_model import PathogenNetwork
sys.path.append("FILEPATH")
from netprop.data import Data
from netprop.dorm_model import DormModel
from netprop.model import Model


def fit_model(df, dorm_models, ref_pathogen):
    """Fit the model"""
    data = Data.load(
        df.query("train").reset_index(),
        obs="log_ratio",
        obs_se="log_ratio_se",
        ref_dorm="pathogen_x",
        alt_dorm="pathogen_y",
        dorm_separator="-"
    )
    model = Model(data, dorm_models, gold_dorm=ref_pathogen)
    print_log_message('Fitting model...')
    model.fit_model()
    print_log_message('Model finished!')
    return model


def main(model_version, infectious_syndrome, holdout, ref_pathogen):
    print_log_message(f"Running for holdout {holdout}")
    model_dir = "FILEPATH"
    df = pd.read_csv("FILEPATH")
    with open("FILEPATH", 'rb') as file:
        dorm_models = pickle.load(file)
    if holdout != 'no_holdout':
        df = df.assign(
           train=lambda d: d['iso3'] != holdout,
           test=lambda d: d['iso3'] == holdout
        )
    else:
        df = df.assign(train=True, test=True)
    assert len(df.query("test")) > 0, "Your holdout resulted in no test data"
    model = fit_model(df, dorm_models, ref_pathogen)
    out_file = "FILEPATH"
    with open(out_file, 'wb') as file:
        pickle.dump(model, file)
    print('Model saved, exiting...')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Model worker for pathogens')
    parser.add_argument('model_version', type=str)
    parser.add_argument('infectious_syndrome', type=str)
    parser.add_argument('holdout', type=str)
    parser.add_argument('ref_pathogen', type=str)
    args = parser.parse_args()
    main(
        args.model_version, args.infectious_syndrome, args.holdout, args.ref_pathogen
    )
