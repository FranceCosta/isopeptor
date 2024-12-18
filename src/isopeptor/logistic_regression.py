#! /usr/env/python
# -*- coding: utf-8 -*-

import joblib
import numpy as np
import os
from pathlib import Path
#import warnings
#import warnings
#warnings.filterwarnings("ignore")

CLF = None

def get_model():
    global CLF
    if CLF is None:
        model_path = Path(__file__).parent / "resources" / "models" / "logistic_regression.pkl"
        if not os.path.isfile(model_path):
            raise ValueError("Model file not found.")
        CLF = joblib.load(model_path)
    return CLF

def predict(rmsd:float, r_asa:float)->float:
    """

        Predicts isopeptide bond probability using rmsd and r_asa

        Args:
            rmsd (float): rmsd with template
            r_asa (float): relative solvent accessible surface

    """
    clf = get_model()
    prob_isopep = clf.predict_proba(np.array([[rmsd, r_asa]]))[:,1]
    return round(prob_isopep[0], 3)