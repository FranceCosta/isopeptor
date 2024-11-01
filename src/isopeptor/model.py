#! /usr/env/python
# -*- coding: utf-8 -*-

import joblib
import numpy as np
import os
from pathlib import Path
import warnings

CLF = None

def get_model():
    global CLF
    if CLF is None:
        model_path = Path(__file__).parent / "resources" / "model" / "model.pkl"
        if not os.path.isfile(model_path):
            raise ValueError("Model file not found.")
        with warnings.catch_warnings(action="ignore"):
            CLF = joblib.load(model_path)
    return CLF

def _predict(rmsd:float, r_asa:float)->float:
    """

        Predicts isopeptide bond probability using rmsd and r_asa

        Args:
            rmsd: rmsd with template
            r_asa: relative solvent accessible surface

    """
    clf = get_model()
    prob_isopep = clf.predict_proba(np.array([[rmsd, r_asa]]))[:,1]
    return round(prob_isopep[0], 3)