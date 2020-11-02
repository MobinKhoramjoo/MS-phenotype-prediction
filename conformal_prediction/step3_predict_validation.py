import numpy as np
from sklearn.svm import SVC
import pandas as pd

import matplotlib.pyplot as plt

# Some plotting code for CP metrics
from pharmbio.cp import plotting

# general utility code
import common_utils
# import * from .common_utils
# from . import common_utils

seed = 19680801

# Load discovery data
X_shuff, y_shuff = common_utils.load_discovery_data()
X_scaled, scaler = common_utils.scale_fit(X_shuff)

# Set seed for further analysis
np.random.seed(seed=seed)

# Instantiate TCP predictor
tcp = common_utils.get_TCP()

# Fit the TCP
tcp.fit(X_scaled,y_shuff)

# Make the predictions - validation data
# ===========================================================================================

# Load validation data
X_valid, y_valid, sample_id = common_utils.load_validation_data()
X_valid_scaled = common_utils.scale(X_valid, scaler)

# Predict
pvals = tcp.predict(X_valid_scaled)

# Save ID, True label (y), p[RRMS], p[PMS]
label0 = common_utils.label_encoder.classes_[0]
label1 = common_utils.label_encoder.classes_[1]
results = pd.DataFrame(data={
    'Sample_ID': sample_id,
    'True_label': y_valid,
    'p['+label0+']':pvals[:,0],
    'p['+label1+']':pvals[:,1]})

# Save the predictions
results.to_csv("results/validation_set_predictions.csv", sep=',')

# Make predictions for transitioners
# ===========================================================================================

# Load transitioners
X_transitioners, sample_id_transitioners = common_utils.load_transitioners()
X_transitioners_scaled = common_utils.scale(X_transitioners, scaler)

pvals_transitioners = tcp.predict(X_transitioners_scaled)

results_trans = pd.DataFrame(data={
    'Sample_ID': sample_id_transitioners,
    'True_label': ['Transition' for i in range(len(X_transitioners))],
    'p['+label0+']':pvals_transitioners[:,0],
    'p['+label1+']':pvals_transitioners[:,1]})

# Save the predictions
results_trans.to_csv("results/transitioners_predictions.csv", sep=',')

# P0 vs P1 plot
markers = [common_utils.marker_dict[label0], common_utils.marker_dict[label1], common_utils.marker_dict['Transitioner']]
colors = common_utils.cm_rrms_pms + ['black']
p1p0_fig = plotting.plot_pvalues(np.hstack((y_valid, np.full(len(pvals_transitioners), 2))), 
    np.vstack((pvals, pvals_transitioners)),
    x_label="p-value ["+ label1+"]", y_label="p-value ["+ label0+"]", cols=[1,0], # Note switched the x / y axes 
    markers=markers, labels=[label0,label1,'Transitioner'], cm=colors)
p1p0_fig.savefig('results/p0p1_plot.pdf')

# Calibration plot
calib_fig = plotting.plot_calibration_curve(y_valid, pvals, labels=[label0, label1], sign_step=0.005, cm=common_utils.cm_rrms_pms)
calib_fig.savefig('results/calibration_plot_discovery_vs_validation.pdf')

# Efficiency plot
label_fig = plotting.plot_label_distribution(y_valid, pvals, figsize=(10,6))
label_fig.savefig('results/label_distribution.pdf')

