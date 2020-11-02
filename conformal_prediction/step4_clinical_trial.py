import numpy as np
import pandas as pd
from sklearn.utils import shuffle

# general code for utils
import common_utils

# Some plotting code for CP metrics
from pharmbio.cp import plotting

seed = 19680801
include_both_cohorts = True

all_data = pd.read_csv(common_utils.data_file_path)

def load_data_for_trial(subject_id):
    '''
    Returns the X_train, y_train, X_subject,collection-points (X-train and y-train shuffled)
    '''
    training_data = all_data.loc[
        (all_data['Subject'] != subject_id) & \
            (all_data['Condition']!= 'ctrl') & \
                (all_data['Condition'] != 'Transition')]
    # If only discovery cohort should be used
    if not include_both_cohorts:
        training_data = training_data.loc[ training_data['Cohort'] == 'MS3' ] 
    
    subj_data = all_data.loc[all_data['Subject'] == subject_id]
    # Train data
    X_train, y_train = common_utils._get_X_y(training_data)
    X_shuff, y_shuff = shuffle(X_train,y_train, random_state=seed)
    X_scaled, scaler = common_utils.scale_fit(X_shuff)
    # Test data
    X_subject, _ = common_utils._get_X_y(subj_data)
    X_subject_scaled = common_utils.scale(X_subject, scaler)
    collection_t = subj_data['Collection'].to_numpy()
    return X_scaled, y_shuff, X_subject_scaled, collection_t

trail_data = all_data.loc[ all_data['Trial'] == 1 ] 
subjects = np.unique(trail_data['Subject'])

# Set seed for further analysis
np.random.seed(seed=seed)

# Instantiate TCP predictor
tcp = common_utils.get_TCP()

p_lab0 = "p["+common_utils.label_encoder.classes_[0]+"]"
p_lab1 = "p["+common_utils.label_encoder.classes_[1]+"]"
treated_results = pd.DataFrame(columns=["Subject", "Collection", p_lab0, p_lab1])

for subj in subjects:
    X, y, X_subj, ts = load_data_for_trial(subj)

    tcp.fit(X, y)
    pvals = tcp.predict(X_subj)
    for i in range(X_subj.shape[0]):
        treated_results = treated_results.append({"Subject": subj,
        'Collection': ts[i], 
        p_lab0 : pvals[i,0], 
        p_lab1 : pvals[i,1]
        },
        ignore_index=True)

# Save results
treated_results.to_csv("results/treated_patients_predictions.csv", sep=",", index=False)

## Create a calibration-plot for when merging the two cohorts
# ===========================================================================================
all_data_collection_1 = all_data.loc [ (all_data['Collection'] == 1) & (all_data['Condition'] != 'Transition')]
if not include_both_cohorts:
    all_data_collection_1 = all_data_collection_1.loc [ all_data_collection_1['Cohort'] == 'MS3' ]

all_subjects = np.unique(all_data_collection_1['Subject'])

calib_data = pd.DataFrame(columns=["Subject", "y_label", p_lab0, p_lab1])

for subj in all_subjects:
    train = all_data_collection_1[all_data_collection_1['Subject'] != subj]
    test = all_data_collection_1[all_data_collection_1['Subject'] == subj]
    X, y = common_utils._get_X_y(train)
    X, scaler = common_utils.scale_fit(X)
    
    X_test, y_test = common_utils._get_X_y(test)
    X_test_scaled = common_utils.scale(X_test, scaler)

    tcp.fit(X, y)
    pval = tcp.predict(X_test_scaled)[0]

    calib_data = calib_data.append({'Subject': subj, 'y_label': y_test[0], p_lab0: pval[0], p_lab1 : pval[1]}, ignore_index=True)

# Save for generating calibration curve
calib_data.to_csv("results/MS3+MS1_calib_data.csv")

label0 = common_utils.label_encoder.classes_[0]
label1 = common_utils.label_encoder.classes_[1]

# Calibration plot
calib_fig = plotting.plot_calibration_curve(calib_data['y_label'].to_numpy(), calib_data[[p_lab0, p_lab1]].to_numpy(), labels=[label0, label1], sign_step=0.005, cm=common_utils.cm_rrms_pms)
calib_fig.savefig('results/MS3+MS1_calib_plot.pdf')
