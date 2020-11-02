
import numpy as np
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV, KFold, LeaveOneOut

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.ticker as mtick
import pandas as pd
import math
import time

# general code for loading the datasets
import common_utils


num_eval_points=30

# Load discovery MS data - shuffling is done within that method
X_shuff, y_shuff = common_utils.load_discovery_data()

X_scaled, _ = common_utils.scale_fit(X_shuff)

# ===========================================================================================
# Grid search parameters 
# Use the Grid search parameters suggested by LIBSVM
C_range = np.logspace(-5,18,base=2,num=num_eval_points)
gamma_range = np.logspace(-25,2,base=2,num=num_eval_points) 
#print("C-range: " + str(C_range))
#print("gamma-range: " + str(gamma_range))

parameters = {
    'C':C_range, 
    'gamma':gamma_range }

tic = time.perf_counter()
grid_search = GridSearchCV(SVC(kernel='rbf',probability=True), parameters, cv=LeaveOneOut(), return_train_score=False)
grid_search.fit(X_scaled, y_shuff)
toc = time.perf_counter()

print(f"Finished gridsearch in {toc - tic:0.4f} seconds")

params = grid_search.best_params_
print("Best parameters: " + str(params))
#Best parameters: {'C': 22.899455339537813, 'gamma': 0.003304438107833851}
print("Best Grid search score: " + str(grid_search.best_score_))
#Best Grid search score: 0.9142857142857143

# Print the results
pd.DataFrame(grid_search.cv_results_).to_csv("results/grid_search_res.csv")


# Utility function to move the midpoint of a colormap to be around
# the values of interest.
# Taken from: https://scikit-learn.org/stable/auto_examples/svm/plot_rbf_parameters.html

class MidpointNormalize(Normalize):

    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

# Viz the scores
scores = grid_search.cv_results_['mean_test_score'].reshape(len(C_range),
                                                     len(gamma_range))
fig = plt.figure(figsize=(10, 10))
plt.subplots_adjust(left=.2, right=0.95, bottom=0.15, top=0.95)
plt.imshow(scores, interpolation='nearest', cmap=plt.cm.hot,
    norm=MidpointNormalize(vmin=0.6, midpoint=0.82, clip=True))
plt.xlabel('gamma')
plt.ylabel('C')
plt.colorbar()

plt.xticks(np.arange(len(gamma_range)), ["{:.2e}".format(v) for v in gamma_range], rotation=45)
plt.yticks(np.arange(len(C_range)), ["{:.2e}".format(v) for v in C_range]) 
plt.title('Validation accuracy')
plt.tight_layout()
plt.show()
