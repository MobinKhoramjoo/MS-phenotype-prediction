
import numpy as np
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV, KFold, LeaveOneOut

import pandas as pd
import math
import time

# general code for loading the datasets
import common_utils

# There is an issue with the linear kernel which doesn't reach convergence
# We cap it at 1000 iterations - which produces warnings
from sklearn.exceptions import ConvergenceWarning
import warnings
warnings.filterwarnings("ignore", category=ConvergenceWarning)


# Load discovery MS data - shuffling is done within that method
X_shuff, y_shuff = common_utils.load_discovery_data()

X_scaled, _ = common_utils.scale_fit(X_shuff)

# ===========================================================================================
# Grid search parameters (using the LIBSVM default grid-search range from grid.py)
C_range = np.logspace(-5,15,base=2,num=6)
gamma_range = np.logspace(-15,3,base=2,num=10)

rbf_params =  {'kernel':['rbf'], 'C':C_range, 'gamma':gamma_range}
linear_params = {'kernel':['linear'], 'C':C_range}
poly_params = {'kernel':['poly'], 'C':C_range, 'degree': [2,3,4]}
sigmoid_params = {'kernel':['sigmoid'], 'C':C_range, 'coef0':[0,1]}
parameters = [rbf_params, linear_params, poly_params, sigmoid_params]

tic = time.perf_counter()
grid_search = GridSearchCV(SVC(max_iter=1000), parameters, cv=LeaveOneOut(), return_train_score=False)
grid_search.fit(X_scaled, y_shuff)
toc = time.perf_counter()

print(f"Finished gridsearch in {toc - tic:0.4f} seconds")

params = grid_search.best_params_
print("Best parameters: " + str(params))
#Best parameters: {'C': 22.899455339537813, 'gamma': 0.003304438107833851}
print("Best Grid search score: " + str(grid_search.best_score_))
#Best Grid search score: 0.9142857142857143

# Print the results
pd.DataFrame(grid_search.cv_results_).to_csv("results/kernel_grid_search_res.csv")
