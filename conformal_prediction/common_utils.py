
import numpy as np
from sklearn.utils import shuffle
import sklearn.preprocessing as pre
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from nonconformist.cp import TcpClassifier
from nonconformist.nc import NcFactory, ClassifierNc


RNG_SEED = 12013251
# Optimal SVM parameters found in Step1
C_val = 22.899455339537813
gamma_val = 0.003304438107833851

# Path to the data file
data_file_path = "../../data/data_v20201019.csv"

label_encoder = pre.LabelEncoder()
label_encoder.fit(['RRMS','PMS'])

def _get_X_y(df):
    X = df[df.columns[6:]].to_numpy() # 6 columns with other info
    y_txt = df['Condition'].to_numpy()
    
    # Convert to numeric response labels
    y = label_encoder.transform(y_txt)
    return X, y

# Discovery data
# ===========================================================================================

def load_discovery_data():
    '''
    Load the discovery dataset and make a random shuffle of the records
    '''
    
    all_data = pd.read_csv(data_file_path)
    #all_data.columns # check header names
    #all_data.head() # print first rows

    # Subset with training data
    training_data = all_data.loc[(all_data['Cohort']=='MS3') & \
        (all_data['Collection']==1) & \
            (all_data['Condition']!= 'ctrl') & \
                (all_data['Condition'] != 'Transition')]
    
    # Verify selection was OK
    #print(training_data['Condition'].value_counts())
    # 35 SPMS and 35 RRMS

    X, y = _get_X_y(training_data)

    # Shuffle the data
    X_shuff, y_shuff = shuffle(X,y, random_state=RNG_SEED)

    return X_shuff, y_shuff

# Validation data
# ===========================================================================================

def load_validation_data():
    '''
    Load the validation dataset no shuffling needed.
    Returns:
    (X-matrix, y-vector, sample-id vector)
    '''
    
    all_data = pd.read_csv(data_file_path)

    validation_data = all_data.loc[(all_data['Cohort']=='MS1') & \
        (all_data['Collection']==1) & \
            (all_data['Condition']!= 'ctrl') & \
                (all_data['Condition'] != 'Transition')]
    
    # validation_data.head()
    # print(validation_data['Condition'].value_counts())
    # 26 RRMS & 16 SPMS

    X, y = _get_X_y(validation_data)

    sample_id = validation_data['SampleID'].to_numpy()

    return X, y, sample_id

# Transitioning patients
# ===========================================================================================
def load_transitioners():
    '''
    Returns a tuple X, SampleID (matrix, vector)
    '''
    all_data = pd.read_csv(data_file_path)

    transitioners = all_data.loc[(all_data['Condition'] == 'Transition')]
    
    X = transitioners[transitioners.columns[6:]].to_numpy() # 6 columns with other info
    sample_id = transitioners['SampleID'].to_numpy()

    return X, sample_id

def scale_fit(X):
    '''
    Returns the scaled X-matrix and the StandardScaler that can be re-used for other scaling
    '''
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    return X_scaled, scaler

def scale(X, scaler):
    '''
    Scale with a pre-fitted scaler
    Returns:
    The scaled X matrix
    '''
    return scaler.transform(X)

def get_TCP():
    '''
    Convenience method for instantiating a instance of TcpClassifier with the right settings
    Returns:
    a TcpClassifier instance
    '''
    model = SVC(gamma=gamma_val,kernel='rbf',C=C_val, probability=True)	# Create the underlying model
    nc = NcFactory.create_nc(model)	# Create a default nonconformity function
    return TcpClassifier(nc, condition=lambda x: x[1]) # the condition arguments makes it a Mondrian TCP

# Colors for the plots
cm_dict = {'PMS': '#CC3333', 'RRMS': '#3C568E'}
cm_rrms_pms = [cm_dict[label_encoder.inverse_transform([0])[0]], cm_dict[label_encoder.inverse_transform([1])[0]]]
marker_dict = {'PMS' : "^", 'RRMS' : "o", "Transitioner": "*"}
# cm_overall_rrms_pms = ['black'] + cm_rrms_pms
