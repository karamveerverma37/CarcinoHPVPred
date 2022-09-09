import pandas as pd
import numpy as np
from pandas import read_csv
import joblib
import sys
np.random.seed(5)
#read dataset
features=sys.argv[1]
#print(features)
df1 = pd.read_csv(features, sep=',')

x_validation=df1[['E1_RDK_AG', 'E2_PKNC_ATC', 'E2_PDNC_TC', 'L2_PDNC_AA', 'L2_PKNC_TAG', 'E2_PKNC_TCA', 'E1_A3', 'E1_DDON_A', 'E1_PDNC_GA', 'E2_Entropy', 'E1_PKNC_AAA', 'E6_ENT_NL_A', 'E2_PKNC_AAA', 'E1_PKNC_GGA', 'E1_PKNC_CTT', 'E2_RDK_AG']]


# data normalization
from sklearn.preprocessing import Normalizer
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler

scaler_filename = "scaler.save"
scaler = joblib.load(scaler_filename) 

#transformer = Normalizer()
#x_validation=transformer.fit_transform(x_validation)
x_validation=scaler.transform(x_validation)
#x_validation=MinMaxScaler().fit_transform(x_validation)
filename1 = 'SVM_scaled.sav'
loaded_model = joblib.load(filename1)

predictions = loaded_model.predict(x_validation)
prob = loaded_model.predict_proba(x_validation)

if 'carc' in predictions:
	prob=prob[0:,0]
else:
	prob=prob[0:,1]

print("Our model predicts this HPV as ",predictions," with probability of","%.2f" % (prob*100))
