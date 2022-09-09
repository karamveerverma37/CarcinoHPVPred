import pandas as pd
import numpy as np
from pandas import read_csv
import joblib
import sys
np.random.seed(5)
#read dataset
features=sys.argv[1]
model_name=sys.argv[2]
#print(features)
#print(model_name)

df1 = pd.read_csv(features, sep=',')

x_validation=df1[['E6_PDNC_AA', 'E6_PKNC_ACC', 'E6_PKNC_ATC', 'E6_PKNC_TGG', 'E6_PKNC_TGT', 'E6_A3','E6_ENT_NL_A','E6_CDK_TA', 'E6_DDON_T']]
KNN_acc=97.4
LR_acc=94.3
svm_acc=94.5


# data scaling
from sklearn.preprocessing import StandardScaler
scaler_filename = "E6_models/scaler_E6.save"
scaler = joblib.load(scaler_filename) 
x_validation=scaler.transform(x_validation)
#load models

if(model_name == 'knn'):
	filename1 = 'E6_models/KNN_E6_scaled.sav'
	mn1='kNN'
	loaded_model1 = joblib.load(filename1)
	predictions1 = loaded_model1.predict(x_validation)
	prob1 = loaded_model1.predict_proba(x_validation)
	if 'carc' in predictions1:
		prob1=prob1[0:,0]
		predictions1='Carcinogenic'
	else:
		prob1=prob1[0:,1]
		predictions1='Non-carcinogenic'
	print("Our model",mn1, "predicts this HPV as ",predictions1," with probability of","%.2f" % (prob1*100))
elif(model_name == 'svm'):
	filename2= 'E6_models/SVM_E6_scaled.sav'
	mn2='Support vector machine'
	loaded_model2 = joblib.load(filename2)
	predictions2 = loaded_model2.predict(x_validation)
	prob2 = loaded_model2.predict_proba(x_validation)
	if 'carc' in predictions2:
		prob2=prob2[0:,0]
		predictions2='Carcinogenic'
	else:
		prob2=prob2[0:,1]
		predictions2='Non-carcinogenic'
	print("Our model",mn2, "predicts this HPV as ",predictions2," with probability of","%.2f" % (prob2*100))

elif(model_name == 'logistic_regression'):
	filename3 = 'E6_models/LR_E6_scaled.sav'
	mn3='Logistic Regression'
	loaded_model3 = joblib.load(filename1)
	predictions3 = loaded_model1.predict(x_validation)
	prob3 = loaded_model1.predict_proba(x_validation)
	if 'carc' in predictions3:
		prob3=prob3[0:,0]
		predictions3='Carcinogenic'
	else:
		prob3=prob3[0:,1]
		predictions3='Non-carcinogenic'
	print("Our model",mn3, "predicts this HPV as ",predictions3," with probability of","%.2f" % (prob1*100))

elif(model_name == 'all'):
	filename1 = 'E6_models/KNN_E6_scaled.sav'
	filename2 = 'E6_models/LR_E6_scaled.sav'
	filename3= 'E6_models/SVM_E6_scaled.sav'
	mn1='kNN'
	mn2='Logistic Regression'
	mn3='Support vector machine'
	loaded_model1 = joblib.load(filename1)
	loaded_model2 = joblib.load(filename2)
	loaded_model3 = joblib.load(filename3)

	predictions1 = loaded_model1.predict(x_validation)
	prob1 = loaded_model1.predict_proba(x_validation)
	predictions2 = loaded_model2.predict(x_validation)
	prob2 = loaded_model2.predict_proba(x_validation)
	predictions3 = loaded_model3.predict(x_validation)
	prob3 = loaded_model3.predict_proba(x_validation)

	if 'carc' in predictions1:
		prob1=prob1[0:,0]
		predictions1='Carcinogenic'
	else:
		prob1=prob1[0:,1]
		predictions1='Non-carcinogenic'

	if 'carc' in predictions2:
		prob2=prob2[0:,0]
		predictions2='Carcinogenic'
	else:
		prob2=prob2[0:,1]
		predictions2='Non-carcinogenic'
	if 'carc' in predictions3:
		prob3=prob3[0:,0]
		predictions3='Carcinogenic'
	else:
		prob3=prob3[0:,1]
		predictions3='Non-carcinogenic'

	if((predictions1 is predictions2) and (predictions1 is predictions3)):
		#print("Our model",filename1, "predicts this HPV as ",predictions1," with probability of","%.2f" % (prob1*100))
		print("Our models",mn1,", ",mn2," and ",mn3, " predicts this HPV as ",predictions1," with a confidance score of","%.2f" % ((((prob1+prob2+prob3)*100)+KNN_acc+LR_acc+svm_acc)/6))
	else:
		print("Our model",mn1, "predicts this HPV as ",predictions1," with probability of ","%.2f" % (prob1*100))
		print("Our model",mn2, "predicts this HPV as ",predictions2," with probability of ","%.2f" % (prob2*100))
		print("Our model",mn3, "predicts this HPV as ",predictions3," with probability of ","%.2f" % (prob3*100))


#preds=[predictions1,predictions2]
#pred=numpy.array(list)
#numpy.unique(pred)
