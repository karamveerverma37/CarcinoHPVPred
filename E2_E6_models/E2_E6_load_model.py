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
df1 = pd.read_csv(features, sep=',')

x_validation=df1[['E2_CDK_TA', 'E2_CDK_CT', 'E2_RDK_TA', 'E2_Entropy', 'E2_PKNC_ATC', 'E2_PKNC_TGA', 'E6_PDNC_AA', 'E6_PKNC_ACC', 'E6_PKNC_ATC', 'E6_PKNC_TGG', 'E6_PKNC_TGT', 'E6_A3']]
LR_acc=100
svm_acc=100
KNN_acc=99.8


# data scaling
from sklearn.preprocessing import StandardScaler
scaler_filename = "E2_E6_models/scaler_E2_E6.save"
scaler = joblib.load(scaler_filename) 
x_validation=scaler.transform(x_validation)
#load models


if(model_name == 'logistic_regression'):
	filename1 = 'E2_E6_models/LR_E2_E6_scaled.sav'
	mn1='Logistic Regression'
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
	filename2= 'E2_E6_models/SVM_E2_E6_scaled.sav'
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

elif(model_name == 'knn'):
	filename3= 'E2_E6_models/KNN_E2_E6_scaled.sav'
	mn3='k-nearest neighbors'
	loaded_model3 = joblib.load(filename3)
	predictions3 = loaded_model3.predict(x_validation)
	prob3 = loaded_model3.predict_proba(x_validation)
	if 'carc' in predictions3:
		prob3=prob3[0:,0]
		predictions3='Carcinogenic'
	else:
		prob3=prob3[0:,1]
		predictions3='Non-carcinogenic'
	print("Our model",mn3, "predicts this HPV as ",predictions3," with probability of","%.2f" % (prob3*100))


elif(model_name == 'all'):
	filename1 = 'E2_E6_models/LR_E2_E6_scaled.sav'
	filename2= 'E2_E6_models/SVM_E2_E6_scaled.sav'
	filename3= 'E2_E6_models/KNN_E2_E6_scaled.sav'
	mn1='Logistic Regression'
	mn2='Support vector machine'
	mn3='k-nearest neighbors'
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
		print("Our models",mn1,", ",mn2," and ",mn3, " predicts this HPV as ",predictions1," with a confidance score of","%.2f" % ((((prob1+prob2+prob3)*100)+LR_acc+svm_acc+KNN_acc)/6))
	else:
		print("Our model",mn1, "predicts this HPV as ",predictions1," with probability of","%.2f" % (prob1*100))
		print("Our model",mn2, "predicts this HPV as ",predictions2," with probability of","%.2f" % (prob2*100))
		print("Our model",mn3, "predicts this HPV as ",predictions3," with probability of","%.2f" % (prob3*100))



#preds=[predictions1,predictions2]
#pred=numpy.array(list)
#numpy.unique(pred)

