#import the libraries needed
import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier as rfc
from sklearn import svm
from sklearn.datasets import make_classification
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.model_selection import train_test_split

#load the data
data=sc.read("Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection.h5ad")
data
pd.options.display.max_columns=None
pd.DataFrame(data.obs)

#Train a:
#RandomForestClassifier

x = data.X
y = data.obs['Status']
genes=data.var_names
#split train and test data
xtrain,xtest,ytrain,ytest = train_test_split(x,y,test_size=0.2)
#train the rfc
random_forest_c =rfc(n_estimators=100)
random_forest_c.fit(xtrain,ytrain)
#make predictons
predictionsx=random_forest_c.predict(xtest)
print("Accuracy:",metrics.accuracy_score(ytest, predictionsx))
print("Matthews correlation coefficient ", metrics.matthews_corrcoef(ytest, predictionsx))
#rank the genes with higher Gini importance value
ranking = np.argsort(random_forest_c.feature_importances_)[::-1]
genes[ranking[:20]]
# plt.scatter(test,predictionsrf)
# plt.xlabel('True values')
# plt.ylabel('Predictions')
# plt.show()


# #Support Vector Machine
# #divide the dataset and train
x= data.X.cluster_ID = 0
xtrainsvm,xtestsvm,ytrainsvm,ytestsvm = train_test_split(x,y,test_size=0.2)
support_vector_machine = svm.SVC()
##### STOPPED: DO NOT RUN THE CODE BELOW
# THE TRAINING OF A SVM ON ALL THE DATASET TAKES WAY TOO MUCH TIME
support_vector_machine.fit(xtrainsvm,ytrainsvm)
#This step took > 40 minutes (and it has not finished)
predictionsvm=support_vector_machine.predict(xtestsvm)
print("Accuracy:",metrics.accuracy_score(ytestsvm, predictionsvm))
print("Matthews correlation coefficient ", metrics.matthews_corrcoef(ytestsvm, predictionsvm))
ranking = np.argsort(support_vector_machine.coef_[0] ** 2)[::-1])
genes[ranking[:20]]
# # plt.scatter(testsvm,predictionsvm)
# plt.xlabel('True values')
# plt.ylabel('Predictions')
# plt.show()
