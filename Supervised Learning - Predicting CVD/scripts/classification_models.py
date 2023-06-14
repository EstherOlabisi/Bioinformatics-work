import numpy as np
import pandas as pd
from numpy.random import seed, randn
from numpy import mean, std, min, max, percentile
from pandas.api.types import CategoricalDtype
import matplotlib.pyplot as plt
from sklearn.feature_selection import f_classif, SelectKBest, RFE
from sklearn.model_selection import train_test_split, ShuffleSplit, RepeatedKFold, cross_val_score, cross_val_predict
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn import linear_model
from sklearn.neighbors import KNeighborsClassifier
from sklearn import model_selection
from sklearn.metrics import balanced_accuracy_score, f1_score, precision_score, recall_score
from sklearn import tree
import argparse
import sys
import os


# ===================== #
# Data pre-processing #
# ===================== #
# This section developed by: Esther Olabisi-Adeniyi
# April 21, 2023


# ============================= #
# Define command line arguments #
# ============================= #
parser = argparse.ArgumentParser(description='Data processing and classification models')
parser.add_argument('--in_file', '-in', action="store", dest='in_file', required=True, help='Name of csv input file.')
parser.add_argument('--kfold', '-k', action="store", dest='kfold', default=5, required=False, help='Number of folds for cross-validation')
parser.add_argument('--num_repeats', '-n', action="store", dest='num_repeats', default=3, required=False, help='Number of repeated CV runs')

#handle user errors
try:
  args = parser.parse_args()
except:
  parser.print_help()
  sys.exit(0)

# define the location of the dataset
#filename='data_missing.csv'
filename    = args.in_file
kf          = int(args.kfold)
num_repeats  = int(args.num_repeats)

####READ FILE####
#encode the file such that non-standard characters are removed
df = pd.read_csv(filename, header=0, sep=";", encoding="utf-8", encoding_errors="ignore")


# ============================= #
# Data pre-processing #
# ============================= #
df = df.iloc[:, [1,7,9,10,11,12]]  #select columns of interest

#confirm non-ascii characters are absent
#str_df = df.to_string()
#if str_df.isascii() == True:
# print("---------\n")
# print ('\n-There are no non-standard characters\n')
#else: 
#print("---------\n")
#print ('\n-There is at least one non-standard character')
#print("---------\n")


####CHECK FOR MISSING VALUES####
#print("Number of missing values per column")
missing_v = pd.isnull(df).sum()
#print(missing_v)
#print("---------\n")


####DROP DUPLICATE ROWS####
#We will drop duplicate rows. The columns are all unique so those are fine as is.  
#print ('-We started with:', df.shape)
df.drop_duplicates(inplace=True)
#print ('-After dropping duplicate rows, we now have:', df.shape, '\n')  
#print("---------\n")


####CHECK FOR CONFOUNDING RECORDS####
#Subset features and check for duplicates among them because this would mean we have records with the same inputs but different outputs
X = df.drop(['cardio'], axis=1)
Y = df['cardio']
dups = X[X.duplicated()]
dups = dups.index.tolist()   
#drop them
df.drop(index=dups, inplace=True)  
#print('-After dropping confounding records, we now have:', df.shape)
#print("---------\n")


####CHECK FOR OUTLIERS IN 'AGE'####
column_to_check = df['age']
# calculate the inter-quartile range
q25, q75 = percentile(column_to_check, 25), percentile(column_to_check, 75)
iqr = q75 - q25
#print('Percentiles: 25th=%.3f, 75th=%.3f, IQR=%.3f' % (q25, q75, iqr))

# calculate the outlier cutoff: k=2
cut_off = iqr * 2
lower, upper = q25 - cut_off, q75 + cut_off

# view outliers 
data_outliers = [x for x in column_to_check if x < lower or x > upper]
#print('Number of identified outliers: %d' % len(data_outliers))
#print('Outliers: ', data_outliers, "\n")


####ADJUST VARIABLE TYPES###
#make binary variables dummies
pd.get_dummies(df, columns= ['smoke', 'alco', 'active', 'cardio'], drop_first = True, prefix_sep='')

#make cholesterol an ordinal categorical variable 
cdtype = CategoricalDtype(
  categories=[1, 2, 3],
  ordered=True)
df['cholesterol'] = df['cholesterol'].astype(cdtype)

#check number of classes
#print(df['cardio'].value_counts())

####SPLIT DATASET####
data = df.values
X = data[:,:-1]  #features
y = data[:,-1]  #target


# split into train-test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.80, random_state=2)


#print("X_train, X_test size: ", X_train.shape, X_test.shape)
#print("y_train, y_test size: ", y_train.shape, y_test.shape)
#print("\nX_train:\n", X_train, "\n")
#print("X_test:\n", X_test, "\n")
#print("y_test:\n", y_test, "\n")
#print("y_train:\n", y_train, "\n")


# ============================= #
# Feature Selection #
# ============================= #

#selectKBest to choose best 3 features
f_sel = SelectKBest(score_func=f_classif, k=3)
fit=f_sel.fit(X_train,y_train)
#print("Select KBest results:", fit.scores_)
#the best features
#print(fit.transform(X_train))
#print("---------\n")

#apply RFE to choose best 3 features
#build logisitic regression model
model = LogisticRegression(solver='lbfgs')


#apply recursive feature elimination based on the LR model
rfe = RFE(model, n_features_to_select = 3)
fit2 = rfe.fit(X_train, y_train)

#format results
#print("RFE results:")
#print("Num Features: %d" % fit2.n_features_)
#print("Selected Features: %s" % fit2.support_)
#print("Feature Ranking: %s" % fit2.ranking_)
mask = fit2.support_
arr = np.empty((1, 1), dtype=bool)
arr.fill(0)
mask = np.append(mask,arr)
selected_features = df.columns[mask]
#print("Selected features: ", selected_features.values)

#training data with top 3 columns only
X_train = X_train[:, 0:3]

#export arrays to create learning curve plots
np.savetxt("X_train.csv", X_train, delimiter=",", fmt='%.f')
np.savetxt("y_train.csv", y_train, delimiter=",", fmt='%.f')


# ===================== #
# Classification models #
# ===================== #
# This section developed by:   Dan Tulpan, dtulpan@uoguelph.ca
# Date:     March 16, 2021
# Modified by: Esther Olabisi-Adeniyi, April 21, 2023


method_name = {
    "ridgec": "Ridge Classification",
    "dtc": "Decision Tree Classification",
    "gbc": "Gradient Boosting Classification",
    "knn": "K-Nearest Neighbour Classification",
    "rf": "Random Forest Classification",
}

scoring = { 'Accuracy':'accuracy',
            'Balanced accuracy':'balanced_accuracy',
            'Precision':'precision_macro',
            'Recall':'recall_macro',
            'F1-score':'f1_macro'
        }

# Set classifier model
# Return: classifier and method name (for printing purposes)
def set_classifier(method):
    if (method == "ridgec"):
        classifier = RidgeClassifier()

    elif (method == "dtc"):
        classifier = tree.DecisionTreeClassifier()

    elif (method == "gbc"):
        classifier = GradientBoostingClassifier(n_estimators=100)

    elif (method == "knn"):
        classifier = KNeighborsClassifier(n_neighbors=5)

    elif (method == "rf"):
        classifier = RandomForestClassifier(n_estimators=100, random_state=0)

    else:
        print("\nError: Invalid method name:" + method + "\n")
        parser.print_help()
        sys.exit(0)
    return classifier, method_name[method]


####EVALUATE CLASSIFIERS USING KFOLD CROSS VALIDATION####
def eval_model(classifier, num_sp, num_rep):
    kfold = RepeatedKFold(n_splits=num_sp, n_repeats=num_rep, random_state=1)

    num_characters = 20
    print("Model".ljust(num_characters),":", method_name[method])
    print("K-folds".ljust(num_characters),":", num_sp)
    print("Num splits".ljust(num_characters),":", num_rep)

    for name,score in scoring.items():
        results = model_selection.cross_val_score(classifier, X_train, y_train, cv=kfold, scoring=score, n_jobs=-1)
        print(name.ljust(num_characters), ": %.3f (%.3f)" % (np.absolute(results.mean()), np.absolute(results.std())))

    #testing evaluation scores
    y_pred = cross_val_predict(classifier, X_test, y_test, cv=num_sp, n_jobs=-1)
    pred_acc = balanced_accuracy_score(y_test, y_pred)
    pred_f1 = f1_score(y_test, y_pred, average='macro')
    pred_prec = precision_score(y_test, y_pred, average='macro')
    pred_rec = recall_score(y_test, y_pred, average='macro')
    print ('[predicted balanced accuracy, f1_score, precision, recall]:', pred_acc, pred_f1, pred_prec, pred_rec)


####CALL CLASSIFIERS####
#define 'k' for num_splits and number of repeats 'n'
kf = 10
num_repeats = 10

method_list = ['ridgec', 'dtc', 'gbc', 'knn', 'rf']
for method in method_list:
  #set classifier
  classifier, name = set_classifier(method)
  #evaluate model
  eval_model(classifier, kf, num_repeats)
  print("---------\n")






