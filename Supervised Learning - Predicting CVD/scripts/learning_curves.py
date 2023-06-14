# ========================================================================= #
# Learning curves for classification 
# reference and inspiration
# source: https://scikit-learn.org/stable/auto_examples/model_selection/plot_learning_curve.html
#
# Code developed and modified by:   Dan Tulpan, dtulpan@uoguelph.ca
# Date:                             March 17, 2021
#
#Further modified by: Esther Olabisi-Adeniyi, April 21, 2023
# How to run:   python3 learning_curves.py -in real_estate.csv  -in1 training_input_file  -in2 training_ouput_File
# ========================================================================= #

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import learning_curve, ShuffleSplit
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import RidgeClassifier
from sklearn import linear_model
from sklearn.neighbors import KNeighborsClassifier
from sklearn import tree
import argparse
import sys
import os



#### Read data from csv file ####
# Assumptions: all data is numeric and has been properly pre-processed
# Return: X and y vectors
def read_data():
    # load the dataset
    X = pd.read_csv(X_train, delimiter=',')
    y = pd.read_csv(y_train, delimiter=',')

    #extract values
    X = X.values
    y = y.values.ravel()

    return X, y


####CREATE A FUNCTION TO SET A CLASSIFIER BASED ON THE OPTION PROVIDED####
def set_classifier(method):
    method_name = {
        "dtc": "Decision Tree Classification",
        "gbc": "Gradient Boosting Classification",
        "knn": "K-Nearest Neighbour Classification",
        "ridgec": "Ridge Classification",
        "rf": "Random Forest Classification"
    }

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


####CREATE A FUNCTION TO SET PLOT PROPERTIES####
def plot_learning_curve(classifier, title, X, y, axes=None, ylim=None, cv=None,
                        n_jobs=None, train_sizes=np.linspace(.1, 1.0, 5)):
    if axes is None:
        _, axes = plt.subplots(1, 2, figsize=(10, 5))

    axes[0].set_title(title)
    if ylim is not None:
        axes[0].set_ylim(*ylim)
    axes[0].set_xlabel("Training examples")
    axes[0].set_ylabel("Score")

    train_sizes, train_scores, test_scores, fit_times, _ = \
        learning_curve(classifier, X, y, cv=cv, n_jobs=n_jobs,
                       train_sizes=train_sizes,
                       return_times=True)
    train_scores_mean = np.mean(train_scores, axis=1)
    train_scores_std = np.std(train_scores, axis=1)
    test_scores_mean = np.mean(test_scores, axis=1)
    test_scores_std = np.std(test_scores, axis=1)
    fit_times_mean = np.mean(fit_times, axis=1)
    fit_times_std = np.std(fit_times, axis=1)

    # Plot learning curve
    axes[0].grid()
    axes[0].fill_between(train_sizes, train_scores_mean - train_scores_std,
                         train_scores_mean + train_scores_std, alpha=0.1,
                         color="r")
    axes[0].fill_between(train_sizes, test_scores_mean - test_scores_std,
                         test_scores_mean + test_scores_std, alpha=0.1,
                         color="g")
    axes[0].plot(train_sizes, train_scores_mean, 'o-', color="r",
                 label="Training score")
    axes[0].plot(train_sizes, test_scores_mean, 'o-', color="g",
                 label="Cross-validation score")
    axes[0].legend(loc="best")

    # Plot n_samples vs fit_times
    axes[1].grid()
    axes[1].plot(train_sizes, fit_times_mean, 'o-')
    axes[1].fill_between(train_sizes, fit_times_mean - fit_times_std,
                         fit_times_mean + fit_times_std, alpha=0.1)
    axes[1].set_xlabel("Training examples")
    axes[1].set_ylabel("fit_times")
    axes[1].set_title("Scalability of the model")

    return plt


# ============================= #
# MAIN PROGRAM #
# ============================= #

#Define command line arguments
parser = argparse.ArgumentParser(description='Classification Model Plots')
parser.add_argument('--in_file1', '-in1', action="store", dest='in_file1', required=True, help='Name of file containing training variables (X).')
parser.add_argument('--in_file2', '-in2', action="store", dest='in_file2', required=True, help='Name of file containing training output (y).')

# handle user errors
try:
    args = parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

# save arguments in separate variables
X_train   = args.in_file1
y_train   = args.in_file2


####READ DATA####
X, y = read_data()


####PLOT AND SAVE FIGURES####
m_list = ['ridgec', 'dtc', 'gbc', 'knn', 'rf']
num_splits = 10
#score curves, each time with 20% data randomly selected as a validation set.
cv = ShuffleSplit(n_splits=num_splits, test_size=0.2, random_state=0) 

#loop through models
for model in m_list:
  classifier, model_name = set_classifier(model)
  title = r"Learning Curves " + model_name
  plot_learning_curve(classifier, title, X, y, axes=None, cv=cv, n_jobs=4)
  #save plots in file
  plt.savefig(model+'.png')
  #plt.show()


#remove intermediate files when finished
os.remove('X_train.csv')
os.remove('y_train.csv')