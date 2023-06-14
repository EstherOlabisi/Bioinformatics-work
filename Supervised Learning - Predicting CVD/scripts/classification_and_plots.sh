#!/bin/bash

#Esther Olabisi-Adeniyi
#Date: April 21, 2023
#Bash script for running cardiovascular disease prediction pipeline 
python3 classification_models.py -in small_cvd.csv -k 10 -n 10
python3 learning_curves.py -in1 X_train.csv -in2 y_train.csv
