# Prediction of cardiovascular disease(CVD) using various classification models 
* Predicting CVD_ML_project_report.pdf: full report about the machine learning project.
* Test_input_file.csv: small input file required to test the scripts described below. 


## /script directory
This directory contains the master shell script to run the Python scripts developed for this project. The Python scripts are contained in the same directory. 
- Classification_and_plots.sh: the master script for executing the Python pipeline. 
- Classification_models.py: processes input data, produces results of the classification models, and generates input data for the next Python script.
- Learning_curves.py: takes the output from classification_models.py as its input and produces learning curve plots in png format. The figures will be stored in your current working directory. 
- To explore what the arguments represents, type or copy and paste the following in your command line: 
    ```
    python3 classification_models.py --help
    python3 learning_curves.py --help
    ``` 


## How to run the master script
- After downloading the master script, Python scripts, and test csv file, you can run the master script by typing the following on your command line: 
```
./classification_and_plots.sh
```
- If you would like to explore the pipeline further than the default inputs, you can change the arguments based on the --help guide as described above.






