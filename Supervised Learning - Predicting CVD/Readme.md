# Prediction of cardiovascular disease(CVD) using various classification models

## /script directory
This directory contains the master shell script to run the python scripts developed for this project. The python scripts are in the same directory. 
- Classification_models.py: processes input data, produces results for the classification models, and generates input data for the next python script.
- Learning_curves.py: takes the output from classification_models.py as its input and produces learning curve plots in png format. The figures will be stored in your current working directory. 
- To explore what the arguments represents, type or copy and paste the following in your command line: 
    ```
    python3 classification_models.py --help
    python3 learning_curves.py --help
    ``` 

## How to run the master script
- After downloading the master script and python scripts, you can run the master script by typing the following on your command line: 
```
./classification_and_plots.sh
```
- If you would like to explore the pipeline further than the default inputs, you can change the arguments using the --help option as described above. For more information on what each argumument of the python scripts represents, check the next section. 






