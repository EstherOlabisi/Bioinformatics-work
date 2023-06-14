## Prediction of cardiovascular disease(CVD) using various classification models


The master script 'classification_and_plots.sh' contains Python commands that execeute the pipeline for CVD prediction. 


----- Running the master script -----
Things to know before you launch the master script that executes the python scripts:

1. You do not need to change or add any arguments or new filename to the master script or python scripts. They should run fine as is.

2. If you would like to explore the program, you can change the python scripts' arguments as preferred. For more information on what each argumument of the python scripts represents, check the next section. 

3. To run the master script, simply type or copy and paste the following on your bash command line: ./classification_and_plots.sh


----- Exploring the python scripts -----

1. classification_models.py: processes input data, produces results for the classification models, and generates input data for the next python script.
    a) To know what each argument represents, type or copy and paste the following in your command line: python3 classification_models.py --help

2. learning_curves.py: takes the output from classification_models.py as its input and produces learning curve plots in png format. The figures will be stored in your current working directory. 
    a) To know what each argument represents, type or copy and paste the following in your command line: python3 learning_curves.py --help

