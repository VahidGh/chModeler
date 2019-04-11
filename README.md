# chModeler

[![Build Status](https://travis-ci.com/VahidGh/chModeler.svg?token=AUxS1pcZ8QLGaaPqajio&branch=master)](https://travis-ci.com/VahidGh/chModeler)

A new approach to the old Hodgkin-Huxley model for:
Making all ion channel models comparable

This package builds models for ion channels from patch clamp data using Binomial Distribution.

The default dataset contains ~200 voltage-gated ion channel data from different species with similar experiment conditions where data have been gathered from literature and digitized for modeling and comparision studies.

This dataset is temporarily available at: http://chopen.herokuapp.com 


## Ion Channel Modeling Using chModeler and Binomial Distribution

For a detailed introduction, formulation, and interactive use-cases please take a look at [this Jupyter Notebook](intro.ipynb).
Run the following commands manually for a faster execution!


## Installation

To install this project, clone it locally with `git` and run `pip3 install .` inside the directory (make sure you are using python >= 3.4).


## Usage
You can run `python3 chModeler.py -h` for command-line usage help.

There are 3 ways to initiate running chModeler:

* Sending command arguments to chModeler directly: 
```bash
# Build model for ion channel No.5 in dataset and save the model in /data directory using default options.
python3 chModeler.py -i 5

# Do not save results and show plots for each fitting steps (True/1 and False/0 can be used alternatively) 
python3 chModeler.py -i 5 -s False -p True

# Use previous models to fit the model with r2 score threshold of 0.98 (increase r2 (e.g. 0.999) for a smaller but faster initial state.) 
python3 chModeler.py -i 5 -ft 2 -r2 0.98

# Use the model located in `data/2ndFit/18_Sh-B1_DROME_10p_2.json` to fit and plot only final results 
python3 chModeler.py -i 5 -ft 2 -s 0 -fp 1 -mf "data\2ndFit\18_Sh-B1_DROME_10p_2.json"
```   
* Reading command arguments from a file: 
```bash 
# Read the arguments from args.txt file
python3 chModeler.py "@args.txt" 
```

* Running the wizard for answering questions and entering the options:
```bash 
# Run the wizard
python3 chModeler.py -w 
```


## TODOs:

 * **Model Validation**: Because these models are based on digitized data from literature, voltage traces are dependent to the related experiment. So. in order to verify each model to a wide range of voltages there is a need to validate final models against a wider range of voltages.
 * **Smarter Initial Parameter Selection**: Right now, previous models for next fittings are being selected based on R2 scores. This can be improved to select models from the same category using some classification technique.
 * **Temperature and Calcium concentration dependency**: Although the *c* constant can cover the calcium concentration dependency, it's needed to be validated.
 * **Standard data support**: This package needs to be updated to support data analysis from patch clamp machines.
 