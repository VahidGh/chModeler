'''
Initiating ion channel model fitting
'''

import os
from sys import argv
import argparse
import pickle
import json
import copy
from datetime import datetime
import matplotlib.pyplot as plt
# import pandas as pd
# import plotly.plotly as py
# import cufflinks as cf
from fitter import analyze_voltage_clamp
from fetcher import get_data_via_api

PATH = os.getcwd() + '/data/'
if not os.path.exists(PATH):
    os.makedirs(PATH)
R2LIMIT = 0
ARGS_FILE = 'args.txt'


def fit_digitize(params=None):

    if params is not None:
        args = args_parse(args=[])
        for key in params.keys():
            args[key] = params[key]
    else:
        args = args_parse()
    print('\nParameters: ', args)

    ds_file = PATH + 'dataset.pickle'
    save_path = PATH
    dataset = pickle.load(open(ds_file, 'rb'), encoding='latin1')

    # if initArgs == '-w':
    if hasattr(args, 'wizard') and args['wizard'] == True:
        fit_type = input('\nPlease specify fitting type (1): \n1: Blind Guess(quick) \n2: Using previous models\n')
        if fit_type == '':
            fit_type = 1
        elif int(fit_type) == 2:
            ds_type = input('\nPlease specify the dataset (1): \n1: Default \n2: From web \n3: From file\n')
            if ds_type == '':
                ds_type = 1
            elif int(ds_type) == 2:
                ds_host = input('\nPlease enter API host ("http://127.0.0.1:8000"):')
                if ds_host == '':
                    ds_host = 'http://127.0.0.1:8000'
                dataset = get_data_via_api(plot=False, host=ds_host)
                pickle.dump(dataset, open(PATH + 'dataset_new.pickle', 'wb'))
                ds_file = PATH + 'dataset_new.pickle'
                dataset = pickle.load(open(ds_file, 'rb'), encoding='latin1')
            elif int(ds_type == 3):
                ds_file = input('\nPlease enter dataset file path ("data/dataset.pickle"):')
                if ds_file == '':
                    ds_file = PATH + 'dataset.pickle'
                    dataset = pickle.load(open(ds_file, 'rb'), encoding='latin1')
            R2LIMIT = input('\nPlease enter r2 score threshold for model selection (0.0):')
            if R2LIMIT == '':
                R2LIMIT = 0.0
            else:
                R2LIMIT = float(R2LIMIT)
            model_file = input(
                '\nPlease enter model file path if there is any specific model to initiate parameters from (""):')

        ion_channel_id = input('\nPlease enter ion channel id:')
        if ion_channel_id == '':
            ion_channel_id = 0
        else:
            ion_channel_id = int(ion_channel_id)
        time_id = input('\nPlease enter ID of time axis in dataset to start modeling from (1):')
        if time_id == '':
            time_id = 1
        else:
            time_id = int(time_id)
        plot = input('\nSpecify if plots shall be shown during each step of fitting (False):')
        if plot == '':
            plot = False
        else:
            plot = eval(plot)
        final_plot = input('\nSpecify if only final plots shall be shown (False):')
        if final_plot == '':
            final_plot = False
        else:
            final_plot = eval(final_plot)
        save = input('\nSpecify if results shall be saved (True):')
        if save == '':
            save = True
        else:
            save = eval(save)
        if save is True:
            save_path = input('\nPlease specify where the results shall be saved (data/):')
            if save_path == '':
                save_path = PATH
    else:
        # args = parser.parse_args()
        fit_type = int(args['fit_type'])
        if fit_type == 2:
            ds_type = int(args['dataset_type'])
            if ds_type == 2:
                ds_host = args['ds_host']
                if ds_host == '':
                    ds_host = 'http://127.0.0.1:8000'
                dataset = get_data_via_api(plot=False, host=ds_host)
                pickle.dump(dataset, open(PATH + 'dataset_new.pickle', 'wb'))
                ds_file = PATH + 'dataset_new.pickle'
                dataset = pickle.load(open(ds_file, 'rb'), encoding='latin1')
            elif ds_type == 3:
                ds_file = args['ds_file']
                if ds_file == '':
                    ds_file = PATH + 'dataset.pickle'
                dataset = pickle.load(open(ds_file, 'rb'), encoding='latin1')
            R2LIMIT = float(args['r2'])
            model_file = args['model_file']
        ion_channel_id = int(args['ion_channel_id'])
        time_id = int(args['time_id'])
        plot = eval(str(args['plot']))
        final_plot = eval(str(args['final_plot']))
        save = eval(str(args['save']))
        if save is True:
            save_path = args['save_path']
            if save_path == '':
                save_path = PATH

    ich = int(ion_channel_id)
    args = {}
    args['current_type'] = 1
    args['onset_id'] = time_id
    args['plot'] = plot
    args['full_plot'] = final_plot
    args['save_plot'] = save

    if args['save_plot'] is True:
        global PATH2
        PATH2 = save_path + str(datetime.now().strftime('%Y%m%d%H%M%S%f')) + '/'
        print(PATH2)
        if not os.path.exists(PATH2):
            os.makedirs(PATH2)
        args['save_path'] = PATH2
    # args['R2LIMIT'] = R2LIMIT

    print('fit type: ', fit_type)
    print('Args: ', args)
    print('Ion channel id: ', dataset[ich]['ion_channel']['id'])
    print('Digitized graph id: ', dataset[ich]['graph']['id'])
    print('\nDataset id: ', ich)

    args_all = [copy.deepcopy(args)] * len(dataset)
    args_all[ich] = args

    # TODO: Correct in dataset
    no_digitize = [27, 54, 63, 89, 142, 143, 144, 145, 147, 148, 152, 153, 154, 200, 201]
    digitize_error = [32, 33, 43, 44, 45, 46, 160, 161, 162, 163]
    init_curr = [112, 122, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 195, 199, 202, 205]

    if dataset[ich]['graph']['id'] in no_digitize:
        # continue
        raise ValueError('No Digitized data available for this id!')
    if dataset[ich]['graph']['id'] in digitize_error:
        for trace in dataset[ich]['traces']:
            trace['x'].append(trace['x'][-1])
            trace['y'].append(trace['y'][-1])
    if dataset[ich]['graph']['id'] in init_curr:
        for trace in dataset[ich]['traces']:
            trace['x'][1] = trace['x'][2]
            trace['y'][1] = trace['y'][2]
    if dataset[ich]['graph']['id'] == 50:
        dataset[ich]['patch_clamp']['holding_potential'] = 0.0

    if fit_type == 1:
        model_data, plot_data = analyze_voltage_clamp(dataset[ich], args=args_all[ich])
    else:
        inits = []
        if model_file:
            print('Found param file: ', model_file)
            with open(model_file) as f:
                resf = json.load(f)
                inits.append(resf)
        else:
            for root, dirs, files in os.walk(PATH):
                for file in files:
                    if file.endswith(".json"):
                        # if file.endswith("_result.json"):
                        # if file.endswith("10p_2.json"):
                        fi = os.path.join(root, file)
                        print('Found param file: ', fi)
                        with open(fi) as f:
                            resf = json.load(f)
                            if resf['error']['r2'] > R2LIMIT:
                                inits.append(resf)
        print('\nTotal initial parameters: ', len(inits))
        args_all[ich]['init_params'] = inits
        model_data, plot_data = analyze_voltage_clamp(dataset[ich], args=args_all[ich])
        # for i, data in enumerate(dataset):
            # if i == ich:
            # if i>-1:
                # print('\nDataset id: ', i)
                # args_all[i]['init_params'] = inits
                # model_data, plot_data = analyze_voltage_clamp(dataset[i], args=args_all[i])
    # if result['error']['r2'] > R2LIMIT and save is True:
    if save is True:
        fname = PATH2 + str(dataset[ich]['graph']['id']) + \
                '_' + dataset[ich]['ion_channel']['channel_name'] + '_model.json'
        with open(fname, 'w') as f:
            # TODO: multiple json should be in [] and seperate by ,
            json.dump(eval(str(model_data)), f, indent=4)
        print('\nResults successfully saved as: ', fname)
    
    # plot_data.show()
    return model_data, plot_data


def args_parse(args=None):

    parser = argparse.ArgumentParser(fromfile_prefix_chars='@',
                                     description='Build models for ion channels from patch clamp data.')
    parser.add_argument('-w', '--wizard', action="store_true",
                        help="Run wizard for modeler options.")
    parser.add_argument('-ft', '--fit_type', type=int, default=1,
                        help="fitting type, 1: Blind Guess(quick), 2: Using previous models")
    parser.add_argument('-dt', '--dataset_type', type=int, default=1,
                        help="dataset, 1: Default existing digitized data, 2: From web, 3: From file")
    parser.add_argument('-dh', '--ds_host', default="http://127.0.0.1:8000",
                        help="dataset host")
    parser.add_argument('-df', '--ds_file', default="data/dataset.pickle",
                        help="dataset file path")
    parser.add_argument('-cw', '--channelworm', action="store_true",
                        help="Select ion channels from ChannelWorm.")
    parser.add_argument('-mf', '--model_file', default="",
                        help="An specific model file path to initiate parameters from.")
    parser.add_argument('-r2', '--r2', type=float, default=0.0,
                        help="r2 score threshold for model selection")
    parser.add_argument('-i', '--ion_channel_id', type=int, default=0,
                        help="Ion Channel ID from dataset to model")
    parser.add_argument('-ti', '--time_id', type=int, default=1,
                        help="ID of time axis in dataset to start modeling from")
    parser.add_argument('-p', '--plot', default=False,
                        help="If TRUE, all step plots will be shown.")
    parser.add_argument('-fp', '--final_plot', default=False,
                        help="If TRUE, final plots will be shown.")
    parser.add_argument('-s', '--save', default=True,
                        help="If TRUE, results will be saved.")
    parser.add_argument('-sp', '--save_path', default='data/',
                        help="Path for saving final results.")    
    
    if args is None:
        return vars(parser.parse_args())
    else:
        return vars(parser.parse_args(args=args))

if __name__ == '__main__':

    args = args_parse()
    fit_digitize(params=args)
