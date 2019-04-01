'''
Fetch all ion channel data from db and return a single dataset
'''

from matplotlib import pyplot as plt
import numpy as np


def get_data_via_api(fig_id=None, adjust=None, plot=False, host=None):

    import requests

    if host is None:
        HOST = 'http://127.0.0.1:8000'
    else:
        HOST = host
    DETAILS_API_URL = HOST + '/api/{}/details/{}/?format=json'
    LIST_API_URL = HOST + '/api/{}/list/?format=json'

    i_type = ['I', 'I_ss', 'I_peak', 'Current', 'Steady-state Current', 'Peak Current']
    t_type = ['T', 'Time']
    v_type = ['V', 'Voltage']
    po_type = ['I_norm', 'Normalized Current', 'G/G_max', 'Po_peak', 'Peak Open Probability', 'Po', 'Open Probability',
               'G', 'Conductance']
    label = None
    sample_type = None
    dataset = []
    interp_range = 1000

    if fig_id:
        graphs = [requests.get(DETAILS_API_URL.format('graph',fig_id)).json()]
    else:
        graphs = requests.get(LIST_API_URL.format('graph')).json()

    for graph in graphs:
        print('\nGraph: {}'.format(graph['id']))
        # Only IT curves for now
        if not (graph['x_axis_type'] in t_type and graph['y_axis_type'] in i_type):
            continue

        graph_data_all = requests.get(LIST_API_URL.format('graph_data')).json()
        graph_data = [g for g in graph_data_all if g['graph'] == graph['id']]
        ion_channel_all = requests.get(LIST_API_URL.format('ion_channel')).json()
        for ic in ion_channel_all:
            if ic['id'] == graph['ion_channel']:
                ion_channel = ic
        # ion_channel = requests.get(DETAILS_API_URL.format('ion_channel',graph['id'])).json()
        print('Ion Channel: {}, {}'.format(ion_channel['id'], ion_channel['channel_name']))
        curator_all = requests.get(LIST_API_URL.format('user')).json()
        for c in curator_all:
            if c['id'] == ion_channel['username']:
                curator = c
        # patch_clamp = requests.get(DETAILS_API_URL.format('patch_clamp',graph['id'])).json()
        patch_clamp_all = requests.get(LIST_API_URL.format('patch_clamp')).json()
        for pc in patch_clamp_all:
            if pc['id'] == graph['patch_clamp']:
                patch_clamp = pc
        print('Patch_clamp: {}'.format(patch_clamp['id']))
        cell_all = requests.get(LIST_API_URL.format('cell')).json()
        for p in cell_all:
            if p['id'] == patch_clamp['cell']:
                cell = p
        print('cell: {}'.format(cell['cell_type']))
        # reference = requests.get(DETAILS_API_URL.format('reference',graph['id'])).json()
        reference_all = requests.get(LIST_API_URL.format('reference')).json()
        for rf in reference_all:
            if rf['id'] == graph['reference']:
                reference = rf
        doi = reference['doi']
        print('Reference: {}, {}'.format(reference['id'], reference['citation'].encode('utf-8')))

        fig_ref = {'fig': graph['figure_ref_address'], 'doi': doi}
        x_var = {'type': graph['x_axis_type'], 'unit': graph['x_axis_unit'], 'toSI': graph['x_axis_toSI']}
        y_var = {'type': graph['y_axis_type'], 'unit': graph['y_axis_unit'], 'toSI': graph['y_axis_toSI']}
        graph_dic = {'fig_ref': fig_ref, 'x_var': x_var, 'y_var': y_var, 'traces': [], 'ion_channel': ion_channel,
                     'patch_clamp': patch_clamp, 'cell': cell, 'graph': graph, 'curator': curator, 'ref': reference}

        # TODO: deactivation (after offset)
        xps = []
        yps = []
        ys = []
        min_xp = 1e100
        max_xp = -1e100
        for ind in range(len(graph_data)):
            obj = graph_data[ind]
            xp, yp = series2array(obj['series_data'])
            for i in range(len(xp)):
                if adjust and 'x' in adjust:
                    xp[i] += adjust['x']
                xp[i] *= x_var['toSI']
                if adjust and 'y' in adjust:
                    yp[i] += adjust['y']
                yp[i] *= y_var['toSI']
            if min(xp) < min_xp:
                min_xp = min(xp)
            if max(xp) > max_xp:
                max_xp = max(xp)
            xps.append(xp)
            yps.append(yp)
        x = np.linspace(min_xp, max_xp, interp_range)

        for ind in range(len(graph_data)):
            obj = graph_data[ind]
            ys.append(np.interp(x, xp, yps[ind]))
            if graph['x_axis_type'] in t_type and graph['y_axis_type'] in i_type:
                sample_type = 'VoltageClamp'
                graph_dic['traces'].append({'vol': int(obj['series_name']) * 1e-3, 'x': xps[ind], 'y': yps[ind],
                                            'y_interp': ys[ind]})
                label = obj['series_name']

            elif graph['x_axis_type'] in t_type and graph['y_axis_type'] in v_type:
                sample_type = 'CurrentClamp'
                graph_dic['traces'].append({'amp': int(obj['series_name']) * 1e-12, 'x': xps[ind], 'y': yps[ind],
                                            'y_interp': ys[ind]})
                label = obj['series_name']

            elif graph['x_axis_type'] in v_type and graph['y_axis_type'] in i_type:
                sample_type = 'IV'
                graph_dic['V'] = xp
                if graph['y_axis_type'] == 'I_peak':
                    graph_dic['I_peak'] = yps[ind]
                else:
                    graph_dic['I'] = yps[ind]
                label = obj['series_name']

            elif graph['x_axis_type'] in v_type and graph['y_axis_type'] in po_type:
                sample_type = 'POV'
                graph_dic['V'] = xp
                if graph['y_axis_type'] == 'Po_peak':
                    graph_dic['PO_peak'] = yps[ind]
                else:
                    graph_dic['PO'] = yps[ind]
                label = obj['series_name']

            if plot:
                plt.plot(xp, yps[ind], 'o')
                plt.plot(x, ys[ind], '--', label=label)

        graph_dic['interpolated_time'] = x
        graph_dic['sample_type'] = sample_type
        dataset.append(graph_dic)
        if plot:
            plt.title('Raw data for %s from Fig.%s, DOI: %s' % (ion_channel['channel_name'], fig_ref['fig'], fig_ref['doi']))
            plt.xlabel('%s (%s)' % (x_var['type'], x_var['unit']))
            plt.ylabel('%s (%s)' % (y_var['type'], y_var['unit']))
            plt.margins(x=0.1, y=0.1)
            plt.legend(bbox_to_anchor=(1.01, 0.25, 0.2, 0), loc=3, mode="expand", borderaxespad=0., fontsize=8)
            # plt.show()
            if graph.file:
                import matplotlib.image as mpimg
                img = mpimg.imread(graph.file)
                pltimg = plt.figure(2)
                plt.imshow(img)
            plt.show()

    return dataset


def series2array(series):

    xy = series.splitlines()
    data = list()
    for row in xy:
        data += [map(float, row.split(','))]
    xp, yp = list(map(list, zip(*data)))

    return xp, yp

