"""
Modeling ion channels from digitized patch clamp data
"""

from matplotlib import pyplot as plt
from matplotlib import gridspec
import matplotlib.image as mpimg
import numpy as np
from scipy.optimize import fmin, curve_fit, least_squares
import copy
from datetime import datetime
import pandas as pd
from scipy.optimize import differential_evolution
np.seterr(over='ignore')
# import plotly.plotly as py
import plotly.graph_objs as go
import plotly.tools as tls
from plotly.offline import plot


def boltzmannerr(params, input_, actual_output, kwargs, full_output=False):

    mu = params[0]
    k = params[1]
    npower = 1
    kpower = 1
    if len(params) < 3:
        g = 1
        e_rev = 0
    else:
        e_rev = params[2]
        if len(params) > 3:
            g = params[3]
    if "area" in kwargs:
        g = kwargs["area"]

    fitted_curve = boltzmann(input_, mu, k, g, e_rev, npower, kpower, kwargs=kwargs)
    if fitted_curve.size == 0:
        sse = 1e9
    else:
        error_vector = fitted_curve - actual_output
        sse = np.sum(error_vector * error_vector) / max(actual_output) ** 2
        # sse = error_vector / max(np.abs(actual_output))
    if full_output:
        return sse, fitted_curve
    else:
        return sse


def boltzmannfit(
    x_in, y_in, x0=None, plot=False, kwargs=None, bounds=(-np.inf, np.inf)
):
    if x0:
        starting = x0
    else:
        starting = [0.01, 0.01, -0.091, 1e-6]
    args = (x_in, y_in, kwargs)

    try:
        if "x_scale" in kwargs:
            x_scale = kwargs["x_scale"]
        else:
            x_scale = 1.0
        # best_fit = least_squares(boltzmannerr, starting, bounds=bounds, args=args, x_scale=x_scale)['x']

        bounds_map = [[bounds[0][i], bounds[1][i]] for i in range(len(bounds[0]))]
        starting = np.asarray(starting)
        dim = starting.shape
        # print(dim)
        # print(starting)
        if x0 is None:
            starting = "latinhypercube"
        elif starting.ndim > 1:
            if dim[0] < 6:
                # print(1, dim[0])
                startings = list(copy.deepcopy(starting))
                for i in range(dim[0], 6):
                    startings.append(starting[0])
                starting = copy.deepcopy(startings)
                # print(len(starting))
        else:
            # print(2, dim)
            startings = []
            for i in range(1, 6):
                startings.append(starting)
            starting = copy.deepcopy(startings)
        # print(starting)
        best_fit = differential_evolution(
            boltzmannerr,
            init=starting,
            bounds=bounds_map,
            args=args,
            updating="deferred",
            workers=-1,
        )["x"]
    except Exception as e:
        print("\nException in boltzmannfit!!")
        print(e)
        best_fit = starting
        pass

    mu = best_fit[0]
    k = best_fit[1]
    npower = 1
    kpower = 1
    if len(best_fit) > 2:
        e_rev = best_fit[2]
        if len(best_fit) > 3:
            g = best_fit[3]
    else:
        g = 1
        e_rev = 0

    if "area" in kwargs:
        g = kwargs["area"]
    if "npower" in kwargs:
        npower = int(kwargs["npower"])
    if "kpower" in kwargs:
        kpower = int(kwargs["kpower"])

    title = "Boltzmann fit for {} curve".format(kwargs["type"])
    x_fit = np.linspace(x_in[0], x_in[-1])
    y_fit = boltzmann(x_fit, mu, k, g, e_rev, npower, kpower, kwargs=kwargs)
    y_fit_xin = boltzmann(x_in, mu, k, g, e_rev, npower, kpower, kwargs=kwargs)

    r2, rmse = rsquare(y_in, y_fit_xin)
    fit_error = {"r2": r2, "rmse": rmse}

    if plot:
        plt.plot(np.asarray(x_in) * 1e3, np.asarray(y_in), "o")
        plt.plot(np.asarray(x_fit) * 1e3, np.asarray(y_fit))
        plt.title(title)
        plt.margins(x=0.1, y=0.1)
        plt.legend()
        plt.show()

    return best_fit, x_fit, y_fit, fit_error


def boltzmann(x, mu, s, g=1, e_rev=0, npower=1, kpower=1, n0=0, rate=1, kwargs=None):

    if "area" in kwargs:
        g = kwargs["area"]
    curr = None
    if s == 0:
        s = 1e-4
    n = npower
    k = kpower
    if kwargs and "type" in kwargs:
        type = kwargs["type"]
        if type != "n/V" and type != "h/V":
            if "g" in kwargs:
                g = kwargs["g"]
            if "npower" in kwargs:
                n = int(kwargs["npower"])
            if "kpower" in kwargs:
                k = int(kwargs["kpower"])
            curr = g * (x - e_rev)

    # if -708 < (mu - x)/k < 708:
    try:
        y = np.divide(rate, (1 + np.exp(np.divide((mu - x), s))))
        if curr is not None:
            pqco = np.math.factorial(int(k)) / (
                np.math.factorial(int(k) - int(n)) * np.math.factorial(int(n))
            )
            peak = pqco * (n / k) ** n * ((k - n) / k) ** (k - n)
            return curr * y * peak
        elif kwargs and "I_peak" in kwargs:
            return kwargs["I_peak"] * y
        else:
            return y
    except OverflowError:
        y = np.divide((mu - x), s)
        if y[y > 0].any():
            return 1e-100 * y
        else:
            if curr is not None:
                return curr
            else:
                return rate * y


# def aberr(params, input_, actual_output, ext_args=None, full_output=False):
#     fitted_curve = alphabeta(input_, params, ext_args)
#     if fitted_curve.size == 0:
#         sse = 1e9
#     else:
#         error_vector = fitted_curve - actual_output
#         # sse = error_vector / max(np.abs(actual_output))
#         sse = np.sum(error_vector * error_vector) / max(actual_output)**2
#     if full_output:
#         return sse, fitted_curve
#     else:
#         return sse


# def abfit(x_in, y_in, x0=None, bounds=(-np.inf, np.inf), ext_args=None, plot=False):
#     if x0:
#         starting = x0
#     else:
#         starting = [1, 0.01, 0.01]
#     args = (x_in, y_in, ext_args)

#     try:
#         # best_fit = least_squares(aberr, starting, bounds=bounds, args=args)['x']
#         bounds_map = [[bounds[0][i],bounds[1][i]] for i in range(len(bounds[0]))]
#         starting = np.asarray(starting)
#         dim = starting.shape
#         # print(dim)
#         # print(starting)
#         if starting.ndim>1:
#             if dim[0] < 6:
#                 # print(1, dim[0])
#                 startings = list(copy.deepcopy(starting))
#                 for i in range(dim[0], 6):
#                     startings.append(starting[0])
#                 starting = copy.deepcopy(startings)
#                 # print(len(starting))
#         else:
#             # print(2, dim)
#             startings = []
#             for i in range(1, 6):
#                 startings.append(starting)
#             starting = copy.deepcopy(startings)
#         # print(starting)
#         best_fit = differential_evolution(aberr, init='latinhypercube', bounds=bounds_map, args=args, updating='deferred', workers=-1)['x']
#     except Exception as e:
#         print('\nException in abfit!!')
#         print(e)
#         best_fit = starting
#         pass

#     x_fit = np.linspace(x_in[0], x_in[-1])
#     y_fit = alphabeta(x_fit, best_fit, ext_args)

#     if plot:
#         plt.plot(np.asarray(x_in)*1e3, np.asarray(y_in), 'o')
#         plt.plot(np.asarray(x_fit)*1e3, np.asarray(y_fit))
#         plt.title('alpha beta curve fit')
#         plt.margins(x=0.1, y=0.1)
#         plt.legend()
#         plt.show()

#     return best_fit, x_fit, y_fit


def alphabeta(v, args, kwargs):

    if "type" in kwargs:
        if kwargs["type"] == "sigmoid":
            r = args[0]
            vh = args[1]
            s = args[2]
            return np.divide(r, (1 + np.exp(-(v - vh) / s)))
        else:
            ra = args[0]
            sa = args[1]
            return np.multiply(ra, np.exp(-sa * v))
    else:
        ra = args[0]
        sa = args[1]
        return np.multiply(ra, np.exp(-sa * v))


def rsquare(y, f):
    tmp = np.logical_not(np.logical_or(np.isnan(y), np.isnan(f)))
    y = y[tmp]
    f = f[tmp]
    if len(y) == 0:
        return 0, np.nan
    # r2 = -(1 - np.sum((y - f) ** 2) / np.sum((y - np.mean(y)) ** 2))
    r2 = max(0, 1 - np.sum((y - f) ** 2) / np.sum((y - np.mean(y)) ** 2))
    rmse = np.sqrt(np.mean((y - f) ** 2))
    return r2, rmse


def pt(t, r, n, k=4, ninf=1):
    pt = 1 - np.exp(-r * t)
    prt = ninf * (pt) ** n
    qrt = (1 - (ninf * pt)) ** (k - n)
    p = (
        (np.math.factorial(k) / (np.math.factorial(k - n) * np.math.factorial(n)))
        * prt
        * qrt
    )
    eq = "((({:.3f}*t)^{} / {}!) * exp(-{:.3f}*t))".format(r * 1e3, n, n, r * 1e3)
    return p, eq


def nt(t, inf, zero, tau):
    if tau == 0:
        tau = 1e-99
    n = inf - ((inf - zero) * np.exp(np.divide(-t, tau / inf)))
    eq = "{:.3f} - (({:.3f} - {:.3f}) * exp(-t/{:.2f})))".format(
        inf, inf, zero, tau * 1e3
    )
    return n, eq


def full_trace(t, params, ext_args, type=1):

    volt = ext_args["volt"]
    v_hold = ext_args["v_hold"]
    trace_size = ext_args["trace_size"]
    n = ext_args["npower"]
    k = ext_args["kpower"]
    if "v_half" in ext_args and len(params) < 7:
        v_half = ext_args["v_half"]
    else:
        v_half = params[9]
    if "slope" in ext_args and len(params) < 7:
        slope = ext_args["slope"]
    else:
        slope = params[10]
    if "e_rev" in ext_args and len(params) < 7:
        e_rev = ext_args["e_rev"]
    else:
        e_rev = params[11]
    if "g" in ext_args and len(params) < 7:
        g = ext_args["g"]
    else:
        g = params[12]

    ra = params[0]
    sa = params[1]
    rb = params[2]
    vb = params[3]
    sb = params[4]
    rc = params[5]

    y = []
    for i in range(len(volt)):
        v = volt[i]
        xx = np.asarray(t[trace_size[i] : trace_size[i + 1]])

        a = alphabeta(v, [ra, sa], {})
        b = alphabeta(v, [rb, vb, sb], {"type": "sigmoid"})
        c = rc
        a0 = alphabeta(v_hold, [ra, sa], {})
        b0 = alphabeta(v_hold, [rb, vb, sb], {"type": "sigmoid"})
        c0 = rc

        n_inf = 1 / (1 + (b / a))
        n0 = 1 / (1 + (b0 / a0))
        n_inf2 = 1 / (1 + (c / a))
        n02 = 1 / (1 + (c0 / a0))

        qt1 = n_inf - ((n_inf - n0) * np.exp(-(a + b) * xx))
        qt2 = n_inf2 - ((n_inf2 - n02) * np.exp(-(a + c) * xx))
        qt = qt1 * qt2
        prt = (1 - qt) ** n
        qrt = qt ** (k - n)
        pqco = np.math.factorial(int(k)) / (
            np.math.factorial(int(k) - int(n)) * np.math.factorial(int(n))
        )
        p = pqco * prt * qrt
        I_peak = 1 / (1 + np.exp(-(v - v_half) / slope))
        I = g * I_peak * p * (v - e_rev)

        if "t_ic" in ext_args:
            I += (
                ext_args["rate_ic"]
                * (v - ext_args["v_hold"])
                * np.exp(-xx / ext_args["t_ic"])
            )
            I[0] = 0

        y.extend(I)

    eq_peak = "1/(1+exp(-(v-{:.2f})/{:.2f}))".format(v_half * 1e3, slope * 1e3)
    eq_q1 = "q1_inf-((q1_inf-q1_0)*exp(-(a+b)*t))"
    eq_q2 = "q2_inf-((q2_inf-q2_0)*exp(-(a+c)*t))"
    eq_a = "{:.2f}*exp(-{:.2f}*v)".format(ra, sa)
    eq_b = "{:.2f}/(1+exp(-(v-{:.2f})/{:.2f}))".format(rb, vb * 1e3, sb * 1e3)
    eq_c = "{:.3f}".format(rc)
    eq_q = "q1*q2"
    eq_p = "1-q"
    eq_pqco = "{}!/(({}-{})!*{}!)".format(k, k, n, n)
    str_eq = (
        "{:.2E}".format(g)
        + " * "
        + eq_peak
        + " * "
        + eq_pqco
        + " * (p^{})".format(n)
        + " * (q^{})".format(k - n)
    )
    str_eq_v = " * (V - {:.2f})".format(e_rev * 1e3)
    str_eq += str_eq_v
    if "t_ic" in ext_args:
        str_eq += " + ({:.2E} * (V - {:.2f}) * exp(-t/{:.2f}))".format(
            ext_args["rate_ic"], ext_args["v_hold"] * 1e3, ext_args["t_ic"] * 1e3
        )

    str_eq = (
        str_eq
        + "\n"
        + "p = "
        + eq_p
        + ", q = "
        + eq_q
        + ", q1 = "
        + eq_q1
        + ", q2 = "
        + eq_q2
        + "\n"
        + "a = "
        + eq_a
        + ", b = "
        + eq_b
        + ", c = "
        + eq_c
    )

    # TODO: Change if expensive
    for j in range(len(y)):
        if np.abs(y[j]) < 1e-50:
            y[j] = 1e-50
        if np.abs(y[j]) > 1e50:
            y[j] = 1e50

    return np.asarray(y), str_eq


def func_full_trace(params, ext_args, type, input_, actual_output, full_output=False):
    y_fit, _ = full_trace(input_, params, ext_args, type)
    if not np.all(np.isfinite(y_fit)):
        sse = 1e50
    else:
        if "onset" in ext_args:
            onset = ext_args["onset"]
            volt = ext_args["volt"]
            trace_size = ext_args["trace_size"]
            out = []
            fit = []
            for i in range(len(volt)):
                yo = copy.deepcopy(actual_output[trace_size[i] : trace_size[i + 1]])
                out.extend(yo[onset[i] :])
                yf = copy.deepcopy(y_fit[trace_size[i] : trace_size[i + 1]])
                fit.extend(yf[onset[i] :])
            error_vector = np.asarray(out) - np.asarray(fit)
        else:
            error_vector = np.asarray(actual_output) - np.asarray(y_fit)
        sse = np.sum(error_vector * error_vector) / max(actual_output) ** 2
        # sse = error_vector / max(np.abs(actual_output))

    if full_output:
        return sse, np.asarray(y_fit)
    else:
        return sse


def fit_full_trace(x_in, y_in, x0, ext_args, type, bounds=(-np.inf, np.inf)):
    starting = x0
    xtol = 1e-8
    ftol = 1e-8
    gtol = 1e-8
    args = (ext_args, type, x_in, y_in)
    try:
        if "x_scale" in ext_args:
            x_scale = ext_args["x_scale"]
        else:
            x_scale = 1.0
        # best_fit = least_squares(fun=func_full_trace, x0=starting, bounds=bounds, args=args,
        #  xtol=xtol, ftol=ftol, gtol=gtol, x_scale=x_scale, verbose=0)['x']
        bounds_map = [[bounds[0][i], bounds[1][i]] for i in range(len(bounds[0]))]
        starting = np.asarray(starting)
        dim = starting.shape
        # print(dim)
        # print(starting)
        if x0 is None:
            starting = "latinhypercube"
        elif starting.ndim > 1:
            # if dim[1] < 4:
            # 2-trace fit
            # starting = 'latinhypercube'
            if dim[0] < 6:
                # print(1, dim[0])
                startings = list(copy.deepcopy(starting))
                for i in range(dim[0], 6):
                    startings.append(starting[0])
                starting = copy.deepcopy(startings)
                # print(len(starting))
        else:
            # print(2, dim)
            startings = []
            for i in range(1, 6):
                startings.append(starting)
            starting = copy.deepcopy(startings)
        # print(starting)
        # print(len(starting),len(starting[0]), np.size(np.asarray(starting),0))
        best_fit = differential_evolution(
            func_full_trace,
            init=starting,
            bounds=bounds_map,
            args=args,
            updating="deferred",
            workers=-1,
        )["x"]

    except Exception as e:
        print("\nException in fit_full_trace!!")
        print(e)
        best_fit = starting
        pass

    I, str_equation = full_trace(x_in, best_fit, ext_args, type)
    sse, fitted_curve = func_full_trace(
        best_fit, ext_args, type, x_in, y_in, full_output=True
    )

    if "onset" in ext_args:
        onset = ext_args["onset"]
        volt = ext_args["volt"]
        trace_size = ext_args["trace_size"]
        fit = []
        out = []
        for i in range(len(volt)):
            yo = y_in[trace_size[i] : trace_size[i + 1]]
            out.extend(yo[onset[i] :])
            yf = fitted_curve[trace_size[i] : trace_size[i + 1]]
            fit.extend(yf[onset[i] :])
        r2, rmse = rsquare(np.asarray(out), np.asarray(fit))
    else:
        r2, rmse = rsquare(y_in, fitted_curve)
    fit_error = {"r2": r2, "rmse": rmse}

    return I, str_equation, best_fit, fit_error


def analyze_voltage_clamp(data, args=None):

    print("-----------------------------------------------------------------")

    if "graph" in data:
        graph = data["graph"]
        print("Graph ID: {}".format(graph["id"]))
        print("Graph File location: {}".format(graph["file"]))
    if "ion_channel" in data:
        ion_channel = data["ion_channel"]
        print(ion_channel["channel_name"])
        print("Animal: {}".format(ion_channel["animal"]))
    if "patch_clamp" in data:
        pc = data["patch_clamp"]
        v_hold = pc["holding_potential"] * 1e-3
        print("Holding potential: ", v_hold)

    if args is None:
        args = {}
    if "current_type" in args:
        current_type = args["current_type"]
    else:
        current_type = None
    # g (s) = G_dens (s/m2) * area (m2)
    if "area" in args:
        area = args["area"]
    else:
        area = 2e-5
    if "npower" in args:
        npower = args["npower"]
    else:
        npower = 1
    if "onset_id" in args:
        onset_id = args["onset_id"]
    else:
        # onset_id = pc['start_time']
        onset_id = 0
    if "interp" in args:
        interp = args["interp"]
    else:
        interp = False
    if "plot" in args:
        plot = args["plot"]
    else:
        plot = True
    if "full_plot" in args:
        full_plot = args["full_plot"]
    else:
        full_plot = True
    if "save_plot" in args:
        save_plot = args["save_plot"]
    else:
        save_plot = False
    if "save_path" in args:
        PATH2 = args["save_path"]
    else:
        PATH2 = "."
    if "R2LIMIT" in args:
        R2LIMIT = args["R2LIMIT"]
    else:
        R2LIMIT = 0
    # eliminate capacitance current effect
    if "peak_tail_bound" in args:
        peak_tail_bound = args["peak_tail_bound"]
    else:
        peak_tail_bound = (-np.inf, np.inf)

    if graph["file"] and plot:
        try:
            img = mpimg.imread(graph["file"])
            # pltimg = plt.figure(2)
            plt.imshow(img)
            plt.show()
        except:
            pass

    cur_type = 1
    if current_type:
        cur_type = current_type
        print("User Current type: ", current_type)
    elif ion_channel["current_type"]:
        current_type = ion_channel["current_type"]
        print("DB Current type: ", ion_channel["current_type"])
        if current_type == "A-type":
            cur_type = 3
    print("Current type: ", cur_type)

    volt = []
    time_trace = []
    series_trace = []
    series = []
    for trace in data["traces"]:
        volt.append(trace["vol"])
        if min(trace["x"]) < 0:
            trace["x"] = [x-min(trace["x"]) for x in trace["x"]]
        time_trace.append(np.asarray(trace["x"]))
        series_trace.append(np.asarray(trace["y"]))
        series.append(np.asarray(trace["y_interp"]))
        if plot:
            plt.plot(trace["x"], trace["y"], "--o", label=trace["vol"])
    volt = np.asarray(volt)
    if min(data["interpolated_time"]) < 0:
        data["interpolated_time"] = [x-min(data["interpolated_time"]) for x in data["interpolated_time"]]
    time = np.asarray(np.asarray(data["interpolated_time"]))
    series = np.asarray(series)
    time_trace = np.asarray(time_trace)
    series_trace = np.asarray(series_trace)
    print(series_trace.shape)

    if plot:
        plt.legend()
        plt.show()

    print("onset id: ", onset_id)
    offsets = []
    onsets = []
    peak_yy = []
    tail_yy = []

    n_trace, n_pt_total = series.shape
    stim_start = time[0]
    stim_end = time[-1]

    ra_v = copy.deepcopy(volt)
    sa_v = copy.deepcopy(volt)
    rb_v = copy.deepcopy(volt)
    vb_v = copy.deepcopy(volt)
    sb_v = copy.deepcopy(volt)
    rc_v = copy.deepcopy(volt)
    a_v = copy.deepcopy(volt)
    b_v = copy.deepcopy(volt)
    c_v = copy.deepcopy(volt)
    act_tau = copy.deepcopy(volt)
    act_tau2 = copy.deepcopy(volt)
    negatives = copy.deepcopy(volt)
    print("Traces before sort: ", n_trace, volt)

    for i in range(n_trace):
        lent = len(series_trace[i])
        if interp:
            s = series[i]
        else:
            s = series_trace[i]
        t = time_trace[i]

        if lent == 8 or lent == 9 or lent == 11:
            # Activation and Deactivation
            off_t = t[6]
        elif lent == 12 or lent == 13 or lent == 14:
            # Activation, Inactivation and Deactivation
            off_t = t[9]
        else:
            # Activation and/or Inactivation
            off_t = t[-1]
        if interp:
            onset = np.abs(series[0] - series_trace[0][onset_id]).argmin()
            offset = np.abs(time - off_t).argmin() - 1
        else:
            onset = onset_id
            offset = np.abs(t - off_t).argmin()
        offsets.append(offset)
        onsets.append(onset)

        if peak_tail_bound[0] == -np.inf:
            on = onset
        else:
            on = peak_tail_bound[0]
        abs_peak = np.max(np.abs(s[on : offset + 1]))
        abs_peaki = np.argmax(np.abs(s[on : offset + 1]))
        if s[on : offset + 1][abs_peaki] < 0:
            negatives[i] = True
            peak_yy.append(-abs_peak)
            tail_yy.append(-np.abs(s[offset]))
        else:
            negatives[i] = False
            peak_yy.append(abs_peak)
            tail_yy.append(s[offset])

    sort_onsets = sorted(zip(volt, onsets), reverse=True)
    onsets = [p[1] for p in sort_onsets]
    sort_offsets = sorted(zip(volt, offsets), reverse=True)
    offsets = [p[1] for p in sort_offsets]
    point_sort_peak = sorted(zip(volt, peak_yy), reverse=True)
    point_sort_tail = sorted(zip(volt, tail_yy), reverse=True)
    negatives_sort2 = sorted(zip(volt, negatives), reverse=True)
    negatives_sort = np.asarray([p[1] for p in negatives_sort2])
    volt_sort = [p[0] for p in point_sort_peak]

    if not np.array_equal(volt, volt_sort):
        time_trace_sort = np.zeros(shape=time_trace.shape)
        series_trace_sort = np.zeros(shape=series_trace.shape)
        series_sort = np.zeros(shape=series.shape)
        for i, vs in enumerate(volt_sort):
            for j, v in enumerate(volt):
                if vs == v:
                    time_trace_sort[i] = time_trace[j]
                    series_trace_sort[i] = series_trace[j]
                    series_sort[i] = series[j]

        volt = np.asarray(volt_sort)
        time_trace = np.asarray(time_trace_sort)
        series_trace = np.asarray(series_trace_sort)
        series = np.asarray(series_sort)

    peak_yy = [p[1] for p in point_sort_peak]
    tail_yy = [p[1] for p in point_sort_tail]

    print("Traces after sort: ", volt)
    print("Negative Currents: ", negatives_sort)
    peak_yy = np.asarray(peak_yy)
    tail_yy = np.asarray(tail_yy)

    n = npower = 3
    k = kpower = 4

    kwargs = {"type": "I/V", "npower": npower, "kpower": kpower}
    if peak_yy[0] != tail_yy[0]:
        kwargs = {"type": "Ipeak/V", "npower": npower, "kpower": kpower}
    kwargs_n = {"type": "n/V"}
    params_inf = ["v_half", "slope", "e_rev", "g"]
    bounds_boltz = [[-0.2, -0.1, -0.2, 0], [0.2, 0.1, 0.2, 1000]]
    # x_scale_inf = [1e-2, 1e-3, 1e-2, 1e-8]
    # kwargs.update({'x_scale': x_scale_inf})

    if "init_params" not in args:
        max_peaksi = np.abs(peak_yy).argmax()
        vmedi = int(len(volt) / 2)
        print("median volt, peak: ", volt[vmedi], (vmedi, peak_yy[vmedi]))
        if peak_yy[vmedi] < 0:
            min_peaksi = np.abs(peak_yy[: max_peaksi + 1]).argmin()
            e_rev0 = volt[min_peaksi] + 1e-3
        else:
            min_peaksi = np.abs(peak_yy).argmin()
            e_rev0 = volt[min_peaksi] + 1e-3
        g0 = (64 / 27) * np.abs(peak_yy[max_peaksi] / (volt[max_peaksi] - e_rev0))
        # x0_peak = [volt[vmedi], 0.009, e_rev0, g0]
        print(
            "max , min peaks: ",
            (max_peaksi, peak_yy[max_peaksi]),
            (min_peaksi, peak_yy[min_peaksi]),
        )
        # x0_peak = None
        # v_halfs0 = [volt[vmedi], -75e-3, -5e-3, 1e-10, 5e-3, 50e-3]
        v_halfs0 = [-75e-3, -5e-3, 5e-3, 50e-3]
        # slopes0 = [-50e-3, -10e-3, -1e-3, 1e-3, 10e-3, 50e-3]
        slopes0 = [-50e-3, -9e-3, 9e-3, 50e-3]
        # e_revs0 = [e_rev0, -75e-3, 1e-10, 50e-3]
        e_revs0 = [-75e-3, 1e-10, 50e-3]
        # gs0 = [g0, 1e-12, 1e-9, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100]
        gs0 = [g0, 1e-12, 1e-9, 1e-6, 1e-3]
        x0_peaks = list(
            np.array(np.meshgrid(v_halfs0, slopes0, e_revs0, gs0)).T.reshape(-1, 4)
        )
        # estimates_peak, x_fit_peak, y_fit_peak, boltzmann_error = boltzmannfit(volt, peak_yy, x0=x0_peaks,
        #    plot=plot, kwargs=kwargs,
        #    bounds=bounds_boltz)
    else:
        if isinstance(args["init_params"], dict):
            ninit = 0
        else:
            ninit = len(args["init_params"])
            print("\nNumber of initial params: ", ninit)
        if ninit > 0:
            # r2s = np.zeros(ninit)
            # boltzmann_errors = [[]]*ninit
            x0_peaks = [[]]*ninit
            # estimates_peaks = [[]]*ninit
            # x_fit_peaks = [[]]*ninit
            # y_fit_peaks = [[]]*ninit
            for i, init_param in enumerate(args["init_params"]):
                x0_peaks[i] = list(init_param['params'].values())[:4]
                # x0_peaks.append(list(init_param["params"].values())[:4])
            # estimates_peak, x_fit_peak, y_fit_peak, boltzmann_error = boltzmannfit(volt, peak_yy,
            #    x0=x0_peaks,
            #    plot=plot,
            #    kwargs=kwargs,
            #    bounds=bounds_boltz)
            # r2s[i] = boltzmann_errors[i]['r2']
            # maxerri = np.argmax(np.asarray(r2s))
            # estimates_peak = estimates_peaks[maxerri]
            # x_fit_peak = x_fit_peaks[maxerri]
            # y_fit_peak = y_fit_peaks[maxerri]
            # boltzmann_error = boltzmann_errors[maxerri]
            # x0_peak = x0_peaks[maxerri]
            # print('\ntop x0 peak, index: ', maxerri, x0_peak)
            # print('\ntop x0 peak, index: ', x0_peak)
        else:
            x0_peaks = list(args['init_params']['params'].values())[:4]
            # x0_peaks.append(list(args["init_params"]["params"].values())[:4])
            # print('x0 peak: ', x0_peak)

    print(
        "N x0 peaks, x0 peaks[0], x0 peaks[-1]: ",
        len(x0_peaks),
        x0_peaks[0],
        x0_peaks[-1],
    )
    estimates_peak, x_fit_peak, y_fit_peak, boltzmann_error = boltzmannfit(
        volt, peak_yy, x0=x0_peaks, plot=plot, kwargs=kwargs, bounds=bounds_boltz
    )
    print("boltzmann fit r2 error: ", boltzmann_error)
    if "npower" in kwargs:
        npower = int(kwargs["npower"])
    if "kpower" in kwargs:
        kpower = int(kwargs["kpower"])
    nv = boltzmann(volt, estimates_peak[0], estimates_peak[1], kwargs=kwargs_n)
    init_args = dict(zip(params_inf, estimates_peak))

    if plot:
        print("I peak params: ", init_args)
        print("nv: ", nv)
        print("npower: ", npower)
        print("kpower: ", kpower)

        plt.figure()
        plt.plot(x_fit_peak, y_fit_peak, label="I_fit")
        plt.plot(volt, peak_yy, "--o", label="I_fit")
        plt.title("I_inf/v, I/v curves")
        plt.legend()

        plt.figure()
        plt.plot(volt, nv, "--.", label="n/v")
        plt.title("n/v, h/v curves")
        plt.legend()
        plt.show()

    ext_args = {
        "v_half": estimates_peak[0],
        "slope": estimates_peak[1],
        "e_rev": estimates_peak[2],
        "g": estimates_peak[3],
    }
    ext_args["v_hold"] = v_hold
    if "ic" in args:
        ext_args["rate_ic"] = args["ic"]["rate_ic"]
        ext_args["t_ic"] = args["ic"]["t_ic"]
        if plot:
            print("ic: ", ext_args)

    # x_scale = [1, 1, 1, 1e-3, 1, 1]
    x00 = 0.01
    r0 = 10
    mins = 1e-3
    maxv = 0.3
    maxr = 500
    params = ["ra"]
    # x0 = [r0]
    # rs0 = [1e-3, 1e-2, 1e-1, 1, 10, 100]
    rs0 = [1e-3, 1, 50]
    min_params = [mins]
    max_params = [maxr]
    params.append("sa")
    # x0.append(-r0)
    # sas0 = [-1e-3, -1e-2, -1e-1, -1, -10, -100, 1e-3, 1e-2, 1e-1, 1, 10, 100]
    sas0 = [-1e-3, -1, -50, 1e-3, 1, 50]
    min_params.append(-maxr)
    max_params.append(maxr)
    params.append("rb")
    # x0.append(10 * r0)
    min_params.append(mins)
    max_params.append(maxr)
    params.append("vb")
    # x0.append(-x00)
    min_params.append(-maxv)
    max_params.append(maxv)
    params.append("sb")
    # x0.append(x00)
    min_params.append(mins)
    max_params.append(maxv)
    params.append("rc")
    # x0.append(x00)
    min_params.append(0)
    max_params.append(maxr)
    # rcs0 = [0, 1e-10, 1e-1, 1, 10, 100]
    rcs0 = [0, 1e-5, 1, 50]

    if "init_params" not in args:
        x0 = list(np.array(np.meshgrid(rs0, sas0, rs0, v_halfs0, slopes0, rcs0)).T.reshape(-1, 6))
        print(
            "N x0, x0[0], x0[-1]: ",
            len(x0),
            x0[0],
            x0[-1],
        )

        # Eliminate noise effect:
        peak5pc = np.where(np.abs(peak_yy) < (5 * np.max(np.abs(peak_yy)) / 100))[0]
        if peak5pc.size == 0:
            peak5pc = [-1]
        peak10pc = np.where(np.abs(peak_yy) < (10 * np.max(np.abs(peak_yy)) / 100))[0]
        if peak10pc.size == 0:
            peak10pc = [-1]
        peak25pc = np.where(np.abs(peak_yy) < (25 * np.max(np.abs(peak_yy)) / 100))[0]
        if peak25pc.size == 0:
            peak25pc = [-1]
        peak95pc = np.where(np.abs(peak_yy) < (95 * np.max(np.abs(peak_yy)) / 100))[0]
        if peak95pc.size == 0:
            peak95pc = [-1]

        xf = []
        yf = []
        trace_size = [0]
        xf_interp = []
        trace_size_interp = [0]
        time_interp = copy.deepcopy(time)
        peak_all = np.max(np.abs(peak_yy))
        tail_all = np.max(np.abs(tail_yy))
        bounds = (min_params, max_params)

        peaks = []
        onsets_interp = []
        offsets_interp = []
        peaks_interp = []

        if plot:
            print("5% peak: ", peak5pc, volt[peak5pc])
            print("10% peak: ", peak10pc, volt[peak10pc])
            print("25% peak: ", peak25pc, volt[peak25pc])
            print("95% peak: ", peak95pc, volt[peak95pc])
            print("Tail currents: ", tail_yy)
            print("Peak currents: ", peak_yy)
            print("Max peak: ", peak_all)
            print("Max tail: ", tail_all)
            # print("Full trace params, x0, bounds: ", params, x0, bounds)
            print("Full trace params, bounds: ", params, bounds)

        plot = False
        for i in range(n_trace):
            if interp:
                trace = series[i, :]
            else:
                trace = series_trace[i, :]
                time = time_trace[i, :]

            offset = offsets[i]
            onset = onsets[i]
            peak = abs(np.max(trace[onset : offset + 1]))
            id_peak = np.argmax(abs(trace[onset : offset + 1])) + onset
            peaks.append(id_peak)
            tail = trace[offset]
            id_tail = offset
            if plot:
                print("\n")
                print("V: ", volt[i] * 1e3)
                print("Offset and Onset ids: ", offset, onset)
                print(id_peak, peak)
                print(id_tail, tail)

            ext_args["v"] = volt[i]
            # if i == 0:
                # x0_1 = copy.deepcopy(x0)

            fit_t = copy.deepcopy(time[: offset + 1])
            fit_y = copy.deepcopy(trace[: offset + 1])

            xf.extend(fit_t)
            yf.extend(fit_y)
            trace_size.append(trace_size[i] + len(fit_t))
            xf_interp.extend(time_interp)
            trace_size_interp.append(trace_size_interp[i] + len(time_interp))
            if onset == 0:
                onset_interp = 0
            else:
                onset_interp = np.abs(time_interp - time[onset]).argmin()
            if offset == len(time) - 1:
                offset_interp = len(time_interp) - 1
            else:
                offset_interp = np.abs(time_interp - time[offset]).argmin() - 1
            peakterp = np.abs(time_interp - time[id_peak]).argmin() + onset_interp
            onsets_interp.append(onset_interp)
            offsets_interp.append(offset_interp)
            peaks_interp.append(peakterp)
            ext_args["npower"] = npower
            ext_args["kpower"] = kpower
            ext_args["onset"] = onset
            ext_args["offset"] = offset
            ext_args["peak"] = id_peak
            ext_args["tail"] = fit_y[-1]
            ext_args["ninf_peak"] = nv[i]
            # ext_args.update({"x_scale": x_scale})

            if plot:
                # print("x0: ", dict(zip(params, x0)))
                print("ext_args: ", ext_args)

        if plot and full_plot is False:
            plt.title(
                "Raw and analysis data for %s" % (data["ion_channel"]["channel_name"])
            )
            plt.margins(x=0.1, y=0.1)
            plt.legend(
                bbox_to_anchor=(1.01, 0.25, 0.2, 0),
                loc=3,
                mode="expand",
                borderaxespad=0.0,
                fontsize=8,
            )
            plt.show()

            plt.plot(volt, a_v, label="alpha")
            plt.plot(volt, b_v, label="beta")
            plt.plot(volt, c_v, label="c")
            plt.legend()
            plt.show()

        plot = args["plot"]
        # x0_2 = copy.deepcopy(x0_1)
        # x0_2 = None
        x0_2 = x0
        for i in range(int(n_trace / 2)):
            v_1 = volt[i]
            i2 = i + int(n_trace / 2)
            v_2 = volt[i2]
            v2 = [v_1, v_2]

            params_2 = copy.deepcopy(params)
            min_bound_2 = copy.deepcopy(min_params)
            max_bound_2 = copy.deepcopy(max_params)
            bounds_2 = [min_bound_2, max_bound_2]
            if plot:
                print("\nSelected voltages and indices: ", [i, i2], v2)
                # print('x0: ', dict(zip(params_2, x0_2)))
            ext_args_2 = copy.deepcopy(ext_args)
            ext_args_2["volt"] = v2

            x = xf[trace_size[i] : trace_size[i + 1]]
            y = yf[trace_size[i] : trace_size[i + 1]]
            x_2 = xf[trace_size[i2] : trace_size[i2 + 1]]
            y_2 = yf[trace_size[i2] : trace_size[i2 + 1]]
            x2 = x + x_2
            y2 = y + y_2
            x2 = np.asarray(x2)
            y2 = np.asarray(y2)
            ext_args_2["onset"] = [onsets[i], onsets[i2]]
            ext_args_2["offset"] = [offsets[i], offsets[i2]]
            ext_args_2["peak"] = [peaks[i], peaks[i2]]
            ext_args_2["trace_size"] = [0, len(x), len(x) + len(x_2)]
            # ext_args_2["x_scale"] = x_scale

            if plot:
                print("bounds: ", bounds_2)
                print("args: ", ext_args_2)

            I2, str_eq_2, best_fit_2, fit_error_2 = fit_full_trace(
                x2, y2, x0_2, ext_args_2, cur_type, bounds=bounds_2
            )

            # x0_2 = best_fit_2
            ra = best_fit_2[0]
            sa = best_fit_2[1]
            rb = best_fit_2[2]
            vb = best_fit_2[3]
            sb = best_fit_2[4]
            rc = best_fit_2[5]

            ra_v[i] = ra_v[i2] = ra
            sa_v[i] = sa_v[i2] = sa
            rb_v[i] = rb_v[i2] = rb
            vb_v[i] = vb_v[i2] = vb
            sb_v[i] = sb_v[i2] = sb
            rc_v[i] = rc_v[i2] = rc

            a1 = alphabeta(v_1, [ra, sa], {})
            a2 = alphabeta(v_2, [ra, sa], {})
            b1 = alphabeta(v_1, [rb, vb, sb], {"type": "sigmoid"})
            b2 = alphabeta(v_2, [rb, vb, sb], {"type": "sigmoid"})
            c1 = rc
            c2 = rc
            a_v[i] = a1
            a_v[i2] = a2
            b_v[i] = b1
            b_v[i2] = b2
            c_v[i] = c1
            c_v[i2] = c2
            act_tau[i] = 1 / (a1 + b1)
            act_tau[i2] = 1 / (a2 + b2)
            act_tau2[i] = 1 / (a1 + c1)
            act_tau2[i2] = 1 / (a2 + c2)

            if plot:
                print("best fit: ", dict(zip(params_2, best_fit_2)))
                print("Equation:\n", str_eq_2)
                plt.figure()
                x = np.asarray(x)
                y = np.asarray(y)
                x_2 = np.asarray(x_2)
                y_2 = np.asarray(y_2)
                I_1 = np.asarray(I2[: len(x)])
                I_2 = np.asarray(I2[len(x) : len(x) + len(x_2)])
                plt.plot(x * 1e3, y * 1e6, "o", label=volt[i] * 1e3)
                plt.plot(x_2 * 1e3, y_2 * 1e6, "o", label=volt[i2] * 1e3)
                plt.plot(x * 1e3, I_1 * 1e6)
                plt.plot(x_2 * 1e3, I_2 * 1e6)

                title = "Raw and analized data for %s" % (
                    data["ion_channel"]["channel_name"]
                )
                plt.title(title + "\nEquations: %s" % (str_eq_2), fontsize=8)
                plt.margins(x=0.1, y=0.1)
                plt.legend(
                    bbox_to_anchor=(1.01, 0.25, 0.2, 0),
                    loc=3,
                    mode="expand",
                    borderaxespad=0.0,
                    fontsize=8,
                )
                # plt.show()

                v_f = np.arange(-0.1, 0.1, 0.001)
                tau_args = {
                    "v_half": ext_args["v_half"],
                    "slope": ext_args["slope"],
                    "v_hold": v_hold,
                }
                plt.figure()
                alpha_f = alphabeta(v_f, args=[ra, sa], kwargs=tau_args)
                beta_f = alphabeta(v_f, [rb, vb, sb], {"type": "sigmoid"})
                c_f = [rc] * len(v_f)
                plt.plot(v_f * 1e3, alpha_f, label="alpha")
                plt.plot(v_f * 1e3, beta_f, label="beta")
                plt.plot(v_f * 1e3, c_f, label="c")
                plt.title("alpha-beta/v")
                plt.margins(x=0.1, y=0.1)
                plt.legend()

                plt.figure()
                tau_f = 1 / (alpha_f + beta_f)
                tau_f2 = 1 / (alpha_f + c_f)
                ninf_f = alpha_f / (alpha_f + beta_f)
                ninf_f2 = alpha_f / (alpha_f + c_f)
                plt.plot(v_f * 1e3, tau_f * 1e3, label="tau")
                plt.plot(v_f * 1e3, tau_f2 * 1e3, label="tau2")
                plt.title("tau/v")
                plt.margins(x=0.1, y=0.1)
                plt.legend()

                plt.figure()
                n = ext_args_2["npower"]
                k = ext_args_2["kpower"]
                v_half = ext_args_2["v_half"]
                slope = ext_args_2["slope"]
                I_inf = 1 / (1 + np.exp(-(v_f - v_half) / slope))
                plt.plot(v_f * 1e3, I_inf, label="I_inf")
                plt.plot(v_f * 1e3, ninf_f, label="ninf1")
                plt.plot(v_f * 1e3, ninf_f2, label="ninf2")
                plt.plot(v_f * 1e3, (ninf_f * ninf_f2) ** (k - n), label="ninf")
                plt.plot(v_f * 1e3, (1 - ninf_f * ninf_f2) ** n, label="1-ninf")
                pqco = np.math.factorial(int(k)) / (
                    np.math.factorial(int(k) - int(n)) * np.math.factorial(int(n))
                )
                p_inf_f = (
                    pqco * (1 - ninf_f * ninf_f2) ** n * (ninf_f * ninf_f2) ** (k - n)
                )
                plt.plot(v_f * 1e3, p_inf_f, label="pinf")
                plt.plot(
                    v_f * 1e3, (1 - ninf_f * ninf_f2) ** n * I_inf, label="ninf*I_inf"
                )
                plt.plot(v_f * 1e3, p_inf_f * I_inf, label="pinf*I_inf")
                plt.title("n_inf")
                plt.margins(x=0.1, y=0.1)
                plt.legend()
                plt.show()

        act_tau = np.asarray(act_tau)
        act_tau2 = np.asarray(act_tau2)
        a_v = np.asarray(a_v)
        b_v = np.asarray(b_v)

        atau_params = params[:2]
        min_atau = min_params[:2]
        max_atau = max_params[:2]
        # tau_args_exp = {'type':'exp'}
        # x0_atau = [ra_v[1], sa_v[1]]
        # x0_atau = [[ra_v[i], sa_v[i]] for i in range(len(ra_v))]
        # best_fit_atau, x_atau, y_atau = abfit(volt, a_v, x0_atau,
        #   bounds=(min_atau, max_atau), ext_args=tau_args_exp, plot=plot)

        btau_params = params[2:5]
        min_btau = min_params[2:5]
        max_btau = max_params[2:5]
        # tau_args_sig = {'type': 'sigmoid'}
        # x0_btau = [[rb_v[i], vb_v[i], sb_v[i]] for i in range(len(rb_v))]
        # best_fit_btau, x_btau, y_btau = abfit(volt, b_v, x0_btau,
        #                                       bounds=(min_btau, max_btau), ext_args=tau_args_sig, plot=plot)

        ctau_params = [params[5]]
        min_ctau = [min_params[5]]
        max_ctau = [max_params[5]]
        # x0_ctau = [1e-1]
        # x0_ctau = rc_v
        # y_ctau = np.asarray([1e-1]*len(y_atau))
        # y_ctau = rc_v
        # best_fit_ctau = [1e-1]
        # best_fit_ctau = rc_v

        # if plot:
        # print('\n')
        # print('alpha initial parameters: ', dict(zip(atau_params, x0_atau)))
        # print('alpha parameters: ', dict(zip(atau_params, best_fit_atau)))
        # print('beta initial parameters: ', dict(zip(btau_params, x0_btau)))
        # print('beta parameters: ', dict(zip(btau_params, best_fit_btau)))
        # print('c initial parameters: ', dict(zip(ctau_params, x0_ctau)))
        # print('c parameters: ', dict(zip(ctau_params, best_fit_ctau)))

        # plt.plot(x_atau,1/(y_atau+y_btau), label='ab_tau_fit')
        # plt.plot(x_atau,1/(y_atau+y_ctau), label='ac_tau_fit')
        # plt.plot(volt,act_tau, '--o', label='tau')
        # plt.plot(volt,act_tau2, '--o', label='tau2')
        # plt.legend()
        # plt.show()

        x0_f_params = []
        # x0_f = []
        x0_f = [
            [ra_v[i], sa_v[i], rb_v[i], vb_v[i], sb_v[i], rc_v[i]]
            for i in range(len(ra_v))
        ]
        min_bound_f = []
        max_bound_f = []

        x0_f_params.extend(atau_params)
        # x0_f.extend(best_fit_atau)
        min_bound_f.extend(min_atau)
        max_bound_f.extend(max_atau)

        x0_f_params.extend(btau_params)
        # x0_f.extend(best_fit_btau)
        min_bound_f.extend(min_btau)
        max_bound_f.extend(max_btau)

        x0_f_params.extend(ctau_params)
        # x0_f.extend(best_fit_ctau)
        min_bound_f.extend(min_ctau)
        max_bound_f.extend(max_ctau)

        # x_scale_f = x_scale
        ext_args_f = {
            "volt": list(volt),
            "v_hold": v_hold,
            "trace_size": list(trace_size),
        }

        bounds_f = (min_bound_f, max_bound_f)
        ext_args_f.update(ext_args)
        xf = np.asarray(xf)
        yf = np.asarray(yf)
        xf_interp = np.asarray(xf_interp)
        ext_args_f["onset"] = list(onsets)
        ext_args_f["offset"] = list(offsets)
        ext_args_f["peak"] = list(peaks)
        # ext_args_f["x_scale"] = list(x_scale_f)

    else:
        # x0_f = x0
        xf = []
        yf = []
        trace_size = [0]
        xf_interp = []
        trace_size_interp = [0]
        time_interp = copy.deepcopy(time)
        peaks = []
        onsets_interp = []
        offsets_interp = []
        peaks_interp = []

        for i in range(n_trace):
            if interp:
                trace = series[i, :]
            else:
                trace = series_trace[i, :]
                time = time_trace[i, :]

            offset = offsets[i]
            onset = onsets[i]
            id_peak = np.argmax(abs(trace[onset : offset + 1])) + onset
            peaks.append(id_peak)

            fit_t = copy.deepcopy(time[: offset + 1])
            fit_y = copy.deepcopy(trace[: offset + 1])

            xf.extend(fit_t)
            yf.extend(fit_y)
            trace_size.append(trace_size[i] + len(fit_t))
            xf_interp.extend(time_interp)
            trace_size_interp.append(trace_size_interp[i] + len(time_interp))
            if onset == 0:
                onset_interp = 0
            else:
                onset_interp = np.abs(time_interp - time[onset]).argmin()
            if offset == len(time) - 1:
                offset_interp = len(time_interp) - 1
            else:
                offset_interp = np.abs(time_interp - time[offset]).argmin() - 1
            peakterp = np.abs(time_interp - time[id_peak]).argmin() + onset_interp
            onsets_interp.append(onset_interp)
            offsets_interp.append(offset_interp)
            peaks_interp.append(peakterp)

        if "bounds" in args["init_params"]:
            bounds_f = args["init_params"]["bounds"]
        # else:
            # mins = 1e-3
            # maxv = 0.3
            # maxr = 500

            # min_params = [mins]
            # max_params = [maxr]
            # min_params.append(-maxr)
            # max_params.append(maxr)
            # min_params.append(mins)
            # max_params.append(maxr)
            # min_params.append(-maxv)
            # max_params.append(maxv)
            # min_params.append(mins)
            # max_params.append(maxv)
            # min_params.append(0)
            # max_params.append(maxr)
            # bounds_f = (min_params, max_params)

        ext_args_f = {
            "volt": list(volt),
            "v_hold": v_hold,
            "trace_size": list(trace_size),
        }
        ext_args_f["onset"] = list(onsets)
        ext_args_f["offset"] = list(offsets)
        ext_args_f["peak"] = list(peaks)
        ext_args_f["npower"] = npower
        ext_args_f["kpower"] = kpower
        # x_scale_f = [1, 1, 1, 1e-3, 1, 1]
        # ext_args_f["x_scale"] = list(x_scale_f)
        ext_args_f.update(ext_args)

        if ninit > 0:
            r2_fs = np.zeros(ninit)
            fit_error_fs = [[]] * ninit
            x0_f = [[]] * ninit
            x0_f_paramss = [[]] * ninit
            best_fit_fs = [[]] * ninit
            str_eq_fs = [[]] * ninit
            I_fs = [[]] * ninit
            for i, init_param in enumerate(args["init_params"]):
                if "bounds" in init_param:
                    bounds_f = init_param["bounds"]
                x0_f[i] = list(init_param["params"].values())[4:]
                # x0_f.append(list(init_param["params"].values())[4:])
                x0_f_paramss[i] = list(init_param["params"].keys())[4:]
            # I_f, str_eq_f, best_fit_f, fit_error_f = fit_full_trace(
                # xf, yf, x0_fs, ext_args_f, cur_type, bounds=bounds_f
            # )
            # r2_fs[i] = fit_error_fs[i]['r2']
            # maxerrfi = np.argmax(np.asarray(r2_fs))
            # best_fit_f = best_fit_fs[maxerrfi]
            # str_eq_f = str_eq_fs[maxerrfi]
            # I_f = I_fs[maxerrfi]
            # fit_error_f = fit_error_fs[maxerrfi]
            # x0_f = x0_fs[maxerrfi]
            # x0_f_params = x0_f_paramss[maxerrfi]
            x0_f_params = x0_f_paramss[0]
            # if 'bounds' in args['init_params'][maxerrfi]:
            # bounds_f = args['init_params'][maxerrfi]['bounds']
            # print('\ntop x0 final, index: ', maxerrfi, x0_f)

        else:
            x0_f = list(args["init_params"]["params"].values())[4:]
            # x0_f.append(list(args["init_params"]["params"].values())[4:])
            x0_f_params = list(args["init_params"]["params"].keys())[4:]
            # I_f, str_eq_f, best_fit_f, fit_error_f = fit_full_trace(
                # xf, yf, x0_f, ext_args_f, cur_type, bounds=bounds_f
            # )
    # else:
        # I_f, str_eq_f, best_fit_f, fit_error_f = fit_full_trace(
            # xf, yf, x0_f, ext_args_f, cur_type, bounds=bounds_f
        # )

    print(
        "N x0_f, x0_f[0], x0_f[-1]: ",
        len(x0_f),
        x0_f[0],
        x0_f[-1],
    )
    
    I_f, str_eq_f, best_fit_f, fit_error_f = fit_full_trace(
        xf, yf, x0_f, ext_args_f, cur_type, bounds=bounds_f
    )

    # if plot or full_plot:
    print("\n")
    print("Final full trace fit:")
    print("Parameters: ", x0_f_params)
    # print('x0: ', dict(zip(x0_f_params, x0_f)))
    print("bounds: ", bounds_f)
    print("args: ", ext_args_f)
    print("------------")

    init_args.update(dict(zip(x0_f_params, best_fit_f)))
    final_args = init_args
    print("Error: ", fit_error_f["r2"])
    tnow = datetime.now()
    ret_val = {
        "time": str(tnow),
        "ion_channel": ion_channel,
        "graph": graph,
        "params": final_args,
        "bounds": list(bounds_f),
        "ext_args": ext_args_f,
        "error": fit_error_f,
        "equation": str_eq_f,
    }
    if fit_error_f["r2"] > R2LIMIT:
        # if plot or full_plot:
        print("best fit: ", final_args)
        print("Error: ", fit_error_f)
        print("Equation:\n", str_eq_f)
        v_full = np.arange(0.1, -0.11, -0.01)
        # t_full = np.arange(0, xf_interp[-1], 1e-4)
        trace_size_full = [0]
        x_full = []
        for vi in range(len(v_full)):
            trace_size_full.append(trace_size_full[vi] + len(time_interp))
            x_full.extend(time_interp)

        if not interp:
            ext_args_f["trace_size"] = list(trace_size_interp)
            ext_args_f["onset"] = list(onsets_interp)
            ext_args_f["offset"] = list(offsets_interp)
            ext_args_f["peak"] = list(peaks_interp)
            I_interp_f = full_trace(xf_interp, best_fit_f, ext_args_f, cur_type)[0]
            ext_args_f["volt"] = list(v_full)
            ext_args_f["trace_size"] = list(trace_size_full)
            I_interp_f2 = full_trace(x_full, best_fit_f, ext_args_f, cur_type)[0]
        else:
            I_interp_f = I_f
            I_interp_f2 = I_f

        if plot or full_plot or save_plot:
            fig = plt.figure()
            # fig, ((ax0, ax1), (ax2, ax3), (ax4, ax5)) = plt.subplots(nrows=3, ncols=2, figsize=(9, 9))
            # fig = plt.figure(figsize=(9,9))
            # gs = gridspec.GridSpec(nrows=3,
            #            ncols=2,
            #            figure=fig,
            #            width_ratios= [1, 1, 1],
            #            height_ratios=[1, 1, 1],
            #            wspace=0.3,
            #            hspace=0.3)
            # ax1 = fig.add_subplot(gs[0, 0])
            ax1 = fig.add_subplot()
            fig1 = plt.figure(1)
        for i in range(len(volt)):
            xx = np.asarray(xf[trace_size[i] : trace_size[i + 1]])
            yy = np.asarray(yf[trace_size[i] : trace_size[i + 1]])
            if not interp:
                x_fit = np.asarray(
                    xf_interp[trace_size_interp[i] : trace_size_interp[i + 1]]
                )
                y_fit = np.asarray(
                    I_interp_f[trace_size_interp[i] : trace_size_interp[i + 1]]
                )
            else:
                x_fit = np.asarray(xf[trace_size[i] : trace_size[i + 1]])
                y_fit = np.asarray(I_interp_f[trace_size[i] : trace_size[i + 1]])
            if plot or full_plot or save_plot:
                plt.plot(xx * 1e3, yy * 1e6, "o", label=volt[i] * 1e3)
                plt.plot(x_fit * 1e3, y_fit * 1e6)

        if not interp:
            ext_args_f["trace_size"] = list(trace_size)
            ext_args_f["onset"] = list(onsets)
            ext_args_f["offset"] = list(offsets)
            ext_args_f["peak"] = list(peaks)
            ext_args_f["volt"] = list(volt)

        if plot or full_plot or save_plot:
            ts = str(tnow.timestamp())

            title = "Raw and analized data for %s" % (
                data["ion_channel"]["channel_name"]
            )
            plt.title(title + "\nEquations: %s" % (str_eq_f), fontsize=8)
            plt.margins(x=0.1, y=0.1)
            plt.legend(
                bbox_to_anchor=(1.01, 0.25, 0.2, 0),
                loc=3,
                mode="expand",
                borderaxespad=0.0,
                fontsize=8,
            )
            if save_plot:
                plt.savefig(
                    PATH2
                    + str(graph["id"])
                    + "_"
                    + ion_channel["channel_name"]
                    + "_sample_"
                    + ts
                    + ".png",
                    bbox_inches="tight",
                    format="png",
                )
                plt.close()
                # save_plotly(fig1, fig1_file)
                # ply = tls.mpl_to_plotly(plt_data[n])
                # fig = go.Figure(data=ply_data)
                # plot(fig, filename='plot_'+str(n)+'.html',auto_open=False)

            v_f = np.arange(-0.1, 0.1, 0.001)
            fig2 = plt.figure(2)
            ax2 = fig.add_subplot()
            # peaks_f = []
            for i in range(len(v_full)):
                if not interp:
                    y_full = np.asarray(
                        I_interp_f2[trace_size_full[i] : trace_size_full[i + 1]]
                    )
                else:
                    y_full = np.asarray(I_interp_f2[trace_size[i] : trace_size[i + 1]])
                # peaks_f.append(max(abs(y_full)))
                plt.plot(time_interp * 1e3, y_full * 1e6, label=round(v_full[i] * 1e3))
            plt.title("Current/Time")
            plt.xlabel("Time(ms)")
            plt.ylabel("Current(uA)")
            plt.margins(x=0.1, y=0.1)
            plt.legend(
                bbox_to_anchor=(1.01, 0.01, 0.2, 0),
                loc=3,
                mode="expand",
                borderaxespad=0.0,
                fontsize=8,
            )
            if save_plot:
                plt.savefig(
                    PATH2
                    + str(graph["id"])
                    + "_"
                    + ion_channel["channel_name"]
                    + "_it_"
                    + ts
                    + ".png",
                    bbox_inches="tight",
                    format="png",
                )
                plt.close()

            fig3 = plt.figure(3)
            ax3 = fig.add_subplot()
            if "v_half" in ext_args_f and len(best_fit_f) < 7:
                v_half = ext_args_f["v_half"]
            else:
                v_half = best_fit_f[9]
            if "slope" in ext_args_f and len(best_fit_f) < 7:
                slope = ext_args_f["slope"]
            else:
                slope = best_fit_f[10]
            if "g" in ext_args_f and len(best_fit_f) < 7:
                g = ext_args_f["g"]
            else:
                g = best_fit_f[12]
            if "e_rev" in ext_args_f and len(best_fit_f) < 7:
                e_rev = ext_args_f["e_rev"]
            else:
                e_rev = best_fit_f[11]

            I_inf = 1 / (1 + np.exp(-(v_f - v_half) / slope))
            # I_peak = g * I_inf * (v_f-e_rev)
            I_peak = boltzmann(v_f, estimates_peak[0], estimates_peak[1], estimates_peak[3], estimates_peak[2], kwargs=kwargs)
            plt.plot(x_fit_peak, y_fit_peak, label="I_fit")
            plt.plot(volt, peak_yy, "--o", label="I_Peak_sample")
            # plt.plot(v_full, peaks_f, label="I_Peak")
            plt.plot(v_f, I_peak, label="I_Peak")
            plt.title("I/v curve")
            plt.margins(x=0.1, y=0.1)
            plt.legend()
            if save_plot:
                plt.savefig(
                    PATH2
                    + str(graph["id"])
                    + "_"
                    + ion_channel["channel_name"]
                    + "_iv_"
                    + ts
                    + ".png",
                    bbox_inches="tight",
                    format="png",
                )
                plt.close()

            fig4 = plt.figure(4)
            ax4 = fig.add_subplot()
            alpha_f = alphabeta(v_f, args=[best_fit_f[0], best_fit_f[1]], kwargs={})
            beta_f = alphabeta(
                v_f,
                args=[best_fit_f[2], best_fit_f[3], best_fit_f[4]],
                kwargs={"type": "sigmoid"},
            )
            c_f = [best_fit_f[5]] * len(v_f)
            plt.plot(v_f * 1e3, alpha_f, label="alpha")
            plt.plot(v_f * 1e3, beta_f, label="beta")
            plt.plot(v_f * 1e3, c_f, label="c")
            plt.title("alpha-beta/v")
            plt.margins(x=0.1, y=0.1)
            plt.legend()
            if save_plot:
                plt.savefig(
                    PATH2
                    + str(graph["id"])
                    + "_"
                    + ion_channel["channel_name"]
                    + "_ab_"
                    + ts
                    + ".png",
                    bbox_inches="tight",
                    format="png",
                )
                plt.close()

            fig5 = plt.figure(5)
            ax5 = fig.add_subplot()
            tau_f = 1 / (alpha_f + beta_f)
            tau_f2 = 1 / (alpha_f + c_f)
            ninf_f = alpha_f / (alpha_f + beta_f)
            ninf_f2 = alpha_f / (alpha_f + c_f)
            plt.plot(v_f * 1e3, tau_f * 1e3, label="tau")
            plt.plot(v_f * 1e3, tau_f2 * 1e3, label="tau2")
            plt.title("tau/v")
            plt.margins(x=0.1, y=0.1)
            plt.legend()
            if save_plot:
                plt.savefig(
                    PATH2
                    + str(graph["id"])
                    + "_"
                    + ion_channel["channel_name"]
                    + "_tau_"
                    + ts
                    + ".png",
                    bbox_inches="tight",
                    format="png",
                )
                plt.close()

            fig6 = plt.figure(6)
            ax6 = fig.add_subplot()

            plt.plot(v_f * 1e3, I_inf, label="peak_inf")
            plt.plot(v_f * 1e3, ninf_f, label="q1_inf")
            plt.plot(v_f * 1e3, ninf_f2, label="q2_inf")
            plt.plot(v_f * 1e3, (ninf_f * ninf_f2) ** (k - n), label="q_inf")
            plt.plot(v_f * 1e3, (1 - ninf_f * ninf_f2) ** n, label="p_inf^n")
            n = ext_args_f["npower"]
            k = ext_args_f["kpower"]
            pqco = np.math.factorial(int(k)) / (
                np.math.factorial(int(k) - int(n)) * np.math.factorial(int(n))
            )
            p_inf_f = pqco * (1 - ninf_f * ninf_f2) ** n * (ninf_f * ninf_f2) ** (k - n)
            plt.plot(v_f * 1e3, p_inf_f, label="pq_inf")
            plt.plot(
                v_f * 1e3, (1 - ninf_f * ninf_f2) ** n * I_inf, label="p_inf*peak_inf"
            )
            plt.plot(v_f * 1e3, p_inf_f * I_inf, label="pq_inf*peak_inf")
            plt.title("n_inf")
            plt.margins(x=0.1, y=0.1)
            plt.legend(
                bbox_to_anchor=(1.01, 0.25, 0.2, 0),
                loc=3,
                mode="expand",
                borderaxespad=0.0,
                fontsize=8,
            )

            if save_plot:
                plt.savefig(
                    PATH2
                    + str(graph["id"])
                    + "_"
                    + ion_channel["channel_name"]
                    + "_ninf_"
                    + ts
                    + ".png",
                    bbox_inches="tight",
                    format="png",
                )
                plt.close()
            else:
                plt.show()

            return ret_val, [fig1,fig2,fig3,fig4,fig5,fig6]
            # return ret_val, [ax1, ax2, ax3, ax4, ax5, ax6]

    return ret_val
