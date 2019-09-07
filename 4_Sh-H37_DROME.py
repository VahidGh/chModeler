"""
fit model to Sh-H37_DROME
"""

import os
import argparse
import numpy as np
import pickle
from datetime import datetime
from matplotlib import pyplot as plt
from fetcher import get_data_via_api
from fitter import *



def fit_vclamp(dataset, args=None):

    print("-----------------------------------------------------------------")

    data = dataset[0]
    if "ion_channel" in data:
        ion_channel = data["ion_channel"]
        print('Ion channel: ',ion_channel)
        # print("Animal: {}".format(ion_channel["animal"]))
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

    cur_type = 1
    if current_type:
        cur_type = current_type
        print("User Current type: ", current_type)
    else:
        print("Current type: ", cur_type)
    
    volt = []
    time_trace = []
    series_trace = []
    series = []
    for graph in dataset:
        # print(len(graph))
        # print(graph[0].keys())
        # print(graph[0]['x_var'])
        # print(graph[0]['y_var'])
        print([graph[0]['x_var']['type'],graph[0]['y_var']['type']])
        # volt = data.prot.volt
        # time = data.prot.time
        # current = data.prot.current


    print("onset id: ", onset_id)
    offsets = []
    onsets = []
    peak_yy = []
    tail_yy = []

    n_trace, n_pt_total = volt.shape
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
    bounds_boltz = [[-0.2, -0.1, -0.2, 0], [0.2, 0.1, 0.2, np.inf]]
    x_scale_inf = [1e-2, 1e-3, 1e-2, 1e-8]
    kwargs.update({"x_scale": x_scale_inf})

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
        x0_peak = [volt[vmedi], 0.009, e_rev0, g0]
        print(
            "max , min peaks: ",
            (max_peaksi, peak_yy[max_peaksi]),
            (min_peaksi, peak_yy[min_peaksi]),
        )
        print("x0 peak: ", x0_peak)
        estimates_peak, x_fit_peak, y_fit_peak, boltzmann_error = boltzmannfit(
            volt, peak_yy, x0=x0_peak, plot=plot, kwargs=kwargs, bounds=bounds_boltz
        )
    else:
        if isinstance(args["init_params"], dict):
            ninit = 0
        else:
            ninit = len(args["init_params"])
            print("\nNumber of initial params: ", ninit)
        if ninit > 0:
            r2s = np.zeros(ninit)
            boltzmann_errors = [[]] * ninit
            x0_peaks = [[]] * ninit
            estimates_peaks = [[]] * ninit
            x_fit_peaks = [[]] * ninit
            y_fit_peaks = [[]] * ninit
            for i, init_param in enumerate(args["init_params"]):
                x0_peaks[i] = list(init_param["params"].values())[:4]
                estimates_peaks[i], x_fit_peaks[i], y_fit_peaks[i], boltzmann_errors[
                    i
                ] = boltzmannfit(
                    volt,
                    peak_yy,
                    x0=x0_peaks[i],
                    plot=plot,
                    kwargs=kwargs,
                    bounds=bounds_boltz,
                )
                r2s[i] = boltzmann_errors[i]["r2"]
            maxerri = np.argmax(np.asarray(r2s))
            estimates_peak = estimates_peaks[maxerri]
            x_fit_peak = x_fit_peaks[maxerri]
            y_fit_peak = y_fit_peaks[maxerri]
            boltzmann_error = boltzmann_errors[maxerri]
            x0_peak = x0_peaks[maxerri]
            print("\ntop x0 peak, index: ", maxerri, x0_peak)

        else:
            x0_peak = list(args["init_params"]["params"].values())[:4]
            print("x0 peak: ", x0_peak)
            estimates_peak, x_fit_peak, y_fit_peak, boltzmann_error = boltzmannfit(
                volt, peak_yy, x0=x0_peak, plot=plot, kwargs=kwargs, bounds=bounds_boltz
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

    if "init_params" not in args:
        x_scale = [1, 1, 1, 1e-3, 1, 1]
        x00 = 0.01
        r0 = 10
        mins = 1e-3
        maxv = 0.3
        maxr = 500
        params = ["ra"]
        x0 = [r0]
        min_params = [mins]
        max_params = [maxr]
        params.append("sa")
        x0.append(-r0)
        min_params.append(-maxr)
        max_params.append(maxr)
        params.append("rb")
        x0.append(10 * r0)
        min_params.append(mins)
        max_params.append(maxr)
        params.append("vb")
        x0.append(-x00)
        min_params.append(-maxv)
        max_params.append(maxv)
        params.append("sb")
        x0.append(x00)
        min_params.append(mins)
        max_params.append(maxv)
        params.append("rc")
        x0.append(x00)
        min_params.append(0)
        max_params.append(maxr)

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
            print("Full trace params, x0, bounds: ", params, x0, bounds)

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
            if i == 0:
                x0_1 = copy.deepcopy(x0)

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
            ext_args.update({"x_scale": x_scale})

            if plot:
                print("x0: ", dict(zip(params, x0)))
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
        x0_2 = copy.deepcopy(x0_1)
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
                print("x0: ", dict(zip(params_2, x0_2)))
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
            ext_args_2["x_scale"] = x_scale

            if plot:
                print("bounds: ", bounds_2)
                print("args: ", ext_args_2)

            I2, str_eq_2, best_fit_2, fit_error_2 = fit_full_trace(
                x2, y2, x0_2, ext_args_2, cur_type, bounds=bounds_2
            )

            x0_2 = best_fit_2
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

        tau_args_exp = {"type": "exp"}
        atau_params = params[:2]
        x0_atau = [ra_v[1], sa_v[1]]
        min_atau = min_params[:2]
        max_atau = max_params[:2]
        best_fit_atau, x_atau, y_atau = abfit(
            volt,
            a_v,
            x0_atau,
            bounds=(min_atau, max_atau),
            ext_args=tau_args_exp,
            plot=plot,
        )

        tau_args_sig = {"type": "sigmoid"}
        btau_params = params[2:5]
        x0_btau = [rb_v[1], vb_v[1], sb_v[1]]
        min_btau = min_params[2:5]
        max_btau = max_params[2:5]
        best_fit_btau, x_btau, y_btau = abfit(
            volt,
            b_v,
            x0_btau,
            bounds=(min_btau, max_btau),
            ext_args=tau_args_sig,
            plot=plot,
        )

        ctau_params = [params[5]]
        min_ctau = [min_params[5]]
        max_ctau = [max_params[5]]
        x0_ctau = [1e-1]
        y_ctau = np.asarray([1e-1] * len(y_atau))
        best_fit_ctau = [1e-1]

        if plot:
            print("\n")
            print("alpha initial parameters: ", dict(zip(atau_params, x0_atau)))
            print("alpha parameters: ", dict(zip(atau_params, best_fit_atau)))
            print("beta initial parameters: ", dict(zip(btau_params, x0_btau)))
            print("beta parameters: ", dict(zip(btau_params, best_fit_btau)))
            print("c initial parameters: ", dict(zip(ctau_params, x0_ctau)))
            print("c parameters: ", dict(zip(ctau_params, best_fit_ctau)))

            plt.plot(x_atau, 1 / (y_atau + y_btau), label="ab_tau_fit")
            plt.plot(x_atau, 1 / (y_atau + y_ctau), label="ac_tau_fit")
            plt.plot(volt, act_tau, "--o", label="tau")
            plt.plot(volt, act_tau2, "--o", label="tau2")
            plt.legend()
            plt.show()

        x0_f_params = []
        x0_f = []
        min_bound_f = []
        max_bound_f = []

        x0_f_params.extend(atau_params)
        x0_f.extend(best_fit_atau)
        min_bound_f.extend(min_atau)
        max_bound_f.extend(max_atau)

        x0_f_params.extend(btau_params)
        x0_f.extend(best_fit_btau)
        min_bound_f.extend(min_btau)
        max_bound_f.extend(max_btau)

        x0_f_params.extend(ctau_params)
        x0_f.extend(best_fit_ctau)
        min_bound_f.extend(min_ctau)
        max_bound_f.extend(max_ctau)

        x_scale_f = x_scale
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
        ext_args_f["x_scale"] = list(x_scale_f)

    if "init_params" in args:
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
        else:
            mins = 1e-3
            maxv = 0.3
            maxr = 500

            min_params = [mins]
            max_params = [maxr]
            min_params.append(-maxr)
            max_params.append(maxr)
            min_params.append(mins)
            max_params.append(maxr)
            min_params.append(-maxv)
            max_params.append(maxv)
            min_params.append(mins)
            max_params.append(maxv)
            min_params.append(0)
            max_params.append(maxr)
            bounds_f = (min_params, max_params)

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
        x_scale_f = [1, 1, 1, 1e-3, 1, 1]
        ext_args_f["x_scale"] = list(x_scale_f)
        ext_args_f.update(ext_args)

        if ninit > 0:
            r2_fs = np.zeros(ninit)
            fit_error_fs = [[]] * ninit
            x0_fs = [[]] * ninit
            x0_f_paramss = [[]] * ninit
            best_fit_fs = [[]] * ninit
            str_eq_fs = [[]] * ninit
            I_fs = [[]] * ninit
            for i, init_param in enumerate(args["init_params"]):
                if "bounds" in init_param:
                    bounds_f = init_param["bounds"]
                x0_fs[i] = list(init_param["params"].values())[4:]
                x0_f_paramss[i] = list(init_param["params"].keys())[4:]
                I_fs[i], str_eq_fs[i], best_fit_fs[i], fit_error_fs[i] = fit_full_trace(
                    xf, yf, x0_fs[i], ext_args_f, cur_type, bounds=bounds_f
                )
                r2_fs[i] = fit_error_fs[i]["r2"]
            maxerrfi = np.argmax(np.asarray(r2_fs))
            best_fit_f = best_fit_fs[maxerrfi]
            str_eq_f = str_eq_fs[maxerrfi]
            I_f = I_fs[maxerrfi]
            fit_error_f = fit_error_fs[maxerrfi]
            x0_f = x0_fs[maxerrfi]
            x0_f_params = x0_f_paramss[maxerrfi]
            if "bounds" in args["init_params"][maxerrfi]:
                bounds_f = args["init_params"][maxerrfi]["bounds"]
            print("\ntop x0 final, index: ", maxerrfi, x0_f)

        else:
            x0_f = list(args["init_params"]["params"].values())[4:]
            x0_f_params = list(args["init_params"]["params"].keys())[4:]
            I_f, str_eq_f, best_fit_f, fit_error_f = fit_full_trace(
                xf, yf, x0_f, ext_args_f, cur_type, bounds=bounds_f
            )
    else:
        I_f, str_eq_f, best_fit_f, fit_error_f = fit_full_trace(
            xf, yf, x0_f, ext_args_f, cur_type, bounds=bounds_f
        )

    # if plot or full_plot:
    print("\n")
    print("Final full trace fit:")
    print("Parameters: ", x0_f_params)
    print("x0: ", dict(zip(x0_f_params, x0_f)))
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
        "x0": list(x0_f),
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

        if plot or full_plot:
            plt.figure()
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
            if plot or full_plot:
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

            v_f = np.arange(-0.1, 0.1, 0.001)
            plt.figure()
            for i in range(len(v_full)):
                if not interp:
                    y_full = np.asarray(
                        I_interp_f2[trace_size_full[i] : trace_size_full[i + 1]]
                    )
                else:
                    y_full = np.asarray(I_interp_f2[trace_size[i] : trace_size[i + 1]])
                plt.plot(time_interp * 1e3, y_full * 1e6, label=round(v_full[i] * 1e3))
            plt.title("Current/Time")
            plt.xlabel("Time(ms)")
            plt.ylabel("Current(uA)")
            plt.margins(x=0.1, y=0.1)
            plt.legend(
                bbox_to_anchor=(1.01, 0.1, 0.2, 0),
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

            plt.figure()
            plt.plot(x_fit_peak, y_fit_peak, label="I_fit")
            plt.plot(volt, peak_yy, "--o", label="I_Peak")
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

            plt.figure()
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

            plt.figure()
            if "v_half" in ext_args_f and len(best_fit_f) < 7:
                v_half = ext_args_f["v_half"]
            else:
                v_half = best_fit_f[9]
            if "slope" in ext_args_f and len(best_fit_f) < 7:
                slope = ext_args_f["slope"]
            else:
                slope = best_fit_f[10]

            I_inf = 1 / (1 + np.exp(-(v_f - v_half) / slope))

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

    return ret_val


def main():

    parser = argparse.ArgumentParser(description="Fit model to Sh-H37_DROME.")
    parser.add_argument(
        "-dp", "--path", default="data/Sh-H37_DROME/", help="Directory containing plot data.",
    )
    parser.add_argument(
        "-p", "--plot", default=False, help="If TRUE, all step plots will be shown."
    )
    parser.add_argument(
        "-fp", "--final_plot", default=False, help="If TRUE, final plots will be shown."
    )
    parser.add_argument(
        "-s", "--save", default=True, help="If TRUE, results will be saved."
    )
    parser.add_argument(
        "-sp", "--save_path", default="data/", help="Path for saving final results."
    )
    args = parser.parse_args()
    print(args)

    fit_args = {}
    path = args.path
    if args.save is True:
        save_path = args.save_path + str(datetime.now().strftime('%Y%m%d%H%M%S')) + '/'
        # print(save_path)
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        fit_args['save_path'] = save_path
    fit_args['plot'] = args.plot
    fit_args['final_plot'] = args.final_plot

    if not os.path.exists(path + 'dataset.pickle'):
        print('path: ', path)
        graph_ids = [11, 12, 13, 14, 15, 16]
        # dataset = get_data_via_api(graph_ids, plot=True)
        data = []
        for graph_id in graph_ids:
            data.append(get_data_via_api(graph_id, plot=True))
        if not os.path.exists(path):
            print('Creating path: ', path)
            os.makedirs(path)
        pickle.dump(data, open(path + 'dataset.pickle', 'wb'))
        print('Dataset dumped at: ', path + 'dataset.pickle')
    else:
        print('Loading Dataset from: ', path + 'dataset.pickle')
        data = pickle.load(open(path + 'dataset.pickle', 'rb'), encoding='latin1')

    fit_args['graph_types'] = [['T', 'I'],
                                ['V', 'I_peak'],
                                ['V', 'G/G_max'],
                                ['V', 'Tau_a'],
                                ['V', 'Tau_i'],
                                ['V', 'I_norm_i']]

    fit_res = fit_vclamp(data, fit_args)
    

if __name__ == "__main__":
    main()
