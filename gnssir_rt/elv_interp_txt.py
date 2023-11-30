"""
written by David Purnell
https://github.com/purnelldj
"""

import numpy as np
from scipy import interpolate
from scipy.optimize import least_squares


def elv_fit(t, elv, kspac, printrmse=False):
    # first just looking for any time where there is a gap in data
    breaks = [
        ind for ind in range(1, len(t)) if t[ind] - t[ind - 1] > kspac / 2
    ]
    breaks = np.append(0, breaks)
    breaks = np.append(breaks, len(t))
    breaks = np.array(breaks, dtype=int)
    tfit = []
    elvfit = []
    if printrmse:
        elvcompare = []
    # then looping through section by section in between the gaps
    for i in range(1, len(breaks)):
        tt = np.array(t[breaks[i - 1] : breaks[i]], dtype=float)
        if len(tt) < 2:
            continue
        elvt = np.array(elv[breaks[i - 1] : breaks[i]], dtype=float)
        knot_s = tt[0] - np.mod(tt[0], kspac)
        knot_e = tt[-1] - np.mod(tt[-1], kspac) + kspac
        if knot_e < knot_s + 3 * kspac:
            continue
        knots = np.linspace(knot_s, knot_e, int((knot_e - knot_s) / kspac + 1))
        nearestf = interpolate.interp1d(tt, elvt, kind="nearest")
        cp_in = nearestf(knots[1:-1])
        cp_in = np.append(elvt[0], cp_in)
        cp_in = np.append(cp_in, elvt[-1])
        elvcompare_t = elvt
        tfit_t = tt

        def resid_spline(cp):
            ffit = interpolate.interp1d(knots, cp, kind="cubic")
            resid = elvcompare_t - ffit(tfit_t)
            return resid

        cp_out = least_squares(resid_spline, cp_in, method="trf")
        ffit_t = interpolate.interp1d(knots, cp_out.x, kind="cubic")
        elvfit_t = ffit_t(tfit_t)
        tfit_t = tfit_t[elvfit_t > 0]
        elvfit_t = elvfit_t[elvfit_t > 0]
        tfit = np.append(tfit, tfit_t)
        elvfit = np.append(elvfit, elvfit_t)
        if printrmse:
            elvcompare = np.append(elvcompare, elvcompare_t)
    elvfit = np.array(elvfit, dtype=float)
    if printrmse:
        elvcompare = np.array(elvcompare, dtype=float)
        rms_fit = np.sqrt(np.mean((elvfit - elvcompare) ** 2))
        print("RMSE of fit and input elevation is = " + str(rms_fit))
    return tfit, elvfit


def elv_interp_array(snrdata, kspac):
    snrdata_interp = np.empty((0, 5))
    for sat in np.unique(snrdata[:, 0]):
        # print('sat now is ' + str(int(sat)) + ' ******')
        tfilter = snrdata[:, 0] == sat
        tt = snrdata[tfilter]
        elv_in = tt[:, 1]
        t_in = tt[:, 3]
        if len(np.unique(elv_in)) < 3:
            continue
        tfit, elvfit = elv_fit(t_in, elv_in, kspac)
        ttt = tt[np.in1d(t_in, tfit)]
        ttt[:, 1] = elvfit
        snrdata_interp = np.vstack((snrdata_interp, ttt))
    return snrdata_interp
