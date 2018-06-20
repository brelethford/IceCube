# -*-coding:utf8-*-

# python
import os
import re

# SciPy
import numpy as np
from numpy.lib.recfunctions import drop_fields

from mstacking.psLLH import StackingPointSourceLLH, StackingMultiPointSourceLLH
from mstacking.ps_model import EnergyLLH#, EnergyDistLLH #,EnergyBDTLLH 

path = "/data/user/coenders/data/MultiYearPointSource/npz/"
#path = "/Users/mhuber/StackingSkylab/MultiYearPointSourceData"

#NOTE: this was used by mhuber for his WHSP sensitivity calculation. Unclear what else it was used for.

hem = np.sin(np.radians(-5.))

def ic40(**kwargs):
    livetime = 375.539
    print("\tLoading IC40 data...")
    exp = np.load(os.path.join(path, "IC40_exp.npy"))
    mc = np.load(os.path.join(path, "IC40_MC.npy"))
    
    #dec_bins = np.linspace(-1., 1., 25 + 1)

    
    dec_bins = np.unique(np.concatenate([
                           np.linspace(-1., -0.25, 10 + 1),
                           np.linspace(-0.25, 0.0, 5 + 1),
                           np.linspace(0.0, 1., 10 + 1),
                           ]))
    
    dec_bins_logE = np.unique(np.concatenate([
                           np.linspace(-1., -0.25, 10 + 1),
                           np.linspace(-0.25, 0.0, 10 + 1),
                           np.linspace(0.0, 1., 10 + 1),
                           ]))
    
    energy_bins = [np.linspace(1., 10., 35 + 1), dec_bins_logE]

    mc = mc[mc["logE"] > 1.]

    llh_model = EnergyLLH(twodim_bins=energy_bins,
                          sinDec_bins=dec_bins)

    if "upscale" in kwargs and kwargs["upscale"] is not None and (kwargs["upscale"] or not type(kwargs["upscale"]) == bool):
        lt = kwargs.pop("livetime", livetime)
        kwargs["upscale"] = (int(kwargs.pop("upscale")), lt)

    kwargs.pop("livetime", None)

    kwargs.setdefault("seed", 2008)

    llh = StackingPointSourceLLH(exp, mc, livetime,
                         llh_model=llh_model,
                         **kwargs
                         )


    print("{0:>80s}".format("[done]"))
    return llh




def ic59(**kwargs):
    livetime = 348.138
    print("\tLoading IC59 data...")
    exp = np.load(os.path.join(path, "IC59_exp.npy"))
    mc = np.load(os.path.join(path, "IC59_corrected_MC.npy"))

    dec_bins = np.unique(np.concatenate([
                         np.linspace(-1., -0.95, 2 + 1),
                         np.linspace(-0.95, -0.25, 25 + 1),
                         np.linspace(-0.25, 0.05, 15 + 1),
                         np.linspace(0.05, 1., 10 + 1),
                         ]))
    dec_bins_logE = np.unique(np.concatenate([
                         np.linspace(-1., 0.05, 20 + 1),
                         np.linspace(0.05, 1., 10 + 1),
                         ]))

    energy_bins = [np.linspace(2., 9.5, 50 + 1), dec_bins]

    mc = mc[mc["logE"] > 1.]

    llh_model = EnergyLLH(twodim_bins=energy_bins,
                          sinDec_bins=dec_bins)

    if "upscale" in kwargs and kwargs["upscale"] is not None and (kwargs["upscale"] or not type(kwargs["upscale"]) == bool):
        lt = kwargs.pop("livetime", livetime)
        kwargs["upscale"] = (int(kwargs.pop("upscale")), lt)

    kwargs.pop("livetime", None)

    kwargs.setdefault("seed", 2009)

    llh = StackingPointSourceLLH(exp, mc, livetime,
                         llh_model=llh_model,
                         **kwargs
                         )

    print("{0:>80s}".format("[done]"))
    return llh


def ic79(**kwargs):
    livetime = 315.506
    print("\tLoading IC79 data...")

    decorr = kwargs.pop("no_mese", False)

    if decorr:
        print("\t\tRemove MESE from exp and MC")
        exp = np.load(os.path.join(path, "IC79_noMESE_exp.npy"))
        mc = np.load(os.path.join(path, "IC79_noMESE_corrected_MC.npy"))
    else:
        exp = np.load(os.path.join(path, "IC79_exp.npy"))
        mc = np.load(os.path.join(path, "IC79_corrected_MC.npy"))

    dec_bins = np.unique(np.concatenate([
                         np.linspace(-1., -0.75, 10 + 1),
                         np.linspace(-0.75, 0., 15 + 1),
                         np.linspace(0., 1., 20 + 1)
                         ]))

    energy_bins = [np.linspace(2., 9., 40 + 1), dec_bins]

    mc = mc[mc["logE"] > 1.]

    llh_model = EnergyLLH(twodim_bins=energy_bins,
                          sinDec_bins=dec_bins)

    if "upscale" in kwargs and kwargs["upscale"] is not None and (kwargs["upscale"] or not type(kwargs["upscale"]) == bool):
        lt = kwargs.pop("livetime", livetime)
        kwargs["upscale"] = (int(kwargs.pop("upscale")), lt)

    kwargs.pop("livetime", None)

    kwargs.setdefault("seed", 2010)

    llh = StackingPointSourceLLH(exp, mc, livetime,
                         llh_model=llh_model,
                         **kwargs
                         )

    print("{0:>80s}".format("[done]"))
    return llh


def ic79b(**kwargs):
    livetime = 315.506
    print("\tLoading IC79 spline data...")

    decorr = kwargs.pop("no_mese", False)

    exp = np.load(os.path.join(path, "IC79b_exp.npy"))
    if decorr:
        print("\t\tRemove MESE events from MC, none in exp anyways")
        mc = np.load(os.path.join(path, "IC79b_noMESE_corrected_MC.npy"))
    else:
        mc = np.load(os.path.join(path, "IC79b_corrected_MC.npy"))

    dec_bins = np.unique(np.concatenate([
                         np.linspace(-1., 1., 50),
                         ]))

    energy_bins = [np.linspace(2., 9., 40 + 1), dec_bins]

    mc = mc[mc["logE"] > 1.]

    llh_model = EnergyLLH(twodim_bins=energy_bins,
                          sinDec_bins=dec_bins)

    if "upscale" in kwargs and kwargs["upscale"] is not None and (kwargs["upscale"] or not type(kwargs["upscale"]) == bool):
        lt = kwargs.pop("livetime", livetime)
        kwargs["upscale"] = (int(kwargs.pop("upscale")), lt)

    kwargs.pop("livetime", None)

    kwargs.setdefault("seed", 2010)

    llh = StackingPointSourceLLH(exp, mc, livetime,
                         llh_model=llh_model,
                         **kwargs
                         )

    print("{0:>80s}".format("[done]"))
    return llh


def ic86_I(**kwargs):
    livetime = 332.61
    print("\tLoading IC86 data...")

    decorr = kwargs.pop("no_mese", False)

    if decorr:
        print("\t\tRemove MESE from exp and MC")
        exp = np.load(os.path.join(path, "IC86_noMESE_exp.npy"))
        mc = np.load(os.path.join(path, "IC86_noMESE_corrected_MC.npy"))
    else:
        exp = np.load(os.path.join(path, "IC86_exp.npy"))
        mc = np.load(os.path.join(path, "IC86_corrected_MC.npy"))

    sinDec = kwargs.pop("sinDec", [-1., 1.])
    exp = exp[(exp["sinDec"] > sinDec[0]) & (exp["sinDec"] < sinDec[-1])]
    mc = mc[(mc["sinDec"] > sinDec[0]) & (mc["sinDec"] < sinDec[-1])]

    dec_bins = np.unique(np.concatenate([
                         np.linspace(-1., -0.2, 10 + 1),
                         np.linspace(-0.2, hem, 4 + 1),
                         np.linspace(hem, 0.2, 5 + 1),
                         np.linspace(0.2, 1., 10),
                         ]))

    dec_bins = dec_bins[(dec_bins >= sinDec[0]) & (dec_bins <= sinDec[1])]
    dec_bins = np.unique(np.concatenate([sinDec, dec_bins]))

    energy_bins = [np.linspace(1., 10., 40 + 1), dec_bins]

    mc = mc[mc["logE"] > 1.]

    llh_model = EnergyLLH(twodim_bins=energy_bins,
                          sinDec_bins=dec_bins)

    if "upscale" in kwargs and kwargs["upscale"] is not None and (kwargs["upscale"] or not type(kwargs["upscale"]) == bool):
        lt = kwargs.pop("livetime", livetime)
        kwargs["upscale"] = (int(kwargs.pop("upscale")), lt)

    kwargs.pop("livetime", None)

    kwargs.setdefault("seed", 2011)

    llh = StackingPointSourceLLH(exp, mc, livetime,
                         llh_model=llh_model,
                         **kwargs
                         )

    print("{0:>80s}".format("[done]"))
    return llh


def ic86_2012_precut(**kwargs):
    livetime = 331.88
    print("\tLoading IC86-II with precut applied...")
    exp = np.load(os.path.join(path, "IC86-2012-precut_exp.npy"))
    mc = np.load(os.path.join(path, "IC86-2012-precut_MC.npy"))

    sinDec = kwargs.pop("sinDec", [-1., 1.])
    exp = exp[(exp["sinDec"] > sinDec[0]) & (exp["sinDec"] <= sinDec[-1])]
    mc = mc[(mc["sinDec"] > sinDec[0]) & (mc["sinDec"] <= sinDec[-1])]

    if "perc" in kwargs:
        perc = kwargs.pop("perc")
        print("\tCutting on percentage {0:8.3%}".format(perc))

        nbins = 400
        h, bins = np.histogram(exp["sinDec"], bins=nbins, range=sinDec)
        NpB = perc * np.sum(h, dtype=np.float) / nbins
        eps = np.zeros_like(h, dtype=np.float)
        eps[h>0] = NpB / h[h>0]
        eps[eps > 1] = 1.
        xdec = np.array([])
        ybdt = np.array([])
        for eps_i, lbin, ubin in zip(eps, bins[:-1], bins[1:]):
            m = (exp["sinDec"] > lbin)&(exp["sinDec"] < ubin)
            if not np.any(m):
                print("No events in", lbin, ubin)
                continue
            bdt_cut = np.percentile(exp["BDT"][m], 100.*(1. - eps_i))
            np.set_printoptions(precision=4)
            xdec = np.append(xdec, (lbin + ubin)/2.)
            ybdt = np.append(ybdt, bdt_cut)

        pol = np.polyfit(xdec, ybdt, 6)

        N1 = len(exp)

        polexp = np.polyval(pol, exp["sinDec"])
        polmc = np.polyval(pol, mc["sinDec"])

        mexp = exp["BDT"] > polexp
        mmc = mc["BDT"] > polmc

        exp = exp[mexp]
        mc = mc[mmc]

        N2 = len(exp)

        print("\tCut: {0:7.2%}".format(float(N2)/N1))

    elif "bdt" in kwargs and "bdt2" in kwargs:
        comp = kwargs.pop("comp", np.logical_or)
        bdt = kwargs.pop("bdt")
        if not callable(bdt):
            print("\tCutting on bdt {0:5.3f}".format(bdt))
            bdt_val = bdt
            bdt = lambda cz: bdt_val * np.ones_like(cz)
        else:
            print("\tCutting bdt from {0:5.3f} to {1:5.3f} in dec range".format(*bdt(np.asarray(sinDec))))
        bdt2 = kwargs.pop("bdt2")
        if not callable(bdt2):
            print("\tCutting on bdt2 {0:5.3f}".format(bdt2))
            bdt_val2 = bdt2
            bdt2 = lambda cz: bdt_val2 * np.ones_like(cz)
        else:
            print("\tCutting bdt2 from {0:5.3f} to {1:5.3f} in dec range".format(*bdt2(np.asarray(sinDec))))
        exp = exp[comp(exp["BDT"] > bdt(exp["sinDec"]),
                       exp["BDT2"] > bdt2(exp["sinDec"]))]
        mc = mc[comp(mc["BDT"] > bdt(mc["sinDec"]),
                mc["BDT2"] > bdt2(mc["sinDec"]))]

    elif "bdt" in kwargs:
        bdt = kwargs.pop("bdt")
        if not callable(bdt):
            print("\tCutting on bdt {0:5.3f}".format(bdt))
            bdt_val = bdt
            bdt = lambda cz: bdt_val * np.ones_like(cz)
        else:
            print("\tCutting bdt from {0:5.3f} to {1:5.3f} in dec range".format(*bdt(np.asarray(sinDec))))
        exp = exp[exp["BDT"] > bdt(exp["sinDec"])]
        mc = mc[mc["BDT"] > bdt(mc["sinDec"])]
    elif "bdt2" in kwargs:
        bdt = kwargs.pop("bdt2")
        if not callable(bdt):
            print("\tCutting on bdt2 {0:5.3f}".format(bdt))
            bdt_val = bdt
            bdt = lambda cz: bdt_val * np.ones_like(cz)
        else:
            print("\tCutting bdt from {0:5.3f} to {1:5.3f} in dec range".format(*bdt(np.asarray(sinDec))))
        exp = exp[exp["BDT2"] > bdt(exp["sinDec"])]
        mc = mc[mc["BDT2"] > bdt(mc["sinDec"])]

    nmin = 10
    nbins = 40
    dec_bins = np.linspace(*sinDec, num=nbins)
    while np.any(np.histogram(exp["sinDec"], bins=dec_bins)[0] < nmin):
        print("Found close-to-empty bins, reduce binsize")
        h, b = np.histogram(exp["sinDec"], bins=dec_bins)
        b = np.array([b[1:], b[:-1]])
        print(b[:, h<nmin])
        if nbins < 2:
            raise ValueError("Need more than 2 bins!")
        nbins -= 1
        dec_bins = np.linspace(*sinDec, num=nbins)

    energy_bins = [np.linspace(0., 9.5, 30 + 1),
                   dec_bins]

    llh_model = EnergyLLH(twodim_bins=energy_bins,
                          sinDec_bins=dec_bins)

    if "upscale" in kwargs and kwargs["upscale"] is not None and (kwargs["upscale"] or not type(kwargs["upscale"]) == bool):
        lt = kwargs.pop("livetime", 332.61)#livetime)
        kwargs["upscale"] = (int(kwargs.pop("upscale")), lt)

    kwargs.pop("livetime", None)

    kwargs.setdefault("seed", 2012)

    llh = StackingPointSourceLLH(exp, mc, livetime,
                         llh_model=llh_model,
                         **kwargs
                         )

    if "pol" in locals():
        llh.pol = pol

    print("{0:>80s}".format("[done]"))
    return llh


def ic86_2012(dataset=[2012, 2013], **kwargs):
    lt12 = 331.358
    lt13 = 359.768
    lt14 = 367.209
    livetime = 0.
    exp = list()

    decorr = kwargs.pop("no_mese", False)

    if 2012 in dataset:
        livetime += lt12
        print("\tLoading IC86-2012 ...")
        if decorr:
            print("\t\tRemove MESE from IC86-2012 exp")
            exp.append(np.load(os.path.join(path, "IC86-2012_noMESE_exp.npy")))
        else:
            exp.append(np.load(os.path.join(path, "IC86-2012_exp.npy")))
    if 2013 in dataset:
        livetime += lt13
        print("\tLoading IC86-2013 ...")
        exp.append(np.load(os.path.join(path, "IC86-2013_exp.npy")))

    if 2014 in dataset:
        livetime += lt14
        print("\tLoading IC86-2014 ...")
        exp.append(np.load(os.path.join(path, "IC86-2014_exp.npy")))


    exp = np.concatenate(exp)
    if decorr:
        print("\t\tRemove MESE from IC86-2012 MC")
        mc = np.load(os.path.join(path, "IC86-2012_noMESE_MC.npy"))
    else:
        mc = np.load(os.path.join(path, "IC86-2012_MC.npy"))

    sinDec = kwargs.pop("sinDec", [-1., 1.])
    exp = exp[(exp["sinDec"] > sinDec[0]) & (exp["sinDec"] < sinDec[-1])]
    mc = mc[(mc["sinDec"] > sinDec[0]) & (mc["sinDec"] < sinDec[-1])]

    exp = drop_fields(exp, ["BDT"])
    mc = drop_fields(mc, ["BDT"])

    dec_bins = np.unique(np.concatenate([
                         np.linspace(-1., -0.93, 4 + 1),
                         np.linspace(-0.93, -0.3, 10 + 1),
                         np.linspace(-0.3, 0.05, 9 + 1),
                         np.linspace(0.05, 1., 18 + 1),
                         ]))
    dec_bins = dec_bins[(dec_bins >= sinDec[0]) & (dec_bins <= sinDec[1])]
    dec_bins = np.unique(np.concatenate([sinDec, dec_bins]))

    energy_bins = [np.linspace(1., 9.5, 50 + 1), dec_bins]


    llh_model = EnergyLLH(twodim_bins=energy_bins,
                          sinDec_bins=dec_bins)

    if "upscale" in kwargs and kwargs["upscale"] is not None and (kwargs["upscale"] or not type(kwargs["upscale"]) == bool):
        lt = kwargs.pop("livetime", livetime)
        kwargs["upscale"] = (int(kwargs.pop("upscale")), lt)

    kwargs.pop("livetime", None)

    kwargs.setdefault("seed", 2012)

    llh = StackingPointSourceLLH(exp, mc, livetime,
                         llh_model=llh_model,
                         **kwargs
                         )

    print("{0:>80s}".format("[done]"))
    return llh


def ic86_2012_bdt(**kwargs):
    livetime = 331.88
    print("\tLoading IC86-II...")
    exp = np.load(os.path.join(path, "IC86-2012_exp.npy"))
    mc = np.load(os.path.join(path, "IC86-2012_MC.npy"))

    sinDec = kwargs.pop("sinDec", [-1., 1.])
    exp = exp[(exp["sinDec"] > sinDec[0]) & (exp["sinDec"] < sinDec[-1])]
    mc = mc[(mc["sinDec"] > sinDec[0]) & (mc["sinDec"] < sinDec[-1])]

    exp = drop_fields(exp, ["BDT2", "perc"])
    mc = drop_fields(mc, ["BDT2", "perc"])

    dec_bins = np.unique(np.concatenate([
                         np.linspace(-1., -0.92, 5 + 1),
                         np.linspace(-0.92, -0.15, 10 + 1),
                         np.linspace(-0.15, 0.01, 10 + 1),
                         np.linspace(0.01, 1., 20 + 1),
                         ]))
    dec_bins = dec_bins[(dec_bins >= sinDec[0]) & (dec_bins <= sinDec[1])]
    dec_bins = np.unique(np.concatenate([sinDec, dec_bins]))

    X = np.concatenate([exp["BDT"], mc["BDT"]])
    bdt_bins = np.percentile(X, [0., 20., 40., 60., 80., 100.])

    energy_bdt_bins = [np.linspace(1., 10., 40 + 1),
                       np.concatenate([[bdt_bins[0] - (bdt_bins[1] - bdt_bins[0])],
                                       bdt_bins,
                                       [bdt_bins[-1] + bdt_bins[-1] - bdt_bins[-2]]]),
                       dec_bins]

    llh_model = EnergyBDTLLH(bins=energy_bdt_bins,
                             sinDec_bins=dec_bins)

    if "upscale" in kwargs and kwargs["upscale"] is not None and (kwargs["upscale"] or not type(kwargs["upscale"]) == bool):
        lt = kwargs.pop("livetime", livetime)
        kwargs["upscale"] = (int(kwargs.pop("upscale")), lt)

    kwargs.pop("livetime", None)

    kwargs.setdefault("seed", 2012)

    llh = StackingPointSourceLLH(exp, mc, livetime,
                         llh_model=llh_model,
                         **kwargs
                         )

    print("{0:>80s}".format("[done]"))
    return llh


def mese(**kwargs):
    livetime = 988.54
    print("\tLoading MESE...")
    exp = np.load(os.path.join(path, "MESE_exp.npy"))
    mc = np.load(os.path.join(path, "MESE_MC.npy"))

    dec_bins = np.unique(np.concatenate([
                         np.linspace(-1., -0.92, 3 + 1),
                         np.linspace(-0.92, hem, 10 + 1),
                         ]))

    energy_bins = [np.linspace(2., 8.5, 40 + 1),
                   np.linspace(-1., hem, 4 + 1),
                   ]

    mc = mc[mc["logE"] > 1.]

    llh_model = EnergyLLH(twodim_bins=energy_bins,
                          sinDec_bins=dec_bins)

    if "upscale" in kwargs and kwargs["upscale"] is not None and (kwargs["upscale"] or not type(kwargs["upscale"]) == bool):
        lt = kwargs.pop("livetime", livetime)
        kwargs["upscale"] = (int(kwargs.pop("upscale")), lt)

    kwargs.pop("livetime", None)

    kwargs.setdefault("mode", "all")

    kwargs.setdefault("seed", 20101112)

    llh = StackingPointSourceLLH(exp, mc, livetime,
                         llh_model=llh_model,
                         **kwargs
                         )

    print("{0:>80s}".format("[done]"))
    return llh

def mese_followup(**kwargs):
    livetime = 988.54
    livetime += 358.402
    livetime += 368.381
    print("\tLoading MESE with 2 follow-up years...")
    exp = np.append(np.load(os.path.join(path, "MESE_exp.npy")),
                         np.load(os.path.join(path, "MESE_followup_exp.npy")))

    

    mc = np.load(os.path.join(path, "MESE_MC.npy"))

    if "dist" in exp.dtype.names:
        exp = drop_fields(exp, ["dist"], usemask=False)
    if "dist" in mc.dtype.names:
        mc = drop_fields(mc, ["dist"], usemask=False)

    sinDec = kwargs.pop("sinDec", [-1., hem])
    exp = exp[(exp["sinDec"] > sinDec[0]) & (exp["sinDec"] < sinDec[-1])]
    mc = mc[(mc["sinDec"] > sinDec[0]) & (mc["sinDec"] < sinDec[-1])]

    dec_bins = np.unique(np.concatenate([np.linspace(-1., -0.93, 4 + 1),
                                        np.linspace(-0.93, hem, 12 + 1),
                                        ]))

    dec_bins = dec_bins[(dec_bins >= sinDec[0]) & (dec_bins <= sinDec[1])]
    dec_bins = np.unique(np.concatenate([sinDec, dec_bins]))
    dec_bins_logE = np.linspace(-1., hem, 4 + 1)
    dec_bins_logE = dec_bins_logE[(dec_bins_logE >= sinDec[0]) & (dec_bins_logE <= sinDec[1])]
    dec_bins_logE = np.unique(np.concatenate([sinDec, dec_bins_logE]))

    energy_bins = [np.linspace(2., 8.5, 67 + 1), np.linspace(-1., hem, 4 + 1),
                                  ]
    mc = mc[mc["logE"] > 1.]
    llh_model = EnergyLLH(twodim_bins=energy_bins, sinDec_bins=dec_bins)

    if "upscale" in kwargs and kwargs["upscale"] is not None and (kwargs["upscale"] or not type(kwargs["upscale"]) == bool):
        lt = kwargs.pop("livetime", livetime)
        kwargs["upscale"] = (int(kwargs.pop("upscale")), lt)

    kwargs.pop("livetime", None)
    kwargs.setdefault("mode", "all")
    kwargs.setdefault("seed", 20101112)
    llh = StackingPointSourceLLH(exp, mc, livetime,
                            llh_model=llh_model,
                            **kwargs)
    
    

    print("{0:>80s}".format("[done]"))
    return llh





def load_4yr(years="IC40",*args, **kwargs):
    print("Loading 4year data...")


    scramble = kwargs.pop("scramble", True)
    upscale = kwargs.pop("upscale", None)
    no_mese = kwargs.pop("no_mese", False)
    # load all years from IC40 to IC86.2011
    llh = StackingMultiPointSourceLLH(**kwargs)

    kwargs["scramble"] = scramble

    if upscale is not None:
        upscale *= 10
    #------------------------------------------------
    if years == "IC40":
        kwargs["upscale"] = upscale
        llh40 = ic40(**kwargs)
        llh.add_sample("IC40", llh40)
    #------------------------------------------------
    elif years == "IC59":
        kwargs["upscale"] = upscale
        llh59 = ic59(**kwargs)
        llh.add_sample("IC59", llh59)
    #------------------------------------------------
    elif years == "IC79":
        kwargs["upscale"] = upscale
        llh79 = ic79(no_mese=no_mese,**kwargs)
        llh.add_sample("IC79", llh79)
    #------------------------------------------------
    elif years == "IC86":
        kwargs["upscale"] = upscale
        llh86 = ic86_I(no_mese=no_mese,**kwargs)
        llh.add_sample("IC86", llh86)

    #------------------------------------------------
    elif years == "IC40+IC59":
        kwargs["upscale"] = upscale
        llh40 = ic40(**kwargs)
        llh.add_sample("IC40", llh40)

        kwargs["upscale"] = upscale + 1 if not upscale is None else None
        llh59 = ic59(**kwargs)
        llh.add_sample("IC59", llh59)

    #-------------------------------------------------
    elif years == "IC40+IC59+IC79":
        kwargs["upscale"] = upscale
        llh40 = ic40(**kwargs)
        llh.add_sample("IC40", llh40)

        kwargs["upscale"] = upscale + 1 if not upscale is None else None
        llh59 = ic59(**kwargs)
        llh.add_sample("IC59", llh59)

        kwargs["upscale"] = upscale + 2 if not upscale is None else None
        llh79 = ic79(no_mese=no_mese,**kwargs)
        llh.add_sample("IC79", llh79)
    
    #-------------------------------------------------
    elif years == "IC40+IC59+IC79+IC86":
        kwargs["upscale"] = upscale
        llh40 = ic40(**kwargs)
        llh.add_sample("IC40", llh40)

        kwargs["upscale"] = upscale + 1 if not upscale is None else None
        llh59 = ic59(**kwargs)
        llh.add_sample("IC59", llh59)

        kwargs["upscale"] = upscale + 2 if not upscale is None else None
        llh79 = ic79(no_mese=no_mese,**kwargs)
        llh.add_sample("IC79", llh79)

        kwargs["upscale"] = upscale + 3 if not upscale is None else None
        llh86 = ic86_I(no_mese=no_mese,**kwargs)
        llh.add_sample("IC86", llh86)

    



    print("{0:>80s}".format("[done]"))
    return llh

def load_3yr_ic79b(*args, **kwargs):
    print("Loading 3year data with splined IC79...")
    scramble = kwargs.pop("scramble", True)
    upscale = kwargs.pop("upscale", None)
    no_mese = kwargs.pop("no_mese", False)


    # load all years from IC40 to IC86.2011
    llh = StackingMultiPointSourceLLH(**kwargs)

    kwargs["scramble"] = scramble

    if upscale is not None:
        upscale *= 10

    kwargs["upscale"] = upscale
    llh40 = ic40(**kwargs)
    llh.add_sample("IC40", llh40)

    kwargs["upscale"] = upscale + 1 if not upscale is None else None
    llh59 = ic59(**kwargs)
    llh.add_sample("IC59", llh59)

    kwargs["upscale"] = upscale + 2 if not upscale is None else None
    llh79b = ic79b(no_mese=no_mese, **kwargs)
    llh.add_sample("IC79b", llh79b)

    print("{0:>80s}".format("[done]"))
    return llh

def load_4yr_ic79b(*args, **kwargs):
    print("Loading 4year data with splined IC79...")
    scramble = kwargs.pop("scramble", True)
    upscale = kwargs.pop("upscale", None)
    no_mese = kwargs.pop("no_mese", False)


    # load all years from IC40 to IC86.2011
    llh = StackingMultiPointSourceLLH(**kwargs)

    kwargs["scramble"] = scramble

    if upscale is not None:
        upscale *= 10

    kwargs["upscale"] = upscale
    llh40 = ic40(**kwargs)
    llh.add_sample("IC40", llh40)

    kwargs["upscale"] = upscale + 1 if not upscale is None else None
    llh59 = ic59(**kwargs)
    llh.add_sample("IC59", llh59)

    kwargs["upscale"] = upscale + 2 if not upscale is None else None
    llh79b = ic79b(no_mese=no_mese, **kwargs)
    llh.add_sample("IC79b", llh79b)


    kwargs["upscale"] = upscale + 3 if not upscale is None else None
    llh86 = ic86_I(no_mese=no_mese, **kwargs)
    llh.add_sample("IC86", llh86)
        

    print("{0:>80s}".format("[done]"))
    return llh


def load_4yr_mese(*args, **kwargs):
    print("Loading 4year plus MESE")

    llh = load_4yr(*args, no_mese=True, **kwargs)

    upscale = kwargs.pop("upscale", None)
    if upscale is not None:
        upscale *= 10

    kwargs["upscale"] = upscale + 4 if not upscale is None else None
    llhMESE = mese(*args, **kwargs)
    llh.add_sample("MESE", llhMESE)

    return llh


def load_4yr_ic79b_mese(*args, **kwargs):
    print("Loading 4year plus MESE")

    llh = load_4yr_ic79b(*args, no_mese=True, **kwargs)

    upscale = kwargs.pop("upscale", None)
    if upscale is not None:
        upscale *= 10

    kwargs["upscale"] = upscale + 4 if not upscale is None else None
    llhMESE = mese(*args, **kwargs)
    llh.add_sample("MESE", llhMESE)

    return llh


def load_6yr(*args, **kwargs):
    print("Loading 6year")
    # load all years from IC40 to IC86.2011
    llh = load_4yr_ic79b(*args, **kwargs)

    upscale = kwargs.pop("upscale", None)
    if upscale is not None:
        upscale *= 10

    kwargs["upscale"] = upscale + 4 if not upscale is None else None

    llh2012 = ic86_2012(**kwargs)
    llh.add_sample("IC86-II", llh2012)

    print("{0:>80s}".format("[done]"))
    return llh


def load_6yr_mese(*args, **kwargs):
    print("Loading 6year and MESE")

    # load all years from IC40 to IC86.2011
    llh = load_4yr_ic79b_mese(*args, **kwargs)

    upscale = kwargs.pop("upscale", None)
    if upscale is not None:
        upscale *= 10

    kwargs["upscale"] = upscale + 5 if not upscale is None else None
    kwargs["livetime"] = 691.

    llh2012 = ic86_2012(no_mese=True, **kwargs)
    llh.add_sample("IC86-II", llh2012)

    print("{0:>80s}".format("[done]"))
    return llh

def load_7yr_nospline_79(*args, **kwargs):  
    print("Loading 7year with unsplined v. of IC79")
    # load all years from IC40 to IC86.2011
    
    llh = load_4yr(years="IC40+IC59+IC79+IC86",*args, **kwargs)
    upscale = kwargs.pop("upscale", None)
    if upscale is not None:
        upscale *= 10
     
    kwargs["upscale"] = upscale + 4 if not upscale is None else None
     
    llh2012 = ic86_2012(dataset=[2012, 2013, 2014], **kwargs)

    llh.add_sample("IC86-121314", llh2012)

    print("{0:>80s}".format("[done]"))
    return llh

def load_7yr_mese_nospline_79(*args, **kwargs):
        print("Loading 7year w/ unsplined 79 and MESE")

        # load through-going events
        llh = load_7yr_nospline_79(*args, no_mese=True, **kwargs)

        upscale = kwargs.pop("upscale", None)
        if upscale is not None:
            upscale *= 10
        
        kwargs["upscale"] = upscale + 5 if not upscale is None else None
        llhMESE = mese_followup(*args, **kwargs)

        llh.add_sample("MESE", llhMESE)
     
        print("{0:>80s}".format("[done]"))
        return llh

def load_7yr(*args, **kwargs):  
    print("Loading 7year")
    # load all years from IC40 to IC86.2011
    
    llh = load_4yr_ic79b(*args, **kwargs)
    
    upscale = kwargs.pop("upscale", None)
    if upscale is not None:
        upscale *= 10
     
    kwargs["upscale"] = upscale + 4 if not upscale is None else None
     
    llh2012 = ic86_2012(dataset=[2012, 2013, 2014], **kwargs)

    llh.add_sample("IC86-121314", llh2012)

    print("{0:>80s}".format("[done]"))
    return llh
    
   
def load_7yr_mese(*args, **kwargs):
        print("Loading 7year and MESE")

        # load through-going events
        llh = load_7yr(*args, no_mese=True, **kwargs)

        upscale = kwargs.pop("upscale", None)
        if upscale is not None:
            upscale *= 10
        
        kwargs["upscale"] = upscale + 5 if not upscale is None else None
        llhMESE = mese_followup(*args, **kwargs)

        llh.add_sample("MESE", llhMESE)
     
        print("{0:>80s}".format("[done]"))
        return llh
        



def monte_carlo(llh, suffix=""):

    if isinstance(llh, StackingMultiPointSourceLLH):
        mc = dict()
        for enum, name in llh._enum.iteritems():
            if re.findall(r"IC86-\d+", name):
                # for 2012ff years only 2012 MC
                name = "IC86-2012"

            if "MESE" in llh._enum.values() and re.findall(r"IC79|IC86", name):
                name += "_noMESE"
            
            fpath = os.path.join(path,"{0:s}{1:s}_corrected_MC.npy".format(
                                      name, suffix))

            if name == "MESE":
                fpath = fpath.replace("_corrected", "")
          
            print("For injection, load {0:s} (sample {1:s})".format(fpath, llh._enum[enum]))

            arr = np.load(fpath)
            if name == "MESE" and not "dist" in llh._sams[enum].exp.dtype.names:
                arr = drop_fields(arr, "dist", usemask=False)

            mc[enum] = drop_fields(arr, "BDT", usemask=False)

    elif isinstance(llh, str):
        fpath = os.path.join(path,"{0:s}{1:s}_corrected_MC.npy".format(
                                             llh, suffix))
        print("For injection, load {0:s}".format(fpath))

        mc = np.load(fpath)
        mc = drop_fields(np.load(fpath), "BDT", usemask=False)

    else:
        raise ValueError("Need StackingMultiPointSourceLLH or string")

    return mc
