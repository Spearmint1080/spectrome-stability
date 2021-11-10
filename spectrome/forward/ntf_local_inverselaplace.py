import numpy as np
from mpmath import *
from sympy import *
from scipy import signal

def ntf_local_inverselaplace(parameters,tt,Q,noisein):
    tau_e = parameters["tau_e"]
    tau_i = parameters["tau_i"]
    gei = parameters[
        "gei"
    ]  # excitatory-inhibitory synaptic conductance as ratio of E-E syn
    gii = parameters[
        "gii"
    ]  # inhibitory-inhibitory synaptic conductance as ratio of E-E syn

    gee = 1
    

    dt = 0.01
    trange = np.arange(0,10,dt)
    scaleterm = 10 #Scale for noise is 1/sqrt(dt)
    nr = np.random.normal(scale = scaleterm, size = len(trange))
    windowed = nr

    def fee(s):
        s1 = float(mp.re(s)) + 1j*float(mp.im(s))
        Fe = np.divide(1 / tau_e ** 2, (s1 + 1 / tau_e) ** 2)
        Fi = np.divide(1 / tau_i ** 2, (s1 + 1 / tau_i) ** 2)
        Hed = (1 + (Fe * Fi * gei)/(tau_e * (s1 + Fi*gii/tau_i)))/(s1 + Fe*gee/tau_e + (Fe * Fi * gei)**2/(tau_e * tau_i * (s1 + Fi * gii / tau_i)))

        nlp = 0
        if noisein == "noise":
            for i in range(len(trange)):
                nlp += np.exp(-s1*trange[i])*Q*windowed[i]*dt
        if noisein == "no noise":
            nlp = Q
        return Hed*nlp

    def fii(s):
        s1 = float(mp.re(s)) + 1j*float(mp.im(s))
        Fe = np.divide(1 / tau_e ** 2, (s1 + 1 / tau_e) ** 2)
        Fi = np.divide(1 / tau_i ** 2, (s1 + 1 / tau_i) ** 2)
        Hid = (1 - (Fe * Fi * gei)/(tau_i * (s1 + Fe*gee/tau_e)))/(s1 + Fi * gii/tau_i + (Fe * Fi * gei)**2/(tau_e * tau_i * (s1 + Fe*gee / tau_e)))

        nlp = 0
        if noisein == "noise":
            for i in range(len(trange)):
                nlp += np.exp(-s1*trange[i])*Q*windowed[i]*dt
        if noisein == "no noise":
            nlp = Q
        return Hid*nlp

    et = [float(invertlaplace(fee,t)) for t in tt]
    it = [float(invertlaplace(fii,t)) for t in tt]

    return et,it