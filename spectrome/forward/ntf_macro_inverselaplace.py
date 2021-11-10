import numpy as np
import mpmath as mp
from ..forward import network_transfer_laplace

def ntf_macro_inverselaplace(brain,parameters,tt,noisein):
    C = brain.reducedConnectome
    D = brain.distance_matrix

    tau_e = parameters["tau_e"]
    tau_i = parameters["tau_i"]
    speed = parameters["speed"]
    gei = parameters[
        "gei"
    ]  # excitatory-inhibitory synaptic conductance as ratio of E-E syn
    gii = parameters[
        "gii"
    ]  # inhibitory-inhibitory synaptic conductance as ratio of E-E syn
    tauC = parameters["tauC"]  #  tauC = 0.5*tau_e
    alpha = parameters["alpha"]
    gee = 1

    # Defining some other parameters used:
    zero_thr = 0.05
    # a = 0.5  # fraction of signal at a node that is recurrent excitatory

    # define sum of degrees for rows and columns for laplacian normalization
    rowdegree = np.transpose(np.sum(C, axis=1))
    coldegree = np.sum(C, axis=0)
    qind = rowdegree + coldegree < 0.2 * np.mean(rowdegree + coldegree)
    rowdegree[qind] = np.inf
    coldegree[qind] = np.inf

    nroi = C.shape[0]

    L1 = np.identity(nroi)
    L2 = np.divide(1, np.sqrt(np.multiply(rowdegree, coldegree)) + np.spacing(1))

    Tau = 0.001 * D / speed
    # Tau = 0

    lpt = np.zeros((nroi,len(tt)))

    # def f(s):
    #     s1 = float(mp.re(s)) + 1j*float(mp.im(s))
    #     print(s1)
    #     model_out = network_transfer_laplace.ntf_laplace(C,L1,L2,Tau,nroi,alpha,tau_e,tau_i,tauC,gei,gee,gii,zero_thr,s1,noisein)
    #     return model_out[0][0]
    # return [float(mp.invertlaplace(f,t)) for t in tt]

    dt = 0.008
    trange = np.arange(0,100,dt)
    # Q = 1 #Scale the noise term appropriately
    scaleterm = 1/np.sqrt(dt) #Scale for noise is 1/sqrt(dt)
    nr = np.random.normal(scale = scaleterm, size = len(trange))

    for i in range(nroi):
        def f(s):
            s1 = float(mp.re(s)) + 1j*float(mp.im(s))
            model_out = network_transfer_laplace.ntf_laplace(C,L1,L2,Tau,nroi,alpha,tau_e,tau_i,tauC,gei,gee,gii,s1,noisein,trange,nr,dt)
            return model_out[i][0]
        lpt[i,:] = [float(mp.invertlaplace(f,t)) for t in tt]

    return lpt