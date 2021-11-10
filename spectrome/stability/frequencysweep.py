""" Estimating stability of macro model by root finding.
    Here, we want to find a root of jw + Fe(w)*(I-alpha*C)/tauC
    This becomes a multiplication of jw + Fe(w)*lambda_i(w)/tauC, multiplication over i where lambda_i are eigenvalues
    Thus, we need to find roots of jw + lambda_i(w)/(tauC*(tau_e*jw+1)**2)
    This gives tauC*jw*(-w**2*tau_e**2 + 2*tau_e*jw + 1). We need this to be zero.
    Separating out real and imaginary parts of this expression gives:
    Re: -2*tau_e*tauC*w**2 + re(lambda_i(w)) = 0
    Im: -tau_e**2*tauC*w**3 + tauC*w + im(lambda_i(w)) = 0
    Since w should be real, this also implies that we will only consider lambda_i(w) 
    with positive real parts (based on the equation for real part of the expression)
"""

import numpy as np

def network_transfer(parameters, brain, orgparameters):
    """Estimating stability of macro model when alpha > 1 using root finding

    Args:
        parameters (dict): w, tauC
        brain (Brain): specific brain to calculate NTF
        orgparameters (dict): original parameters of NTF that won't be changing

    Returns:
        Condensed version of characteristic equations for root finding

    """
    C = brain.reducedConnectome
    D = brain.distance_matrix

    w = parameters[0]
    tauC = parameters[1]

    tau_e = orgparameters["tau_e"]
    speed = orgparameters["speed"]
    alpha = orgparameters["alpha"]

    # define sum of degrees for rows and columns for laplacian normalization
    rowdegree = np.transpose(np.sum(C, axis=1))
    coldegree = np.sum(C, axis=0)
    qind = rowdegree + coldegree < 0.2 * np.mean(rowdegree + coldegree)
    rowdegree[qind] = np.inf
    coldegree[qind] = np.inf

    nroi = C.shape[0]

    Tau = 0.001 * D / speed
    Cc = C * np.exp(-1j * Tau * w)

    # Eigen Decomposition of Complex Laplacian Here
    L1 = np.identity(nroi)
    L2 = np.divide(1, np.sqrt(np.multiply(rowdegree, coldegree)) + np.spacing(1))
    L = L1 - alpha * np.matmul(np.diag(L2), Cc)

    d, v = np.linalg.eig(L)  # decomposition with scipy.linalg.eig
    eig_ind = np.argsort(np.abs(d))  # sorting in ascending order and absolute value
    eig_vec = v[:, eig_ind]  # re-indexing eigen vectors according to sorted index
    eig_val = d[eig_ind]  # re-indexing eigen values with same sorted index

    eigsub = eig_val[eig_val.real>=0]

    ind = np.argmin(np.abs(2*tau_e*tauC*w**2 - eigsub.real))

    return [(2*tau_e*tauC*w**2 - eigsub[ind].real), 
            (-tau_e**2*tauC*w**3 + tauC*w + eigsub[ind].imag)]