import numpy as np
# import mpmath as mp

def ntf_laplace(C,L1,L2,Tau,nroi,alpha,tau_e,tau_i,tauC,gei,gee,gii,s,noisein,trange,nr,dt):
  
    nroi = C.shape[0]
    K = nroi

#   Cortical model
    Fe = np.divide(1 / tau_e ** 2, (s + 1 / tau_e) ** 2)
    Fi = np.divide(1 / tau_i ** 2, (s + 1 / tau_i) ** 2)

    Hed = (1 + (Fe * Fi * gei)/(tau_e * (s + Fi * gii/tau_i)))/(s + Fe * gee/tau_e + (Fe * Fi * gei)**2/(tau_e * tau_i * (s + Fi * gii / tau_i)))
    
    Hid = (1 - (Fe * Fi * gei)/(tau_i * (s + Fe * gee/tau_e)))/(s + Fi * gii/tau_i + (Fe * Fi * gei)**2/(tau_e * tau_i * (s + Fe * gee / tau_e)))

    Q = 1
    windowed = nr
    nlp = 0
    if noisein == "noise":
        for i in range(len(trange)):
            nlp += np.exp(-s*trange[i])*Q*windowed[i]*dt
        Htotal = (Hed + Hid)*nlp # Complete system with local input and noise
    if noisein == "no noise local":
        Htotal = (Hed + Hid)*Q
    if noisein == "no noise":
        nlp = Q
        Htotal = nlp # Only frequency response of macro
    
    Cc = C * np.exp(-Tau * s)

    L = (s + Fe / tauC) * L1 - (alpha * Fe * np.matmul(np.diag(L2), Cc)) / tauC

    Xlocal = Htotal * np.ones((nroi,1))

    model_out2 = np.matmul(np.linalg.inv(L),Xlocal)

    return model_out2

