""" Computing and sorting eigenmodes for alpha and beta band spatial correlations"""
from ..forward import network_transfer_spatialcorrelation_notsorted as nt_ns
import numpy as np
from scipy.stats import pearsonr

def run_local_coupling_forward_Xk(brain, params, freqs, PSD, SC, rois_with_MEG, band):

    """Network Transfer Function for spectral graph model.

    Args:
        brain (Brain): specific brain to calculate NTF
        parameters (dict): parameters for ntf. We shall keep this separate from Brain
        for now, as we want to change and update according to fitting.
        frequency (float): frequency at which to calculate NTF
        PSD: PSD of a subject to compute spatial correlation 
        SC: Number of eigenvectors to include
        rois_with_MEG: rois for which MEG spectra is available
        band: "alpha" or "beta"

    Returns:
        spcorr2 (numpy asarray): spatial correlation of summed eigenmodes
        eigvec_sorted (numpy asarray): sorted eigenmodes
        summed_PSD (numpy asarray): PSD summed for the frequency band of interest
        eig_ind (numpy asarray):  index of sorted eigenmodes
    """

    if band == "alpha":
        freqband = np.where((freqs>=8) & (freqs<=12))[0]
    if band == "beta":
        freqband = np.where((freqs>=13) & (freqs<=25))[0]

    eigvec_ns = np.zeros((len(rois_with_MEG),SC,len(freqband)))

    for i in range(len(freqband)):
        w = 2 * np.pi * freqs[freqband[i]]
        eigenvectors_ns = nt_ns.network_transfer_Xk(
            brain, params, w, rois_with_MEG, SC
        )
        eigvec_ns[:,:,i] = eigenvectors_ns

    spcorr = np.zeros(SC)

    eigvec_ns_summed = np.sum(eigvec_ns[:,:,:],axis = 2)


    summed_PSD = np.sum(PSD[:,freqband], axis = 1)

    for i in range(SC):
        spcorr[i] = pearsonr(summed_PSD, eigvec_ns_summed[:,i])[0]
    
    eig_ind = np.argsort(-spcorr)

    eigvec_sorted = eigvec_ns_summed[:,eig_ind]

    
    spcorr2 = np.zeros(SC)
    for i in range(SC):
        spcorr2[i] = pearsonr(summed_PSD, np.sum(eigvec_sorted[:,0:(i+1)], axis = 1))[0]

    return spcorr2, eigvec_sorted, summed_PSD, eig_ind