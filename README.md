# spectrome-stability

`Spectrome` is a combination of the words "spectrum" and "connectome". This package is the collection of codes that constructed the analysis for the preprint ["Spectral graph theory of brain oscillations - revisited and improved"](https://www.biorxiv.org/content/10.1101/2021.09.28.462078v1). This repository is developed based on the original model's [repository](https://github.com/Raj-Lab-UCSF/spectrome).

The spectral graph model (SGM) is a brain structure-function model that simulates brain activity power spectrum given a structural connectome. The model is linear, low-dimensional, and provides an analytical relationship between the brain's structural and functional patterns.

## Citation:
<!-- The code in this repository is used for the analysis as shown in: Parul Verma, Srikantan Nagarajan, and Ashish Raj. “Spectral Graph Theory of Brain Oscillations - revisited and improved” (https://www.biorxiv.org/content/10.1101/2021.09.28.462078v1). If you found this useful, please cite the following: -->

<!-- ```
@article {verma2021spectral,
	author = {Verma, Parul and Nagarajan, Srikantan and Raj, Ashish},
	title = {Spectral graph theory of brain oscillations -- revisited and improved},
	elocation-id = {2021.09.28.462078},
	year = {2021},
	doi = {10.1101/2021.09.28.462078},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2021/09/30/2021.09.28.462078},
	eprint = {https://www.biorxiv.org/content/early/2021/09/30/2021.09.28.462078.full.pdf},
	journal = {bioRxiv}
}
``` -->

## Abstract:
We explore the stability and dynamic properties of a hierarchical, linearized, and analytic spectral graph model for neural oscillations that integrates the structuring wiring of the brain. Previously we have shown that this model can accurately capture the frequency spectra and the spatial patterns of the alpha and beta frequency bands obtained from magnetoencephalography recordings without regionally varying parameters. Here, we show that this macroscopic model based on long-range excitatory connections exhibits dynamic oscillations with a frequency in the alpha band even without any oscillations implemented at the mesoscopic level. We show that depending on the parameters, the model can exhibit combinations of damped oscillations, limit cycles, or unstable oscillations. We determined bounds on model parameters that ensure stability of the oscillations simulated by the model. Finally, we estimated time-varying model parameters to capture the temporal fluctuations in magnetoencephalography activity. We show that a dynamic spectral graph modeling framework with a parsimonious set of biophysically interpretable model parameters can thereby be employed to capture oscillatory fluctuations observed in electrophysiological data in various brain states and diseases.

## Set-up:

First clone the environment to your computer, either download this repo as a `.zip` file or `git clone https://github.com/Raj-Lab-UCSF/spectrome.git`.

Set up a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) if you do not have all the packages/compatible versions. The list of dependencies is listed in `environment.yml`.

Set-up environment using conda, detailed instructions can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html). Or after cloning this repository, go to the repo by typing `cd spectrome` and then typing:
`conda env create -f environment.yml`

If conda complains about not finding packages/libraries, make sure `conda-forge` is in the list of channels being searched by `conda`.
You may add `conda-forge` to your list of channels with the command: `conda config --add channels conda-forge`.

The default name of the environment is `spectrome`, activate the environment with `source activate spectrome`, and deactivate with `source deactivate` or `conda deactivate`.

If you want to be able to run `spectrome` from anywhere, just add it's path to your PYTHONPATH. For instance, if you downloaded `spectrome` to `~/Documents/spectrome` do `export PYTHONPATH=$PYTHONPATH:~/Documents/spectrome`. You may have to restart your terminal to make sure this change takes effect.

After completing the set-up for conda environment and `spectrome` path, you may go to the `spectrome` folder and type `jupyter notebook` or `jupyter lab` in your terminal to run the Jupyter notebooks.

## Files:
 - `../spectrome/notebooks`: contains three jupyter notebooks, `run_model_example.ipynb` is the basic simulation of frequency spectrums with default parameters for the HCP template connectome. `spatialcorrelation.ipynb` looks at the spatial correlations between the eigenmodes and the empirical spectra, and `reproduce_MEG_modeled_spectra.ipynb` compares optimized modeled spectra with the MEG spectra.

 - `../spectrome/data`: contains intermediate data.
    - `mean80_fibercount/length.csv`: HCP template connectome and distance matrix.
    - `individual_connectomes_reordered.nc`: individual subject's connectivity matrices (N = 36).
    - `individual_psd_reordered_matlab.nc`: individual subject's MEG spectra (N = 36).
    - `individual_wavelet_reordered_smooth.csv`: individual subject's MEG Morlet wavelet decomposition (N = 36)
