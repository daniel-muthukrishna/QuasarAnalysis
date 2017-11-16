import numpy as np
from sklearn.decomposition import PCA


def get_data_spectra(spectra):
    sdssFluxes = []
    names = list(spectra.keys())
    for name in names:
        sdssFlux = spectra[name]['sdssFlux']
        sdssFluxes.append(sdssFlux)

    return np.array(sdssFluxes)


def pca_analysis(spectra):
    sdssFluxes = get_data_spectra(spectra)
    pca = PCA(n_components=5)
    pca.fit(sdssFluxes)




