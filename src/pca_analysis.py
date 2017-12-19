import numpy as np
from sklearn.decomposition import PCA, FastICA
import copy


def get_data_spectra(spectra):
    sdssFluxes = []
    names = list(spectra.keys())
    for name in names:
        sdssFlux = spectra[name]['sdssFlux']
        sdssFluxes.append(sdssFlux)

    return np.array(sdssFluxes)


def pca_analysis(spectra, nComps=5):
    sdssFluxes = get_data_spectra(spectra)
    pca = PCA(n_components=nComps)
    pca.fit(sdssFluxes)
    weights = pca.fit_transform(sdssFluxes)
    comps = pca.components_.transpose()
    reconFluxes = pca.inverse_transform(weights)  # weights.dot(pca.components_) + pca.mean_
    loss = ((sdssFluxes - reconFluxes) ** 2).mean()

    spectraPCA = pca_spectra_dict(spectra, weights, reconFluxes)

    return spectraPCA, comps, pca.mean_


def pca_spectra_dict(spectra, weights, reconFluxes):
    spectraPCA = copy.deepcopy(spectra)
    names = list(spectraPCA.keys())
    numSpectra = len(names)

    for i in range(numSpectra):
        spectraPCA[names[i]]['weights'] = weights[i]
        spectraPCA[names[i]]['reconFlux'] = reconFluxes[i]

    return spectraPCA


def ica_analysis(spectra, nComps=5):
    sdssFluxes = get_data_spectra(spectra)
    ica = FastICA(n_components=nComps)
    ica.fit(sdssFluxes)
    weights = ica.fit_transform(sdssFluxes)
    comps = ica.components_.transpose()
    reconFluxes = ica.inverse_transform(weights)  # weights.dot(ica.mixing_) + ica.mean_
    loss = ((sdssFluxes - reconFluxes) ** 2).mean()

    spectraICA = ica_spectra_dict(spectra, weights, reconFluxes)

    return spectraICA, ica.mixing_, ica.mean_


def ica_spectra_dict(spectra, weights, reconFluxes):
    spectraICA = copy.deepcopy(spectra)
    names = list(spectraICA.keys())
    numSpectra = len(names)

    for i in range(numSpectra):
        spectraICA[names[i]]['weights'] = weights[i]
        spectraICA[names[i]]['reconFlux'] = reconFluxes[i]

    return spectraICA


