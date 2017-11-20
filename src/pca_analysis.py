import numpy as np
from sklearn.decomposition import PCA, FastICA


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
    weights = pca.fit_transform(sdssFluxes)
    comps = pca.components_
    reconFluxes = pca.inverse_transform(weights)  # weights.dot(comps) + pca.mean_
    loss = ((sdssFluxes - reconFluxes) ** 2).mean()

    spectraPCA = pca_spectra_dict(spectra, weights, reconFluxes)

    return spectraPCA, comps


def pca_spectra_dict(spectra, weights, reconFluxes):
    spectraPCA = spectra
    names = list(spectraPCA.keys())
    numSpectra = len(names)

    for i in range(numSpectra):
        spectraPCA[names[i]]['weights'] = weights[i]
        spectraPCA[names[i]]['reconFlux'] = reconFluxes[i]

    return spectraPCA


def ica_analysis(spectra):
    sdssFluxes = get_data_spectra(spectra)
    ica = FastICA(n_components=5)
    ica.fit(sdssFluxes)
    weights = ica.fit_transform(sdssFluxes)
    comps = ica.components_
    reconFluxes = ica.inverse_transform(weights)  # weights.dot(comps) + ica.mean_
    loss = ((sdssFluxes - reconFluxes) ** 2).mean()

    spectraICA = ica_spectra_dict(spectra, weights, reconFluxes)

    return spectraICA, comps


def ica_spectra_dict(spectra, weights, reconFluxes):
    spectraICA = spectra
    names = list(spectraICA.keys())
    numSpectra = len(names)

    for i in range(numSpectra):
        spectraICA[names[i]]['weights'] = weights[i]
        spectraICA[names[i]]['reconFlux'] = reconFluxes[i]

    return spectraICA


