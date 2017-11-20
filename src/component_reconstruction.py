import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import medfilt
from src.get_spectra import get_sdss_dr12_spectrum, load_spectra_filenames


def get_components(componentsFile):
    components = np.loadtxt(componentsFile)

    return components


def get_wavelengths(waveFile):
    wavelengths = np.loadtxt(waveFile)

    return wavelengths


def get_weights(weightsFile):
    weightInfo = np.genfromtxt(weightsFile, dtype=str)
    names = weightInfo[:, 0]
    weights = weightInfo[:, 1:-1].astype(float)
    balFlags = weightInfo[:, -1].astype(int)

    return names, weights, balFlags


def reconstruct_spectra(componentsFile, waveFile, weightsFile):
    comps = get_components(componentsFile)
    waves = get_wavelengths(waveFile)
    names, weights, balFlags = get_weights(weightsFile)

    # components * weights
    reconstructedSpectra = np.transpose(comps.dot(np.transpose(weights)))

    return waves, reconstructedSpectra, names, balFlags, weights, comps


def get_sdss_spectra(names):
    redshifts, fluxes = [], []
    filenamesDict = load_spectra_filenames()
    for name in names[0:]:
        flux, z = get_sdss_dr12_spectrum(name, filenamesDict)
        redshifts.append(z)
        fluxes.append(flux)

    return fluxes, redshifts


def spectra_dict(componentsFile, waveFile, weightsFile):
    reconWave, reconFluxes, names, balFlags, weights, comps = reconstruct_spectra(componentsFile, waveFile, weightsFile)
    sdssFluxes, sdssRedshifts = get_sdss_spectra(names)

    spectraDict = {}
    numSpectra = len(names)

    for i in range(numSpectra):
        spectraDict[names[i]] = {'reconWave': reconWave, 'reconFlux': reconFluxes[i], 'balFlag': balFlags[i],
                                 'sdssWave': reconWave, 'sdssFlux': sdssFluxes[i], 'sdssRedshifts': sdssRedshifts[i],
                                 'weights': weights[i]}

    return spectraDict


def deredshift_spectrum(wave, z):
    return wave/(z + 1.)


def smooth(flux):
    return medfilt(flux, kernel_size=7)


def plot_spectrum(name, spectra, addToTitle=''):
    plt.figure()
    plt.plot(spectra[name]['sdssWave'], smooth(spectra[name]['sdssFlux']), label='DR12_spectrum')
    plt.plot(spectra[name]['reconWave'], spectra[name]['reconFlux'], label='reconstruction')
    plt.legend()
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Flux')
    bal = 'BAL' if spectra[name]['balFlag'] else 'non-BAL'
    plt.title('{0}_{1}_{2}'.format(name, bal, addToTitle))


def reconstruct_each_component(spectra, comps):
    flux = {}

    names = list(spectra.keys())
    for name in names[0:10]:
        flux[name] = {}
        weights = spectra[name]['weights']
        for compNum in range(len(weights)):  # loop though num of comps
            flux[name][compNum] = comps[:, compNum] * weights[compNum]

    return flux


def plot_each_component(spectra, comps, mean=0):
    names = list(spectra.keys())
    numComps = len(spectra[names[0]]['weights'])
    wave = spectra[names[0]]['reconWave']
    flux = reconstruct_each_component(spectra, comps)

    for name in names[0:10]:
        plt.figure()
        plt.plot(wave, smooth(spectra[name]['sdssFlux']), label='sdss')
        plt.plot(wave, spectra[name]['reconFlux'], label='recon')
        bal = 'BAL' if spectra[name]['balFlag'] else 'non-BAL'
        for compNum in range(numComps):
            plt.plot(wave, flux[name][compNum], label=compNum+1)
        plt.plot(wave, np.sum([flux[name][compNum] for compNum in range(numComps)], axis=0) + mean, label='$1\_to\_{0}$'.format(numComps))
        plt.xlabel('Wavelength ($\AA$)')
        plt.ylabel("Median Fractional Difference")
        plt.title("{0}: {1}".format(bal, name))
        plt.legend()

