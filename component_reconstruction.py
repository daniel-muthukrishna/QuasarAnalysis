import numpy as np
import matplotlib.pyplot as plt
import pickle
from get_spectra import get_boss_dr12_spec


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

    return waves, reconstructedSpectra, names, balFlags


def get_sdss_spectra(names):
    waves, fluxes = [], []
    for name in names[0:10]:
        wave, dw, flux, err, mask = get_boss_dr12_spec(name)
        waves.append(wave)
        fluxes.append(flux)

    return waves, fluxes


def spectra_dict(componentsFile, waveFile, weightsFile):
    reconWave, reconFluxes, names, balFlags = reconstruct_spectra(componentsFile, waveFile, weightsFile)
    sdssWaves, sdssFluxes = get_sdss_spectra(names)

    spectraDict = {}
    numSpectra = len(names)

    for i in range(4):
        spectraDict[names[i]] = {'reconWave': reconWave, 'reconFlux': reconFluxes[i], 'balFlag': balFlags[i],
                                 'sdssWave': sdssWaves[i], 'sdssFlux': sdssFluxes[i]}

    return spectraDict


def save_spectra(componentsFile, waveFile, weightsFile):
    spectra = spectra_dict(componentsFile, waveFile, weightsFile)
    with open('spectra.pickle', 'wb') as f:
        pickle.dump(spectra, f, pickle.HIGHEST_PROTOCOL)


def plot_spectrum(name, spectra):
    plt.figure(name)
    plt.plot(spectra[name]['sdssWave'], spectra[name]['sdssFlux'], label='DR12_spectrum')
    plt.plot(spectra[name]['reconWave'], spectra[name]['reconFlux'], label='reconstruction')
    plt.legend()
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Flux')
    bal = 'BAL' if spectra[name]['balFlag'] else 'non-BAL'
    plt.title('{0}_{1}'.format(name, bal))


if __name__ == '__main__':
    spectra1 = spectra_dict(componentsFile='dm_6c_16003000_171024.comp', waveFile='wav_16003000.dat', weightsFile='dm_hbal_weights.dat')
    names1 = list(spectra1.keys())
    for idx in range(3):
        plot_spectrum(names1[idx], spectra1)
    plt.show()
