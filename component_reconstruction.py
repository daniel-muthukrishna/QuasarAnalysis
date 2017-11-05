import numpy as np
import matplotlib.pyplot as plt
import pickle
from get_spectra import get_sdss_dr12_spectrum, load_spectra_filenames
from scipy.stats import pearsonr


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
    redshifts, fluxes = [], []
    filenamesDict = load_spectra_filenames()
    for name in names:
        flux, z = get_sdss_dr12_spectrum(name, filenamesDict)
        redshifts.append(z)
        fluxes.append(flux)

    return fluxes, redshifts


def spectra_dict(componentsFile, waveFile, weightsFile):
    reconWave, reconFluxes, names, balFlags = reconstruct_spectra(componentsFile, waveFile, weightsFile)
    sdssFluxes, sdssRedshifts = get_sdss_spectra(names)

    spectraDict = {}
    numSpectra = len(names)

    for i in range(numSpectra):
        spectraDict[names[i]] = {'reconWave': reconWave, 'reconFlux': reconFluxes[i], 'balFlag': balFlags[i],
                                 'sdssWave': reconWave, 'sdssFlux': sdssFluxes[i], 'sdssRedshifts': sdssRedshifts[i]}

    return spectraDict


def deredshift_spectrum(wave, z):
    return wave/(z + 1.)


def load_spectra(savedSpectra='spectra.pickle'):
    with open(savedSpectra, 'rb') as f:
        spectra = pickle.load(f)

    return spectra

def plot_spectrum(name, spectra):
    plt.figure(name)
    print(len(spectra[name]['sdssWave']), len(spectra[name]['sdssFlux']), len(spectra[name]['reconFlux'])) 
    plt.plot(spectra[name]['sdssWave'], spectra[name]['sdssFlux'], label='DR12_spectrum')
    plt.plot(spectra[name]['reconWave'], spectra[name]['reconFlux'], label='reconstruction')
    plt.legend()
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Flux')
    bal = 'BAL' if spectra[name]['balFlag'] else 'non-BAL'
    plt.title('{0}_{1}'.format(name, bal))


def compare_reconstruction_and_data(name, spectra):
    sdssFlux = spectra[name]['sdssFlux']
    reconFlux = spectra[name]['reconFlux']
    pearsonCorr = pearsonr(sdssFlux, reconFlux)[0]

    return pearsonCorr


def compare_all(spectra):
    pearsonVals = []
    names = list(spectra.keys())
    for i in range(len(names)):
        pearson = compare_reconstruction_and_data(names[i], spectra)
        pearsonVals.append(pearson)
        bal = 'BAL' if spectra[names[i]]['balFlag'] else 'non-BAL'
        print(pearson, bal)
        if pearson < 0.5:
            print("################## {0} ".format(names[i]))
            #plot_spectrum(names[i], spectra)

    meanPearson = np.mean(pearsonVals)
    stdPearson = np.std(pearsonVals)
    print("Average Pearson: {0}, {1}".format(meanPearson, stdPearson))
    
    #for name, pearson in zip(names, pearsonVals):
        

if __name__ == '__main__':
    spectra1 = load_spectra()
    compare_all(spectra1)
    #names1 = list(spectra1.keys())
    #for idx in range(3):
     #   plot_spectrum(names1[idx], spectra1)
    plt.show()
