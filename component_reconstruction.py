import numpy as np
import matplotlib.pyplot as plt
import pickle
from get_spectra import get_sdss_dr12_spectrum, load_spectra_filenames
from scipy.stats import pearsonr, chisquare
from scipy.signal import medfilt
from sklearn.cluster import KMeans, k_means


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

    return waves, reconstructedSpectra, names, balFlags, weights


def get_sdss_spectra(names):
    redshifts, fluxes = [], []
    filenamesDict = load_spectra_filenames()
    for name in names[0:]:
        flux, z = get_sdss_dr12_spectrum(name, filenamesDict)
        redshifts.append(z)
        fluxes.append(flux)


    return fluxes, redshifts


def spectra_dict(componentsFile, waveFile, weightsFile):
    reconWave, reconFluxes, names, balFlags, weights = reconstruct_spectra(componentsFile, waveFile, weightsFile)
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


def load_spectra(savedSpectra='spectra.pickle'):
    with open(savedSpectra, 'rb') as f:
        spectra = pickle.load(f)

    return spectra


def plot_spectrum(name, spectra, addToTitle=''):
    plt.figure(name)
    plt.plot(spectra[name]['sdssWave'], medfilt(spectra[name]['sdssFlux'], kernel_size=7), label='DR12_spectrum')
    plt.plot(spectra[name]['reconWave'], spectra[name]['reconFlux'], label='reconstruction')
    plt.legend()
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Flux')
    bal = 'BAL' if spectra[name]['balFlag'] else 'non-BAL'
    plt.title('{0}_{1}_{2}'.format(name, bal, addToTitle))


def compare_reconstruction_and_data(name, spectra):
    sdssFlux = medfilt(spectra[name]['sdssFlux'], kernel_size=7)
    reconFlux = spectra[name]['reconFlux']
    pearsonCorr = round(pearsonr(sdssFlux, reconFlux)[0], 3)
    chi2 = round(chisquare(sdssFlux, reconFlux)[0], 0)

    return pearsonCorr, chi2


def compare_all(spectra):
    pearsonVals = []
    removeNames = []
    names = list(spectra.keys())
    for name in names:
        pearson, chi2 = compare_reconstruction_and_data(name, spectra)
        pearsonVals.append(pearson)
        bal = 'BAL' if spectra[name]['balFlag'] else 'non-BAL'
        if pearson < 0.66 or chi2 > 300:
            removeNames.append(name)
            print(pearson, bal, chi2)
            print("################## {0} ".format(name))
            # plot_spectrum(name, spectra, addToTitle='PC={0}_{1}'.format(pearson, chi2))

    meanPearson = np.mean(pearsonVals)
    stdPearson = np.std(pearsonVals)
    print("Average Pearson: {0}, {1}".format(meanPearson, stdPearson))
    
    return removeNames


def clustering(spectra, removeNames):
    kmeans = KMeans(n_clusters=2)
    names = list(spectra.keys())
    weightsArray = []
    for name in names:
        if name not in removeNames:
            weightsArray.append(spectra[name]['weights'])
    weightsArray = np.array(weightsArray)
    kmeans.fit(weightsArray)


if __name__ == '__main__':
    spectra1 = load_spectra()
    removeNames1 = compare_all(spectra1)
    clustering(spectra1, removeNames1)
    # names1 = list(spectra1.keys())
    # for idx in range(10):
    #     plot_spectrum(names1[idx], spectra1)
    plt.show()
