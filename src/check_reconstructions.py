import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, chisquare
from src.component_reconstruction import smooth


def compare_reconstruction_and_data(name, spectra):
    """ Calculate pearson correlation between data and reconstruction """
    sdssFlux = smooth(spectra[name]['sdssFlux'])
    reconFlux = spectra[name]['reconFlux']
    pearsonCorr = round(pearsonr(sdssFlux, reconFlux)[0], 3)
    chi2 = 0  # round(chisquare(sdssFlux, reconFlux)[0], 0)

    return pearsonCorr, chi2


def compare_all(spectra):
    """ Remove spectra with a low pearson correlation """
    pearsonVals = []
    removeNames = []
    names = list(spectra.keys())
    for name in names:
        pearson, chi2 = compare_reconstruction_and_data(name, spectra)
        pearsonVals.append(pearson)
        bal = 'BAL' if spectra[name]['balFlag'] else 'non-BAL'
        if pearson < 0.66 or chi2 > 300:
            removeNames.append(name)
            # print(name, pearson, bal, chi2)
            # plot_spectrum(name, spectra, addToTitle='PC={0}_{1}'.format(pearson, chi2))

    meanPearson = np.mean(pearsonVals)
    stdPearson = np.std(pearsonVals)
    print("Average Pearson: {0}, {1}".format(meanPearson, stdPearson))

    return removeNames


def frac_diff(data, recon):
    fracDiff = (data - recon)/recon
    return fracDiff


def frac_diff_all(wave, spectra, removeNames, saveDir):
    fracDiffs = {'All': [], 'BAL': [], 'nonBAL': []}
    lossSquared = []
    names = list(spectra.keys())
    for name in names:
        if name not in removeNames:
            sdssFlux = spectra[name]['sdssFlux']
            reconFlux = spectra[name]['reconFlux']
            fracDiff = frac_diff(sdssFlux, reconFlux)
            lossSquared.append((sdssFlux-reconFlux)**2)
            fracDiffs['All'].append(fracDiff)
            if spectra[name]['balFlag']:
                fracDiffs['BAL'].append(fracDiff)
            else:
                fracDiffs['nonBAL'].append(fracDiff)
    loss = np.mean(lossSquared)
    print("Reconstruction Loss is: {0}".format(loss))
    for key in fracDiffs.keys():
        pass # plot_frac_diffs(wave, fracDiffs[key], name=key, saveDir=saveDir)

    plot_dict_of_frac_diffs(wave, fracDiffs, title='all\_bal\_nonBal', saveDir=saveDir, loss="loss={}".format(round(loss, 5)))

    return loss


def frac_diff_clusters(wave, clusters, saveDir):
    clusterNames = list(clusters.keys())
    fracDiffs = dict((key, []) for key in clusterNames)

    for clusterName in clusterNames:
        sdssFluxes = clusters[clusterName]['sdssFluxes']
        reconFluxes = clusters[clusterName]['reconFluxes']
        for i in range(len(sdssFluxes)):
            fracDiff = frac_diff(sdssFluxes[i], reconFluxes[i])
            fracDiffs[clusterName].append(fracDiff)

    for key in fracDiffs.keys():
        pass # plot_frac_diffs(wave, fracDiffs[key], name='Cluster_%s' % key, saveDir)

    plot_dict_of_frac_diffs(wave, fracDiffs, title='Clusters', saveDir=saveDir)


def plot_frac_diffs(wave, fracDiffs, name, saveDir):
    medianFracDiffs = np.median(fracDiffs, axis=0)
    plt.figure()
    plt.plot(wave, medianFracDiffs, label='Median')
    # plt.plot(np.mean(fracDiffs, axis=0), label='Mean')
    plt.axhline(0, color='k')
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel("Fractional Difference")
    plt.legend()
    plt.savefig("{0}/Fractional_Difference_{1}.png".format(saveDir, name).replace('\\', ''))


def plot_dict_of_frac_diffs(wave, fracDiffsDict, title, saveDir, loss=''):
    plt.figure()
    plt.title("{0} {1}".format(title, loss))
    for key, fracDiffs in fracDiffsDict.items():
        medianFracDiffs = np.median(fracDiffs, axis=0)
        plt.plot(wave, medianFracDiffs, label=key)
    plt.axhline(0, color='k')
    # reconFlux = 0.03 * (reconFlux - min(reconFlux)) / (max(reconFlux) - min(reconFlux))
    # plt.plot(reconFlux)
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel("Median Fractional Difference")
    plt.legend()
    plt.savefig("{0}/Fractional_Difference_{1}.png".format(saveDir, title).replace('\\', ''))
