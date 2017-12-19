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


def frac_diff(data, recon, noise):
    fracDiff = ((data - recon)/noise)**2
    return fracDiff


def frac_diff_all(wave, spectra, removeNames, saveDir, title=''):
    fracDiffs = {'All': [], 'BAL': [], 'nonBAL': []}
    lossSquared = {'All': [], 'BAL': [], 'nonBAL': []}
    chi2 = {'All': [], 'BAL': [], 'nonBAL': []}
    names = list(spectra.keys())
    # color = iter(plt.cm.rainbow(np.linspace(0, 1, 18)))
    plt.figure()
    for name in names:
        if name not in removeNames:
            sdssFlux = spectra[name]['sdssFlux']
            reconFlux = spectra[name]['reconFlux']
            noise = spectra[name]['sdssFluxErr']

            fracDiff = frac_diff(sdssFlux, reconFlux, noise)
            fracDiffs['All'].append(fracDiff)
            lossSquared['All'].append((sdssFlux-reconFlux)**2)
            chi2['All'].append(chisquare(sdssFlux, reconFlux)[0])
            if spectra[name]['balFlag']:
                fracDiffs['BAL'].append(fracDiff)
                lossSquared['BAL'].append((sdssFlux - reconFlux) ** 2)
                chi2['BAL'].append(chisquare(sdssFlux, reconFlux)[0])
            else:
                fracDiffs['nonBAL'].append(fracDiff)
                lossSquared['nonBAL'].append((sdssFlux - reconFlux) ** 2)
                chi2['nonBAL'].append(chisquare(sdssFlux, reconFlux)[0])
            # c=next(color)
            # plt.plot(wave, fracDiff, label='{0}_BAL={1}'.format(name, spectra[name]['balFlag']), color=c)
            # plt.plot(wave, sdssFlux, color=c)
            # plt.plot(wave, reconFlux, color='k')
    loss = {key: np.median(val) for key, val in lossSquared.items()}
    lossErr = {key: np.std(val) for key, val in lossSquared.items()}
    chi2Median = {key: np.median(val) for key, val in chi2.items()}
    chi2Err = {key: np.std(val) for key, val in chi2.items()}
    print("Reconstruction Loss is: {0}".format(loss))
    print("Chi2 is: {0}".format(chi2Median))
    for key in fracDiffs.keys():
        pass  # plot_frac_diffs(wave, fracDiffs[key], name=key, saveDir=saveDir)

    plot_dict_of_frac_diffs(wave, fracDiffs, title=title, saveDir=saveDir, loss=loss, chi2=chi2Median, lossErr=lossErr, chi2Err=chi2Err)

    return loss


def frac_diff_clusters(wave, clusters, saveDir):
    clusterNames = list(clusters.keys())
    fracDiffs = dict((key, []) for key in clusterNames)

    for clusterName in clusterNames:
        sdssFluxes = clusters[clusterName]['sdssFluxes']
        reconFluxes = clusters[clusterName]['reconFluxes']
        noise = clusters[clusterName]['sdssFluxesErr']
        for i in range(len(sdssFluxes)):
            fracDiff = frac_diff(sdssFluxes[i], reconFluxes[i], noise[i])
            fracDiffs[clusterName].append(fracDiff)

    for key in fracDiffs.keys():
        pass  # plot_frac_diffs(wave, fracDiffs[key], name='Cluster_%s' % key, saveDir)

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


def plot_dict_of_frac_diffs(wave, fracDiffsDict, title, saveDir, loss=None, chi2=None, lossErr=None, chi2Err=None):
    fig, ax = plt.subplots(1, sharex=True)
    for key, fracDiffs in fracDiffsDict.items():
        lossVal = ": loss={0}".format(round(loss[key], 7)) if loss is not None else ''
        lossValErr = "$\pm${0}".format(round(lossErr[key], 1)) if lossErr is not None else ''
        chi2Val = "_chi2={0}".format(round(chi2[key], 1)) if chi2 is not None else ''
        chi2ValErr = "$\pm${0}".format(round(chi2Err[key], 0)) if chi2Err is not None else ''
        medianFracDiffs = np.median(fracDiffs, axis=0)
        stdFracDiffs = np.std(fracDiffs, axis=0)
        rms = np.sqrt(np.median(fracDiffs, axis=0))
        ax.plot(wave, rms, label="{0}".format(key), zorder=10)
        # ax[1].plot(wave, medianFracDiffs, zorder=10)
        # plt.fill_between(wave, medianFracDiffs-stdFracDiffs, medianFracDiffs+stdFracDiffs, alpha=0.3)
    ax.set_title(title)
    ax.set_ylabel(r"$\sqrt{Median \ ((spec - recon)/specNoise)^2}$")
    # ax.axhline(0, color='k')
    ax.grid()
    ax.legend()
    # ax[1].set_ylabel(r"Median of all $((data - recon)/snr)^2$")
    # ax[1].grid()
    plt.xlabel('Wavelength ($\AA$)')
    plt.savefig("{0}/new_statistic{1}.png".format(saveDir, title).replace('\\', ''))
