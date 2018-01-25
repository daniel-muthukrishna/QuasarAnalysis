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


def dispersion(data, recon, noise):
    dispersionStat = ((data - recon)/noise)**2
    return dispersionStat


def mask_spectrum(flux, fluxErr, recon):
    # Set all flux zeros to nan
    maskIndexes = np.where((fluxErr <= 0))[0]
    maskedFlux = np.copy(flux)
    maskedFlux[maskIndexes] = np.nan

    # maskIndexes1 = np.where((flux <= 0))[0]
    # maskedFlux1 = np.copy(flux)
    # maskedFlux1[maskIndexes1] = recon[maskIndexes1]
    # fluxErr[maskIndexes1] = 0.05
    # flux = np.copy(maskedFlux1)
    # maskIndexes = np.where(fluxErr < 0)[0]
    # maskedFlux = np.copy(flux)
    # maskedFlux[maskIndexes] = np.nan

    return maskedFlux, maskIndexes


def dispersion_all(wave, spectra, removeNames, saveDir, title=''):
    dispersionStats = {'All': [], 'BAL': [], 'nonBAL': []}
    lossSquared = {'All': [], 'BAL': [], 'nonBAL': []}
    chi2 = {'All': [], 'BAL': [], 'nonBAL': []}
    names = list(spectra.keys())
    # color = iter(plt.cm.rainbow(np.linspace(0, 1, 18)))
    i=0
    count = 0
    for name in names:
        if name not in removeNames:
            sdssFlux = spectra[name]['sdssFlux']
            reconFlux = spectra[name]['reconFlux']
            noise = spectra[name]['sdssFluxErr']
            sdssFlux, maskIndexes = mask_spectrum(sdssFlux, noise, reconFlux)

            dispersionStat = dispersion(sdssFlux, reconFlux, noise)
            dispersionStats['All'].append(dispersionStat)
            lossSquared['All'].append((sdssFlux-reconFlux)**2)
            chi2['All'].append(chisquare(sdssFlux, reconFlux)[0])
            if spectra[name]['balFlag']:
                dispersionStats['BAL'].append(dispersionStat)
                lossSquared['BAL'].append((sdssFlux - reconFlux) ** 2)
                chi2['BAL'].append(chisquare(sdssFlux, reconFlux)[0])
            else:
                dispersionStats['nonBAL'].append(dispersionStat)
                lossSquared['nonBAL'].append((sdssFlux - reconFlux) ** 2)
                chi2['nonBAL'].append(chisquare(sdssFlux, reconFlux)[0])
            if np.isnan(np.min(sdssFlux)):
                count+=1
            # print(name, wave[maskIndexes], sdssFlux[maskIndexes], spectra[name]['sdssFlux'][maskIndexes], noise[maskIndexes])
            # plt.figure()
            # plt.plot(wave, dispersionStat, label=r"$((spec - recon)/specNoise)^2$".format(name, spectra[name]['balFlag']))
            # plt.plot(wave, sdssFlux, label='spec')
            # plt.plot(wave, reconFlux, color='k', label='recon')
            # plt.plot(wave, noise, label='specNoise')
            # plt.xlabel('Wavelength')
            # plt.legend()
            # plt.show()
    print("Count is: ", count)

    loss = {key: np.median(val) for key, val in lossSquared.items()}
    lossErr = {key: np.std(val) for key, val in lossSquared.items()}
    chi2Median = {key: np.median(val) for key, val in chi2.items()}
    chi2Err = {key: np.std(val) for key, val in chi2.items()}
    print("Reconstruction Loss is: {0}".format(loss))
    print("Chi2 is: {0}".format(chi2Median))
    for key in dispersionStats.keys():
        pass  # plot_dispersions(wave, dispersionStats[key], name=key, saveDir=saveDir)

    plot_dict_of_dispersions(wave, dispersionStats, title=title, saveDir=saveDir, loss=loss, chi2=chi2Median, lossErr=lossErr, chi2Err=chi2Err)

    return loss


def dispersion_clusters(wave, clusters, saveDir):
    clusterNames = list(clusters.keys())
    dispersionStats = dict((key, []) for key in clusterNames)

    for clusterName in clusterNames:
        sdssFluxes = clusters[clusterName]['sdssFluxes']
        reconFluxes = clusters[clusterName]['reconFluxes']
        noise = clusters[clusterName]['sdssFluxesErr']
        for i in range(len(sdssFluxes)):
            sdssFlux = mask_spectrum(sdssFluxes[i], noise[i], reconFluxes[i])
            dispersionStat = dispersion(sdssFlux, reconFluxes[i], noise[i])
            dispersionStats[clusterName].append(dispersionStat)

    for key in dispersionStats.keys():
        pass  # plot_dispersions(wave, dispersionStats[key], name='Cluster_%s' % key, saveDir)

    plot_dict_of_dispersions(wave, dispersionStats, title='Clusters', saveDir=saveDir)


def plot_dispersions(wave, dispersionStats, name, saveDir):
    mediandispersionStats = np.median(dispersionStats, axis=0)
    plt.figure()
    plt.plot(wave, mediandispersionStats, label='Median')
    # plt.plot(np.mean(dispersionStats, axis=0), label='Mean')
    plt.axhline(0, color='k')
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel("Fractional Difference")
    plt.legend()
    plt.savefig("{0}/Fractional_Difference_{1}.png".format(saveDir, name).replace('\\', ''))


def plot_dict_of_dispersions(wave, dispersionStatsDict, title, saveDir, loss=None, chi2=None, lossErr=None, chi2Err=None):
    fig, ax = plt.subplots(1, sharex=True)
    for key, dispersionStats in dispersionStatsDict.items():
        lossVal = ": loss={0}".format(round(loss[key], 7)) if loss is not None else ''
        lossValErr = "$\pm${0}".format(round(lossErr[key], 1)) if lossErr is not None else ''
        chi2Val = "_chi2={0}".format(round(chi2[key], 1)) if chi2 is not None else ''
        chi2ValErr = "$\pm${0}".format(round(chi2Err[key], 0)) if chi2Err is not None else ''
        mediandispersionStats = np.median(dispersionStats, axis=0)
        stddispersionStats = np.std(dispersionStats, axis=0)
        dispersionStats = np.array(dispersionStats)
        rms = np.sqrt(np.nanmean(dispersionStats, axis=0))
        ax.plot(wave, rms, label="{0}".format(key), zorder=10)
        # ax[1].plot(wave, mediandispersionStats, zorder=10)
        # plt.fill_between(wave, mediandispersionStats-stddispersionStats, mediandispersionStats+stddispersionStats, alpha=0.3)
    ax.set_title(title)
    ax.set_ylabel(r"$\sqrt{\frac{\sum{ ((spec - recon)/specNoise)^2}}{N}}$")
    # ax.axhline(0, color='k')
    ax.grid()
    ax.legend()
    # ax[1].set_ylabel(r"Median of all $((data - recon)/snr)^2$")
    # ax[1].grid()
    plt.xlabel('Wavelength ($\AA$)')
    plt.savefig("{0}/new_statistic_abs_{1}.png".format(saveDir, title).replace('\\', ''))
