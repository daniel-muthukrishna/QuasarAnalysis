import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def get_mags(spectra, removeNames=[]):
    filters = ['u', 'g', 'r', 'i', 'z']
    magsArray, magsErrArray, weightsArray = [], [], []
    names = list(spectra.keys())
    numComps = len(spectra[names[0]]['weights'])
    weightLabels = ['w{0}'.format(i+1) for i in range(numComps)]
    for name in names:
        if name not in removeNames:
            bal = 'BAL' if spectra[name]['balFlag'] else 'non-BAL'
            if spectra[name]['mags'][0] != 0:
                magsArray.append(dict(zip(filters, spectra[name]['mags'])))
                magsErrArray.append(dict(zip(filters, spectra[name]['magsErr'][0])))
                weightsArray.append(spectra[name]['weights'])

    magsDF = pd.DataFrame(magsArray)
    magsErrDF = pd.DataFrame(magsErrArray)
    weightsDF = pd.DataFrame(data=np.array(weightsArray), columns=weightLabels)

    return magsDF, magsErrDF, weightsDF


def plot_mags(spectra, removeNames=[], saveDir=''):
    magDF, magErrDF, weightsDF = get_mags(spectra, removeNames)

    fig, ax = plt.subplots(nrows=len(magDF.columns), ncols=len(weightsDF.columns), sharex='col', sharey='row', figsize=(18,15))
    fig.subplots_adjust(wspace=0, hspace=0)
    for i, f in enumerate(magDF.columns):
        for j, w in enumerate(weightsDF.columns):
            ax[i, j].scatter(weightsDF[w], magDF[f], alpha=0.3)
            if i == len(magDF.columns) - 1:
                ax[i, j].set_xlabel(w)
            if j == 0:
                ax[i, j].set_ylabel(f)
    plt.savefig('{0}/magnitudes_vs_weights'.format(saveDir), bbox_inches='tight')

    # for f in magDF.columns:
    #     for w in weightsDF.columns:
    #         title = f + '_mag vs ' + w
    #         plt.figure(title)
    #         plt.title(title)
    #         plt.xlabel(w)
    #         plt.ylabel(f)
    #         plt.scatter(weightsDF[w], magDF[f], alpha=0.4)
    #         # plt.errorbar(weightsDF[w], magDF[f], yerr=magErrDF[f], fmt='o')
    #         # plt.savefig("%s/title" % saveDir)

