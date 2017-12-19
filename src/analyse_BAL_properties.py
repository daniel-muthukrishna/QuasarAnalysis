import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def get_mags(spectra, removeNames=[], balProps=None):
    filters = ['u', 'g', 'r', 'i', 'z']
    magsArray, magsErrArray, weightsArray, otherParams = [], [], [], []
    names = list(spectra.keys())
    numComps = len(spectra[names[0]]['weights'])
    weightLabels = ['w{0}'.format(i+1) for i in range(numComps)]
    for name in names:
        if name not in removeNames:  # name not in DR12Q.fits
            bal = 'BAL' if spectra[name]['balFlag'] else 'non-BAL'
            if spectra[name]['mags'][0] != 0:
                magsArray.append(dict(zip(filters, spectra[name]['mags'])))
                magsErrArray.append(dict(zip(filters, spectra[name]['magsErr'])))
                weightsArray.append(spectra[name]['weights'])
                otherParams.append({'CIV_EW': balProps[name]['ewCIV'], 'CIV_EW_Err': balProps[name]['ewCIVErr'],
                                    'CIII_EW': balProps[name]['ewCIII'], 'CIII_EW_Err': balProps[name]['ewCIIIErr'],
                                    'z_CIV': balProps[name]['zCIV'], 'z_CIII': balProps[name]['zCIII'],
                                    'z_PCA': balProps[name]['zPCA']})

    magsDF = pd.DataFrame(magsArray)
    magsErrDF = pd.DataFrame(magsErrArray)
    weightsDF = pd.DataFrame(data=np.array(weightsArray), columns=weightLabels)
    otherParams = pd.DataFrame(otherParams)

    return magsDF, magsErrDF, weightsDF, otherParams


class PlotProperties(object):
    def __init__(self, spectra, balProps, removeNames=[], saveDir=''):
        self.magDF, self.magErrDF, self.weightsDF, self.otherParams = get_mags(spectra, removeNames, balProps)
        self.saveDir = saveDir

    def mags_vs_comps(self):
        fig, ax = plt.subplots(nrows=len(self.magDF.columns), ncols=len(self.weightsDF.columns), sharex='col', sharey='row',
                               figsize=(18, 15))
        fig.subplots_adjust(wspace=0, hspace=0)
        for i, f in enumerate(self.magDF.columns):
            for j, w in enumerate(self.weightsDF.columns):
                ax[i, j].scatter(self.weightsDF[w], self.magDF[f], alpha=0.3)
                if i == len(self.magDF.columns) - 1:
                    ax[i, j].set_xlabel(w)
                if j == 0:
                    ax[i, j].set_ylabel(f)
        plt.savefig('{0}/magnitudes_vs_weights'.format(self.saveDir), bbox_inches='tight')

        # for f in magDF.columns:
        #     for w in weightsDF.columns:
        #         title = f + '_mag vs ' + w
        #         plt.figure(title)
        #         plt.title(title)
        #         plt.xlabel(w)
        #         plt.ylabel(f)
        #         plt.scatter(self.weightsDF[w], sel.fmagDF[f], alpha=0.4)
        #         # plt.errorbar(self.weightsDF[w], self.magDF[f], yerr=self.magErrDF[f], fmt='o')
        #         # plt.savefig("%s/title" % self.saveDir)


    def civ_vs_mags(self):
        pass


