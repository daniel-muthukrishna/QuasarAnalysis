import os
import matplotlib.pyplot as plt
from src.load_spectra import load_spectra
from src.check_reconstructions import compare_all, frac_diff_all, frac_diff_clusters
from src.analyse_weights import clustering, analyse_clusters, plot_weights, plot_clusters
from src.pca_analysis import pca_analysis, ica_analysis
from src.component_reconstruction import plot_each_component
from src.analyse_BAL_properties import plot_mags


class AnalyseSpectraComponents(object):
    def __init__(self, spectra, saveDir, title=''):
        self.spectra = spectra
        self.names = list(spectra.keys())
        self.wave = spectra[self.names[0]]['reconWave']
        self.removeNames = self.remove_names()
        self.saveDir = saveDir
        self.title = title
        if not os.path.exists(saveDir):
            os.makedirs(saveDir)
        print(title)

    def remove_names(self):
        removeNames = compare_all(self.spectra)
        return removeNames

    def get_clusters(self):
        # clustersk2 = clustering(self.spectra, self.removeNames, numClusters=2, plotClusters=True, saveDir=self.saveDir)
        clustersk3 = clustering(self.spectra, self.removeNames, numClusters=3, plotClusters=True, saveDir=self.saveDir)
        # clustersk4 = clustering(self.spectra, self.removeNames, numClusters=4, plotClusters=True, saveDir=self.saveDir)
        # clustersk5 = clustering(self.spectra, self.removeNames, numClusters=5, plotClusters=True, saveDir=self.saveDir)
        analyse_clusters(clustersk3, self.wave, saveDir=self.saveDir)
        return clustersk3

    def plot_weights(self):
        plot_weights(self.spectra, self.removeNames, saveDir=self.saveDir)

    def calc_frac_diff(self, clusters=None):
        loss = frac_diff_all(self.wave, self.spectra, self.removeNames, saveDir=self.saveDir, title=self.title)
        if clusters is not None:
            frac_diff_clusters(self.wave, clusters, saveDir=self.saveDir)

    def do_all(self):
        clusters = self.get_clusters()
        self.plot_weights()
        self.calc_frac_diff(clusters=clusters)


def main():
    # ORIGINAL COMPS ANALYSIS
    spectra, comps = load_spectra()
    compsAnalysis = AnalyseSpectraComponents(spectra, saveDir='Figures/OriginalComps', title='Original 6 comps')
    compsAnalysis.do_all()
    # plot_mags(spectra, saveDir='Figures/OriginalComps')
    # plot_each_component(spectra, comps)

    # MY 5-Comp PCA COMPS ANALYSIS
    spectraPCA, comps, pcaMean = pca_analysis(spectra, nComps=5)
    myPcaCompsAnalysis = AnalyseSpectraComponents(spectraPCA, saveDir='Figures/MyPCA', title='PCA 5 comps')
    myPcaCompsAnalysis.do_all()
    # plot_each_component(spectraPCA, comps, pcaMean)

    # MY 5-Comp ICA COMPS ANALYSIS
    spectraICA, comps, icaMean = ica_analysis(spectra, nComps=5)
    myIcaCompsAnalysis = AnalyseSpectraComponents(spectraICA, saveDir='Figures/MyICA', title='Fast-ICA 5 comps')
    myIcaCompsAnalysis.do_all()
    # plot_each_component(spectraICA, comps, icaMean)

    plt.show()


if __name__ == '__main__':
    main()
