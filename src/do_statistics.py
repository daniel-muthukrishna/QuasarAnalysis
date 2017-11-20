import os
import matplotlib.pyplot as plt
from src.load_spectra import load_spectra
from src.check_reconstructions import compare_all, frac_diff_all, frac_diff_clusters
from src.analyse_weights import clustering, analyse_clusters, plot_weights, plot_clusters
from src.pca_analysis import pca_analysis, ica_analysis
from src.component_reconstruction import plot_each_component


class AnalyseSpectraComponents(object):
    def __init__(self, spectra, saveDir):
        self.spectra = spectra
        self.names = list(spectra.keys())
        self.wave = spectra[self.names[0]]['reconWave']
        self.removeNames = self.remove_names()
        self.saveDir = saveDir
        if not os.path.exists(saveDir):
            os.makedirs(saveDir)

    def remove_names(self):
        removeNames = compare_all(self.spectra)
        return removeNames

    def get_clusters(self):
        clusters = clustering(self.spectra, self.removeNames, plotClusters=False, saveDir=self.saveDir)
        analyse_clusters(clusters, self.wave, saveDir=self.saveDir)
        return clusters

    def plot_weights(self):
        plot_weights(self.spectra, self.removeNames, saveDir=self.saveDir)

    def calc_frac_diff(self, clusters=None):
        frac_diff_all(self.wave, self.spectra, self.removeNames, saveDir=self.saveDir)
        if clusters is not None:
            frac_diff_clusters(self.wave, clusters, saveDir=self.saveDir)

    def do_all(self):
        clusters = self.get_clusters()
        self.plot_weights()
        self.calc_frac_diff(clusters)


def main():
    # ORIGINAL COMPS ANALYSIS
    spectra, comps = load_spectra()
    compsAnalysis = AnalyseSpectraComponents(spectra, saveDir='Figures/OriginalComps')
    compsAnalysis.do_all()
    # plot_each_component(spectra, comps)

    # MY PCA COMPS ANALYSIS
    spectraPCA, comps = pca_analysis(spectra)
    myPcaCompsAnalysis = AnalyseSpectraComponents(spectraPCA, saveDir='Figures/MyPCA')
    myPcaCompsAnalysis.do_all()
    # plot_each_component(spectra, comps)

    # MY ICA COMPS ANALYSIS
    spectraICA, comps = ica_analysis(spectra)
    myIcaCompsAnalysis = AnalyseSpectraComponents(spectraICA, saveDir='Figures/MyICA')
    myIcaCompsAnalysis.do_all()
    # plot_each_component(spectra, comps)

    plt.show()


if __name__ == '__main__':
    main()
