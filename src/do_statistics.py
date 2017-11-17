import matplotlib.pyplot as plt
from src.load_spectra import load_spectra
from src.check_reconstructions import compare_all, frac_diff_all, frac_diff_clusters
from src.analyse_weights import clustering, analyse_clusters, plot_weights, plot_clusters
from src.pca_analysis import pca_analysis
from src.component_reconstruction import plot_each_component


def main():
    spectra = load_spectra()
    names = list(spectra.keys())
    wave = spectra[names[0]]['reconWave']
    removeNames = compare_all(spectra)
    clusters = clustering(spectra, removeNames, plotClusters=False)
    plot_weights(spectra, removeNames)
    analyse_clusters(clusters, wave)
    frac_diff_all(wave, spectra)
    frac_diff_clusters(wave, clusters)
    pca_analysis(spectra)
    plot_each_component(spectra, componentsFile='data_files/given/dm_6c_16003000_171024.comp')

    # for idx in range(10):
    #     plot_spectrum(names1[idx], spectra1)

    plt.show()


if __name__ == '__main__':
    main()
