import matplotlib.pyplot as plt
from src.load_spectra import load_spectra
from src.check_reconstructions import compare_all, frac_diff_all, frac_diff_clusters
from src.analyse_weights import clustering, analyse_clusters, plot_weights, plot_clusters
from src.pca_analysis import pca_analysis


def main():
    spectra = load_spectra()
    # names = list(spectra.keys())
    # wave = spectra[names[0]]['reconWave']
    # removeNames = compare_all(spectra)
    # clusters = clustering(spectra, removeNames, plotClusters=False)
    # plot_weights(spectra, removeNames)
    # analyse_clusters(clusters, wave)
    # frac_diff_all(spectra)
    # frac_diff_clusters(clusters)
    pca_analysis(spectra)

    # for idx in range(10):
    #     plot_spectrum(names1[idx], spectra1)

    plt.show()


if __name__ == '__main__':
    main()
