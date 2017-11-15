import pickle
import matplotlib.pyplot as plt
from src.remove_bad_reconstructions import compare_all
from src.analyse_weights import clustering, analyse_clusters, plot_weights


def load_spectra(savedSpectra='data_files/created/spectra.pickle'):
    with open(savedSpectra, 'rb') as f:
        spectra = pickle.load(f)

    return spectra


def main():
    spectra = load_spectra()
    names = list(spectra.keys())
    wave = spectra[names[0]]['reconWave']
    removeNames = compare_all(spectra)
    plot_weights(spectra, removeNames)
    clusters = clustering(spectra, removeNames)
    analyse_clusters(clusters, wave)
    # for idx in range(10):
    #     plot_spectrum(names1[idx], spectra1)
    plt.show()


if __name__ == '__main__':
    main()
