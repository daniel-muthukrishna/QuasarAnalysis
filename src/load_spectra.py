import pickle


def load_spectra(savedSpectra='data_files/created/spectra.pickle'):
    with open(savedSpectra, 'rb') as f:
        spectra = pickle.load(f)

    return spectra
