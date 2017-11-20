import pickle
from src.component_reconstruction import get_components

def load_spectra(savedSpectra='data_files/created/spectra.pickle'):
    with open(savedSpectra, 'rb') as f:
        spectra = pickle.load(f)
    comps = get_components(componentsFile='data_files/given/dm_6c_16003000_171024.comp')

    return spectra, comps
