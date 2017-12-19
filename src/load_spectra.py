import pickle
import json
from src.component_reconstruction import get_components


def load_spectra(savedSpectra='data_files/created/spectra.pickle', savedProperties='data_files/created/balProperties.json'):
    with open(savedSpectra, 'rb') as f:
        spectra = pickle.load(f)

    # if savedProperties is not None:
    #     with open(savedProperties, 'r') as f:
    #         balProperties = json.load(f)
    #     for name in spectra.keys():
    #         try:
    #             spectra[name]['snr'] = balProperties[name.replace('SDSSJ', '')]['snr']
    #         except KeyError:
    #             print(name)

    comps = get_components(componentsFile='data_files/given/dm_6c_16003000_171024.comp')

    return spectra, comps
