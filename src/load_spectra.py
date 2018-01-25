import pickle
import json
from src.component_reconstruction import get_components


def load_spectra(savedSpectra='data_files/created/spectra2.pickle', savedProperties='data_files/created/balProperties.json'):
    with open(savedSpectra, 'rb') as f:
        spectra = pickle.load(f)

    comps = get_components(componentsFile='data_files/given/dm_6c_16003000_171024.comp')

    if savedProperties is not None:
        nameNotInBalProps = []
        with open(savedProperties, 'r') as f:
            balProperties = json.load(f)
        for name in spectra.keys():
            try:
                spectra[name]['snr'] = balProperties[name]['snr']
            except KeyError:
                print(name)
                nameNotInBalProps.append(name)

        return spectra, comps, balProperties, nameNotInBalProps
    else:
        return spectra, comps
