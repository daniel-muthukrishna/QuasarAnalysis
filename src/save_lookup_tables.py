import json
import pickle
import os
import numpy as np
from src.component_reconstruction import spectra_dict

MAIN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def save_spectra(componentsFile, waveFile, weightsFile):
    spectra = spectra_dict(componentsFile, waveFile, weightsFile)
    with open(os.path.join(MAIN_DIR, 'data_files/created/spectra.pickle'), 'wb') as f:
        pickle.dump(spectra, f, pickle.HIGHEST_PROTOCOL)


def save_filepaths(filepathsFile):
    filepathsDict = {}
    filepathsInfo = np.genfromtxt(filepathsFile, dtype=str)
    for name, icaFile, sdssFile in filepathsInfo:
        filepathsDict[name] = {'ica': icaFile, 'sdss': sdssFile}

    with open(os.path.join(MAIN_DIR, 'data_files/created/spectraFilepaths.json'), 'w') as f:
        json.dump(filepathsDict, f, sort_keys=True)

    return filepathsDict


if __name__ == '__main__':
    save_filepaths(os.path.join(MAIN_DIR, 'data_files/given/dm_hbal_files.dat'))
    save_spectra(componentsFile=os.path.join(MAIN_DIR, 'data_files/given/dm_6c_16003000_171024.comp'),
                 waveFile=os.path.join(MAIN_DIR, 'data_files/given/wav_16003000.dat'),
                 weightsFile=os.path.join(MAIN_DIR, 'data_files/given/dm_hbal_weights.dat'))
