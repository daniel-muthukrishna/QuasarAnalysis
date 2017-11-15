import json
import pickle
import numpy as np
from src.component_reconstruction import spectra_dict


def save_spectra(componentsFile, waveFile, weightsFile):
    spectra = spectra_dict(componentsFile, waveFile, weightsFile)
    with open('data_files/created/spectra.pickle', 'wb') as f:
        pickle.dump(spectra, f, pickle.HIGHEST_PROTOCOL)


def save_filepaths(filepathsFile):
    filepathsDict = {}
    filepathsInfo = np.genfromtxt(filepathsFile, dtype=str)
    for name, icaFile, sdssFile in filepathsInfo:
        filepathsDict[name] = {'ica': icaFile, 'sdss': sdssFile}

    with open('data_files/created/spectraFilepaths.json', 'w') as f:
        json.dump(filepathsDict, f, sort_keys=True)

    return filepathsDict


if __name__ == '__main__':
    save_filepaths('data_files/given/dm_hbal_files.dat')
    save_spectra(componentsFile='data_files/given/dm_6c_16003000_171024.comp',
                 waveFile='data_files/given/wav_16003000.dat',
                 weightsFile='data_files/given/dm_hbal_weights.dat')
