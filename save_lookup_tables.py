import pickle
import json
import numpy as np
from component_reconstruction import spectra_dict


def save_spectra(componentsFile, waveFile, weightsFile):
    spectra = spectra_dict(componentsFile, waveFile, weightsFile)
    with open('spectra.pickle', 'wb') as f:
        pickle.dump(spectra, f, pickle.HIGHEST_PROTOCOL)


def save_filepaths(filepathsFile):
    filepathsDict = {}
    filepathsInfo = np.genfromtxt(filepathsFile, dtype=str)
    for name, icaFile, sdssFile in filepathsInfo:
        filepathsDict[name] = {'ica': icaFile, 'sdss': sdssFile}

    with open('spectraFilepaths.json', 'w') as f:
        json.dump(filepathsDict, f, sort_keys=True)

    return filepathsDict


if __name__ == '__main__':
    save_filepaths('dm_hbal_files.dat')
    save_spectra(componentsFile='dm_6c_16003000_171024.comp', waveFile='wav_16003000.dat', weightsFile='dm_hbal_weights.dat')