import json
import pickle
import os
import numpy as np
from astropy.io import fits
from src.component_reconstruction import spectra_dict

MAIN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def save_spectra(componentsFile, waveFile, weightsFile):
    spectra = spectra_dict(componentsFile, waveFile, weightsFile)
    with open(os.path.join(MAIN_DIR, 'data_files/created/spectra.pickle'), 'wb') as f:
        pickle.dump(spectra, f, pickle.HIGHEST_PROTOCOL)


def save_filepaths(filePathsFile, savePath='data_files/created/spectraFilepaths.json'):
    filepathsDict = {}
    filepathsInfo = np.genfromtxt(filePathsFile, dtype=str)
    for name, icaFile, sdssFile in filepathsInfo:
        filepathsDict[name] = {'ica': icaFile, 'sdss': sdssFile}

    savePathDir = os.path.join(MAIN_DIR, os.path.dirname(savePath))
    if not os.path.exists(savePathDir):
        os.makedirs(savePathDir)
    with open(os.path.join(MAIN_DIR, savePath), 'w') as f:
        json.dump(filepathsDict, f, sort_keys=True)

    return filepathsDict


def save_BAL_properties(filePath, savePath='data_files/created/balProperties.json'):
    balProperties = {}
    with fits.open(filePath, memmap=False) as hdulist:
        data = hdulist[1].data
        names = data['SDSS_NAME']
        snr = data['SNR_SPEC']
        zCIV = data['Z_CIV']
        zCIII = data['Z_CIII']
        zPCA = data['Z_PCA']
        ewCIV = data['REWE_CIV']
        ewCIVErr = data['ERR_REWE_CIV']
        ewCIII = data['REWE_CIII']
        ewCIIIErr = data['ERR_REWE_CIII']
        fwhmCIII = data['FWHM_CIII']

    for i in range(len(names)):
        balProperties['SDSSJ'+names[i]] = {'snr': snr[i], 'ewCIV': ewCIV[i], 'ewCIVErr': ewCIVErr[i], 'ewCIII': ewCIII[i],
                                   'ewCIIIErr': ewCIIIErr[i], 'zCIV': zCIV[i], 'zCIII': zCIII[i], 'zPCA': zPCA[i]}

    with open(os.path.join(MAIN_DIR, savePath), 'w') as f:
        json.dump(balProperties, f, sort_keys=True)


if __name__ == '__main__':
    # save_filepaths(os.path.join(MAIN_DIR, 'data_files/given/dm_hbal_files.dat'))
    # save_spectra(componentsFile=os.path.join(MAIN_DIR, 'data_files/given/dm_6c_16003000_171024.comp'),
    #              waveFile=os.path.join(MAIN_DIR, 'data_files/given/wav_16003000.dat'),
    #              weightsFile=os.path.join(MAIN_DIR, 'data_files/given/dm_hbal_weights.dat'))
    save_BAL_properties(os.path.join(MAIN_DIR, 'data_files/given/DR12Q.fits'))
