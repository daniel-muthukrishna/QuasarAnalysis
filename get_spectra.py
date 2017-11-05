# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 15:18:10 2015

@author: lc585
Edited 15-Dec-2016 by mjt91 to include a mask array for DR12

"""
from __future__ import division

from astropy.table import Table
import numpy as np
import json
import warnings
from astropy.utils.exceptions import AstropyWarning
from astropy.io import fits
from get_wavelength import get_wavelength
import os


def get_sdss_dr7_spec(name):
    """
    Given SDSS DR7 name will return wavelength, dw, flux and err.
    Uses Paul's better spectra

    """

    # Open master list
    f = open('/data/vault/phewett/DR7/ordn_splist_dr7_master.lis', 'r')

    SpecFiles, SDSSNames = [], []
    # Loop over lines and extract variables of interest
    for line in f:
        line = line.strip()  # removes '\n' from eol
        columns = line.split()
        SpecFiles.append(columns[0])
        SDSSNames.append(columns[1])

    f.close()

    SpecFiles = np.array(SpecFiles)
    SDSSNames = np.array(SDSSNames)

    if name[:5] != 'SDSSJ':
        name = 'SDSSJ' + name

    i = np.where(SDSSNames == name)[0]

    if len(i) != 0:
        hdulist = fits.open(SpecFiles[i][0])
        data = hdulist[0].data
        hdr = hdulist[0].header
        hdulist.close()

        wavelength, dw = get_wavelength(hdr)
        flux = data[0, :].flatten()
        err = data[1, :].flatten()

        # create the mask array
        mask = []
        for j in range(len(wavelength)):
            if 5574 < wavelength[j] < 5583:
                mask.append(0)
            elif 6297 < wavelength[j] < 6306:
                mask.append(0)
            else:
                mask.append(1)

        return wavelength, dw, flux, err, mask

    else:
        print('No spectrum found')
        return None, None, None, None, None


def get_boss_dr12_spec(name):
    """
    Given DR12 name will return spectrum
    """

    t = Table.read('/data/vault/phewett/ICAtest/DR12exp/Spectra/boss_dr12_name_file.lis',
                   format='ascii',
                   names=['name', 'loc'])

    if name[:5] != 'SDSSJ':
        name = 'SDSSJ' + name

    hdulist = fits.open(t[t['name'] == name]['loc'].data[0])

    data = hdulist[1].data
    hdr = hdulist[0].header
    hdulist.close()

    wavelength = 10 ** np.array([j[1] for j in data])
    dw = wavelength * (10 ** hdr['COEFF1'] - 1.0)
    flux = np.array([j[0] for j in data])
    err = np.sqrt(abs(np.array(
        [j[6] for j in data])))  # square root of sky #abs added by mjt91 in an attempt to remove error messages

    # create the mask array
    mask = []
    for j in range(len(wavelength)):
        if 5574 < wavelength[j] < 5583:
            mask.append(0)
        elif 6297 < wavelength[j] < 6306:
            mask.append(0)
        else:
            mask.append(1)

    return wavelength, dw, flux, err, mask


def load_spectra_filenames(spectraFile='spectraFilenames.json'):
    with open(spectraFile, 'r') as f:
        filenamesDict = json.load(f)

    return filenamesDict


def get_sdss_dr12_spectrum(name, filenamesDict, whichFiles='icaSpectra'):
    """
    :param name:
    :param filenamesDict:
    :param whichFiles: 'icaSpectra' or 'sdssSpectra' are the two possible arguments
    :return:
    """

    if whichFiles == 'icaSpectra':
        filepath = filepathsDict[name]['ica']
    elif whichFiles == 'sdssSpectra':
        filepath = filepathsDict[name]['sdss']

    hdulist = fits.open(filepath)

    data = hdulist[1].data
    hdr = hdulist[0].header
    hdulist.close()

    wavelength = 10 ** np.array([j[1] for j in data])
    dw = wavelength * (10 ** hdr['COEFF1'] - 1.0)
    flux = np.array([j[0] for j in data])
    err = np.sqrt(abs(np.array(
        [j[6] for j in data])))  # square root of sky #abs added by mjt91 in an attempt to remove error messages

    # create the mask array
    mask = []
    for j in range(len(wavelength)):
        if 5574 < wavelength[j] < 5583:
            mask.append(0)
        elif 6297 < wavelength[j] < 6306:
            mask.append(0)
        else:
            mask.append(1)

    return wavelength, dw, flux, err, mask


if __name__ == '__main__':
    wavelength, dw, flux, err, mask = get_boss_dr12_spec('000000.66+145828.8')

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.plot(wavelength, mask * flux)
    plt.show()
