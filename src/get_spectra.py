import json
from astropy.io import fits


def load_spectra_filenames(spectraFile='spectraFilepaths.json'):
    with open(spectraFile, 'r') as f:
        filenamesDict = json.load(f)

    return filenamesDict


def get_sdss_dr12_spectrum(name, filepathsDict, whichFiles='icaSpectra'):
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
    else:
        print('Invalid whichFiles argument!')
        return

    with fits.open(filepath, memmap=False) as hdulist:
        z = hdulist[0].header['Z_ICA']
        flux = hdulist[0].data[0]

    return flux, z


if __name__ == '__main__':
    wavelength, dw, flux, err, mask = get_boss_dr12_spec('000000.66+145828.8')

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.plot(wavelength, mask * flux)
    plt.show()
