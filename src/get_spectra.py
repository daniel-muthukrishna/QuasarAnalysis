import json
from astropy.io import fits


def load_spectra_filenames(spectraFile='../data_files/created/spectraFilepaths.json'):
    with open(spectraFile, 'r') as f:
        filenamesDict = json.load(f)

    return filenamesDict


def get_sdss_dr12_spectrum(name, filepathsDict):
    """
    :param name:
    :param filenamesDict:
    :param whichFiles: 'icaSpectra' or 'sdssSpectra' are the two possible arguments
    :return:
    """

    filepathICA = filepathsDict[name]['ica']
    filepathSDSS = filepathsDict[name]['sdss']

    with fits.open(filepathICA, memmap=False) as hdulist:
        z = hdulist[0].header['Z_ICA']
        flux = hdulist[0].data[0]

    with fits.open(filepathSDSS, memmap=False) as hdulist:
        data = hdulist[2].data
        mag = data['PSFMAG'][0]

    return flux, z, mag


if __name__ == '__main__':
    wavelength, dw, flux, err, mask = get_boss_dr12_spec('000000.66+145828.8')

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.plot(wavelength, mask * flux)
    plt.show()
