import json
from astropy.io import fits


def load_spectra_filenames(spectraFile='../data_files/created/spectraFilepaths.json'):
    with open(spectraFile, 'r') as f:
        filenamesDict = json.load(f)

    return filenamesDict


def get_sdss_dr12_spectrum(name, filepathsDict):
    filepathICA = filepathsDict[name]['ica']
    filepathSDSS = filepathsDict[name]['sdss']

    with fits.open(filepathICA, memmap=False) as hdulist:
        z = hdulist[0].header['Z_ICA']
        flux = hdulist[0].data[0]

    otherInfo = {}
    with fits.open(filepathSDSS, memmap=False) as hdulist:
        data = hdulist[2].data
        otherInfo['mag'] = data['PSFMAG'][0]
        otherInfo['magErr'] = data['PSFMAGERR'][0]

    return flux, z, otherInfo


if __name__ == '__main__':
    filenamesDict = load_spectra_filenames()
    flux, z, otherInfo = get_sdss_dr12_spectrum('000000.66+145828.8')

