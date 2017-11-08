import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import pearsonr, chisquare
from sklearn.cluster import KMeans
import corner
from chainconsumer import ChainConsumer
from component_reconstruction import plot_spectrum, smooth


def load_spectra(savedSpectra='spectra.pickle'):
    with open(savedSpectra, 'rb') as f:
        spectra = pickle.load(f)

    return spectra


def compare_reconstruction_and_data(name, spectra):
    sdssFlux = smooth(spectra[name]['sdssFlux'])
    reconFlux = spectra[name]['reconFlux']
    pearsonCorr = round(pearsonr(sdssFlux, reconFlux)[0], 3)
    chi2 = 0  # round(chisquare(sdssFlux, reconFlux)[0], 0)

    return pearsonCorr, chi2


def compare_all(spectra):
    pearsonVals = []
    removeNames = []
    names = list(spectra.keys())
    for name in names:
        pearson, chi2 = compare_reconstruction_and_data(name, spectra)
        pearsonVals.append(pearson)
        bal = 'BAL' if spectra[name]['balFlag'] else 'non-BAL'
        if pearson < 0.66 or chi2 > 300:
            removeNames.append(name)
            print(pearson, bal, chi2)
            print("################## {0} ".format(name))
            # plot_spectrum(name, spectra, addToTitle='PC={0}_{1}'.format(pearson, chi2))

    meanPearson = np.mean(pearsonVals)
    stdPearson = np.std(pearsonVals)
    print("Average Pearson: {0}, {1}".format(meanPearson, stdPearson))

    return removeNames


def get_weights(spectra, removeNames, lineType='All'):
    weightsArray = []
    namesShortenedList = []
    weightLabels = ['w1', 'w2', 'w3', 'w4', 'w5', 'w6']
    names = list(spectra.keys())
    for name in names:
        if name not in removeNames:
            if lineType == 'All':
                weightsArray.append(spectra[name]['weights'])
                namesShortenedList.append(name)
            elif lineType == 'BALs':
                if spectra[name]['balFlag'] == 1:
                    weightsArray.append(spectra[name]['weights'])
                    namesShortenedList.append(name)
            elif lineType == 'nonBALs':
                if spectra[name]['balFlag'] == 0:
                    weightsArray.append(spectra[name]['weights'])
                    namesShortenedList.append(name)
            else:
                raise Exception("Invalid lineType argument specified! Arg value must be 'All', 'BALs', or 'nonBALs'")

    weightsArray = np.array(weightsArray)
    weightsDF = pd.DataFrame(data=weightsArray, columns=weightLabels)

    return weightsDF, pd.DataFrame(namesShortenedList)


def plot_weights(spectra, removeNames):
    weightsDF, namesDF = get_weights(spectra, removeNames, lineType='All')
    weightsDFBals, namesDFBals = get_weights(spectra, removeNames, lineType='BALs')
    weightsDFNonBals, namesDFNonBals = get_weights(spectra, removeNames, lineType='nonBALs')
    corner.corner(weightsDF)

    c = ChainConsumer()
    c.add_chain(weightsDFBals.values, parameters=weightsDF.columns.tolist())
    c.add_chain(weightsDFNonBals.values, parameters=weightsDF.columns.tolist())
    # c.add_chain(weightsDF.values, parameters=weightsDF.columns.tolist())

    c.plotter.plot()


def clustering(spectra, removeNames, numClusters=3):
    weightsDF, namesDF = get_weights(spectra, removeNames)
    weightNames = weightsDF.columns.tolist()

    kmeans = KMeans(n_clusters=numClusters)
    kmeans.fit(weightsDF)
    labels = kmeans.predict(weightsDF)
    centroids = kmeans.cluster_centers_
    centroidsDF = pd.DataFrame(data=centroids, columns=weightNames)

    fig, ax = plt.subplots(nrows=6, ncols=6, sharex='col', sharey='row')
    fig.subplots_adjust(wspace=0, hspace=0)
    colmap = {1: 'r', 2: 'g', 3: 'b', 4: 'c'}
    colors = list(map(lambda x: colmap[x + 1], labels))

    for i, label1 in enumerate(weightNames):
        for j, label2 in enumerate(weightNames):
            if label1 != label2 and i >= j:
                ax[i, j].scatter(weightsDF[label1], weightsDF[label2], color=colors, alpha=0.3)
                for idx, centroid in centroidsDF.iterrows():
                    pass  #ax[i, j].scatter(centroid[label1], centroid[label2], color=colmap[idx+1], s=20)
                if i == 5:
                    ax[i, j].set_xlabel(label2)
                if j == 0:
                    ax[i, j].set_ylabel(label1)

    clusters = {}
    for i in range(numClusters):
        clusters[i] = {'weights': [], 'names': [], 'reconFluxes': [], 'sdssFluxes': []}
    for (idx, weightsRow), (idx2, nameDF) in zip(weightsDF.iterrows(), namesDF.iterrows()):
        name = nameDF.values[0]
        clusters[labels[idx]]['weights'].append(weightsRow.values)
        clusters[labels[idx]]['names'].append(name)
        clusters[labels[idx]]['reconFluxes'].append(spectra[name]['reconFlux'])
        clusters[labels[idx]]['sdssFluxes'].append(spectra[name]['sdssFlux'])
    for i in range(numClusters):
        clusters[labels[i]]['weights'] = np.array(clusters[labels[i]]['weights'])
        clusters[labels[i]]['names'] = np.array(clusters[labels[i]]['names'])
        clusters[labels[i]]['reconFluxes'] = np.array(clusters[labels[i]]['reconFluxes'])
        clusters[labels[i]]['sdssFluxes'] = np.array(clusters[labels[i]]['sdssFluxes'])

    c = ChainConsumer()
    for clusterName in clusters.keys():
        c.add_chain(clusters[clusterName]['weights'], parameters=weightNames)
    c.plotter.plot()

    return clusters


def analyse_clusters(clusters, wave):
    clusterNames = list(clusters.keys())
    numClusters = len(clusterNames)

    plt.figure("Clusters_k{0}_avg_spectra_withSDSS".format(numClusters))
    plt.title("Clusters_ k{0}".format(numClusters))
    for clusterName in clusterNames:
        np.savetxt("Clusters_k{0}_c{1}.txt".format(numClusters, clusterName), clusters[clusterName]['names'], fmt="%s")
        reconFluxes = clusters[clusterName]['reconFluxes']
        reconFluxAvg = np.average(reconFluxes, axis=0)
        plt.plot(wave, reconFluxAvg, label="Cluster_{0}_recon".format(clusterName))
        sdssFluxes = clusters[clusterName]['sdssFluxes']
        sdssFluxAvg = np.average(sdssFluxes, axis=0)
        plt.plot(wave, sdssFluxAvg, label="Cluster_{0}_sdss".format(clusterName))
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Flux')
    plt.legend()

if __name__ == '__main__':
    spectra1 = load_spectra()
    names1 = list(spectra1.keys())
    wave1 = spectra1[names1[0]]['reconWave']
    removeNames1 = compare_all(spectra1)
    plot_weights(spectra1, removeNames1)
    clusters1 = clustering(spectra1, removeNames1)
    analyse_clusters(clusters1, wave1)
    # for idx in range(10):
    #     plot_spectrum(names1[idx], spectra1)
    plt.show()
