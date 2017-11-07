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
    weightLabels = ['w1', 'w2', 'w3', 'w4', 'w5', 'w6']
    names = list(spectra.keys())
    for name in names:
        if name not in removeNames:
            if lineType == 'All':
                weightsArray.append(spectra[name]['weights'])
            elif lineType == 'BALs':
                if spectra[name]['balFlag'] == 1:
                    weightsArray.append(spectra[name]['weights'])
            elif lineType == 'nonBALs':
                if spectra[name]['balFlag'] == 0:
                    weightsArray.append(spectra[name]['weights'])
            else:
                raise Exception("Invalid lineType argument specified! Arg value must be 'All', 'BALs', or 'nonBALs'")

    weightsArray = np.array(weightsArray)
    weightsDF = pd.DataFrame(data=weightsArray, columns=weightLabels)

    return weightsDF


def plot_weights(spectra, removeNames):
    weightsDF = get_weights(spectra, removeNames, lineType='All')
    weightsDFBals = get_weights(spectra, removeNames, lineType='BALs')
    weightsDFNonBals = get_weights(spectra, removeNames, lineType='nonBALs')
    corner.corner(weightsDF)

    c = ChainConsumer()
    c.add_chain(weightsDFBals.values, parameters=weightsDF.columns.tolist())
    c.add_chain(weightsDFNonBals.values, parameters=weightsDF.columns.tolist())
    # c.add_chain(weightsDF.values, parameters=weightsDF.columns.tolist())

    c.plotter.plot()

    # for label1 in weightsDF.columns.tolist():
    #     for label2 in weightsDF.columns.tolist():
    #         if label1 != label2:
    #             plt.figure("%s vs %s" % (label1, label2))
    #             plt.scatter(weightsDF[label1], weightsDF[label2])


def clustering(spectra, removeNames, numClusters=3):
    weightsDF = get_weights(spectra, removeNames)
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
        clusters[labels[i]] = []
    for idx, row in weightsDF.iterrows():
        clusters[labels[idx]].append(row.values)
    for i in range(numClusters):
        clusters[labels[i]] = np.array(clusters[labels[i]])

    c = ChainConsumer()
    for clusterName in clusters.keys():
        c.add_chain(clusters[clusterName], parameters=weightNames)
    c.plotter.plot()




if __name__ == '__main__':
    spectra1 = load_spectra()
    removeNames1 = compare_all(spectra1)
    plot_weights(spectra1, removeNames1)
    clustering(spectra1, removeNames1)
    # names1 = list(spectra1.keys())
    # for idx in range(10):
    #     plot_spectrum(names1[idx], spectra1)
    plt.show()