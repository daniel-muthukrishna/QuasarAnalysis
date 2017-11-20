import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import KMeans
import corner
from chainconsumer import ChainConsumer


def get_weights(spectra, removeNames, lineType='All'):
    """ Get the Weights DataFrame for each of the spectra """
    weightsArray = []
    namesShortenedList = []
    names = list(spectra.keys())
    numComps = len(spectra[names[0]]['weights'])
    weightLabels = ['w{0}'.format(i+1) for i in range(numComps)]
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


def plot_weights(spectra, removeNames, saveDir):
    """ Plot the weights parameter space """
    weightsDF, namesDF = get_weights(spectra, removeNames, lineType='All')
    weightsDFBals, namesDFBals = get_weights(spectra, removeNames, lineType='BALs')
    weightsDFNonBals, namesDFNonBals = get_weights(spectra, removeNames, lineType='nonBALs')
    # corner.corner(weightsDF)

    c = ChainConsumer()
    c.add_chain(weightsDFBals.values, parameters=weightsDF.columns.tolist(), name='BALs')
    c.add_chain(weightsDFNonBals.values, parameters=weightsDF.columns.tolist(), name='non-BALs')
    # c.add_chain(weightsDF.values, parameters=weightsDF.columns.tolist(), name='All')
    c.configure(smooth=None)

    c.plotter.plot(filename="{0}/weights_bal_vs_nonBal.png".format(saveDir))

    # for label1 in weightsDF.columns.tolist():
    #     for label2 in weightsDF.columns.tolist():
    #         if label1[1] > label1[1]:
    #             plt.figure("%s vs %s" % (label1, label2))
    #             plt.scatter(weightsDF[label1], weightsDF[label2], alpha=0.3)


def clustering(spectra, removeNames, numClusters=3, plotClusters=True, saveDir='Figures'):
    """ Attempt to cluster the weights """
    weightsDF, namesDF = get_weights(spectra, removeNames)

    kmeans = KMeans(n_clusters=numClusters)
    kmeans.fit(weightsDF)
    labels = kmeans.predict(weightsDF)
    centroids = kmeans.cluster_centers_

    clusters = {}
    for i in range(numClusters):
        clusters[i] = {'weights': [], 'names': [], 'reconFluxes': [], 'sdssFluxes': [], 'balCounts': 0, 'nonBalCounts': 0}
    for (idx, weightsRow), (idx2, nameDF) in zip(weightsDF.iterrows(), namesDF.iterrows()):
        name = nameDF.values[0]
        clusters[labels[idx]]['weights'].append(weightsRow.values)
        clusters[labels[idx]]['names'].append(name)
        clusters[labels[idx]]['reconFluxes'].append(spectra[name]['reconFlux'])
        clusters[labels[idx]]['sdssFluxes'].append(spectra[name]['sdssFlux'])
        if spectra[name]['balFlag']:
            clusters[labels[idx]]['balCounts'] += 1
        else:
            clusters[labels[idx]]['nonBalCounts'] += 1

    for i in range(numClusters):
        clusters[labels[i]]['weights'] = np.array(clusters[labels[i]]['weights'])
        clusters[labels[i]]['names'] = np.array(clusters[labels[i]]['names'])
        clusters[labels[i]]['reconFluxes'] = np.array(clusters[labels[i]]['reconFluxes'])
        clusters[labels[i]]['sdssFluxes'] = np.array(clusters[labels[i]]['sdssFluxes'])
        print("BAL COUNTS", i, clusters[labels[i]]['balCounts'], clusters[labels[i]]['nonBalCounts'])

    if plotClusters is True:
        # saveDir = os.path.join(saveDir, 'cluster_distributions')
        # if not os.path.exists(saveDir):
        #     os.makedirs(saveDir)
        plot_clusters(clusters, weightsDF, centroids, labels, saveDir)

    return clusters


def plot_clusters(clusters, weightsDF, centroids, labels, saveDir):
    weightNames = weightsDF.columns.tolist()
    centroidsDF = pd.DataFrame(data=centroids, columns=weightNames)
    numClusters = len(clusters.keys())

    # fig, ax = plt.subplots(nrows=6, ncols=6, sharex='col', sharey='row')
    # fig.subplots_adjust(wspace=0, hspace=0)
    # colmap = {1: 'r', 2: 'g', 3: 'b', 4: 'c', 5: 'm', 6: 'y'}
    # colors = list(map(lambda x: colmap[x + 1], labels))
    #
    # for i, label1 in enumerate(weightNames):
    #     for j, label2 in enumerate(weightNames):
    #         if label1 != label2 and i >= j:
    #             ax[i, j].scatter(weightsDF[label1], weightsDF[label2], color=colors, alpha=0.3)
    #             for idx, centroid in centroidsDF.iterrows():
    #                 pass  #ax[i, j].scatter(centroid[label1], centroid[label2], color=colmap[idx+1], s=20)
    #             if i == 5:
    #                 ax[i, j].set_xlabel(label2)
    #             if j == 0:
    #                 ax[i, j].set_ylabel(label1)
    # fig.savefig("{0}/clusters_distribution_k{1}.png".format(saveDir, numClusters))

    c = ChainConsumer()
    for clusterName in clusters.keys():
        c.add_chain(np.array(clusters[clusterName]['weights']), parameters=weightNames, name='c{0}'.format(clusterName))
    c.configure(smooth=None)
    c.plotter.plot(filename="{0}/clusters_contours_k{1}.png".format(saveDir, numClusters))



def analyse_clusters(clusters, wave, saveDir):
    """ Plot the reconstructions of the median spectrum from each cluster """
    clusterNames = list(clusters.keys())
    numClusters = len(clusterNames)

    plt.figure()
    plt.title("Cluster\_k{0}".format(numClusters))
    for clusterName in clusterNames:
        np.savetxt("Clusters_k{0}_c{1}.txt".format(numClusters, clusterName), clusters[clusterName]['names'], fmt="%s")
        reconFluxes = clusters[clusterName]['reconFluxes']
        reconFluxAvg = np.median(reconFluxes, axis=0)
        plt.plot(wave, reconFluxAvg, label="$Cluster{0}\_recon$".format(clusterName))
        sdssFluxes = clusters[clusterName]['sdssFluxes']
        sdssFluxAvg = np.median(sdssFluxes, axis=0)
        plt.plot(wave, sdssFluxAvg, label="$Cluster{0}\_sdss$".format(clusterName))
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Flux')
    plt.legend()
    plt.savefig("{0}/Clusters_k{1}_avg_spectra_withSDSS.png".format(saveDir, numClusters))
