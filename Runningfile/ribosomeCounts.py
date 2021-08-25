from __future__ import absolute_import, division, print_function

import os

from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from matplotlib.pyplot import MultipleLocator

def mergefigure11(path, exportpath):

    plt.figure(figsize = (8.5, 11))
    for pathx in path:

        uniqueMoleculeCounts = TableReader(pathx + "UniqueMoleculeCounts")
        ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index('active_ribosome')

        main_reader = TableReader(pathx + "Main")
        initialTime = main_reader.readAttribute("initialTime")
        time = main_reader.readColumn("time") - initialTime
        nActive = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]

        plt.plot(time / 60, nActive)
        plt.plot([time[0] / 60., time[-1] / 60.], [2 * nActive[0], 2 * nActive[0]], "r--")

    plt.xlabel("Time (min)")
    plt.ylabel("Counts")
    plt.title("Active Ribosomes Final:Initial")
    plt.savefig(exportpath + 'ribosomeCounts.pdf')
