"""
Plot NTP counts
"""

from __future__ import absolute_import, division, print_function

import os

from matplotlib import pyplot as plt
from six.moves import cPickle, range

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot


def mergefigure2(path, exportpath):

    plt.figure(figsize=(8.5, 11))

    for pathx in path:
        sim_data = cPickle.load(open('/'.join(pathx.split('/')[:-6])+'/kb'+'/simData.cPickle', 'rb'))

        dntpIDs = sim_data.molecule_groups.dntps
        (dntpCounts,) = read_bulk_molecule_counts(pathx, (dntpIDs,))

        main_reader = TableReader(pathx + "Main")
        initialTime = main_reader.readAttribute("initialTime")
        time = main_reader.readColumn("time") - initialTime

        for idx in range(4):

            plt.subplot(2, 2, idx + 1)

            plt.plot(time / 60., dntpCounts[:, idx], linewidth = 2)
            plt.xlabel("Time (min)")
            plt.ylabel("Counts")
            plt.title(dntpIDs[idx])

        plt.tight_layout()
    plt.savefig(exportpath + 'dntpCounts.pdf')
