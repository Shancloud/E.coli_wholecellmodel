"""
Plot NTP counts
"""

from __future__ import absolute_import, division, print_function

import os

from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot
from six.moves import range


def mergefigure5(path, exportpath):

    plt.figure(figsize=(8.5, 11))

    for pathx in path:
        ntp_ids = ['ATP[c]', 'CTP[c]', 'GTP[c]', 'UTP[c]']
        (ntpCounts,) = read_bulk_molecule_counts(pathx, (ntp_ids,))

        main_reader = TableReader(pathx + "Main")
        initialTime = main_reader.readAttribute("initialTime")
        time = main_reader.readColumn("time") - initialTime

        for idx in range(4):

            plt.subplot(2, 2, idx + 1)

            plt.plot(time / 60., ntpCounts[:, idx], linewidth = 2)
            plt.xlabel("Time (min)")
            plt.ylabel("Counts")
            plt.title(ntp_ids[idx])

        plt.subplots_adjust(hspace = 0.5)

    plt.savefig(exportpath + 'ntpCounts.pdf')
