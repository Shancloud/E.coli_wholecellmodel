"""
Plot amino acid counts
"""

from __future__ import absolute_import, division, print_function

import os

from matplotlib import pyplot as plt
from six.moves import cPickle, range

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot



def mergefigure1(path, exportpath):
    plt.figure(figsize=(8.5, 11))
    n = len(path)

    for pathx in path:
        sim_data = cPickle.load(open('/'.join(pathx.split('/')[:-6])+'/kb'+'/simData.cPickle', 'rb'))

        aaIDs = sim_data.molecule_groups.amino_acids
        (aaCounts,) = read_bulk_molecule_counts(pathx, (aaIDs,))

        main_reader = TableReader(pathx + "Main")
        initialTime = main_reader.readAttribute("initialTime")
        time = main_reader.readColumn("time") - initialTime

        for idx in range(21):

            plt.subplot(6, 4, idx + 1)

            plt.plot(time / 60., aaCounts[:, idx], linewidth = 2)
            plt.xlabel("Time (min)")
            plt.ylabel("Counts")
            plt.title(aaIDs[idx], fontsize=8)
            plt.tick_params(labelsize=8)

        plt.tight_layout()

    plt.savefig(exportpath + 'aaCounts.pdf')
