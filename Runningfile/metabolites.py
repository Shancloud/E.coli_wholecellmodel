"""
Shows fold change of metabolites over the course of the simulation
"""

from __future__ import absolute_import, division, print_function

from six.moves import cPickle
import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.plotting_tools import COLORS_LARGE
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


def mergefigure18(path, exportpath):
    n = len(path)
    plt.figure(figsize=(8.5, 11*n))
    list_time = []
    list_normalizedCounts = []
    list_aa_mask = []
    list_sorted_idx = []
    list_metaboliteNames = []
    listx = []
    listy = []
    for pathx in path:
        with open('/'.join(pathx.split('/')[:-6])+'/kb'+'/simData.cPickle', 'rb') as f:
            sim_data = cPickle.load(f)
        aa_ids = sim_data.molecule_groups.amino_acids

        # Listeners used
        enzymeKineticsdata = TableReader(os.path.join(pathx, "EnzymeKinetics"))
        main_reader = TableReader(os.path.join(pathx, "Main"))

        # Metabolite data
        metaboliteNames = np.array(enzymeKineticsdata.readAttribute("metaboliteNames"))
        metaboliteCounts = enzymeKineticsdata.readColumn("metaboliteCountsFinal")
        normalizedCounts = metaboliteCounts / metaboliteCounts[1, :]
        aa_mask = np.array([m in aa_ids for m in metaboliteNames])

        # Highlight outliers
        mean_final = normalizedCounts[-1, ~aa_mask].mean()
        std_final = normalizedCounts[-1, ~aa_mask].std()
        highlighted = (
            (normalizedCounts[-1, :] > mean_final + 3*std_final)
            | (normalizedCounts[-1, :] < mean_final - 3*std_final)
        )

        # Sort amino acids for labeling
        sorted_idx = np.argsort(normalizedCounts[-1, aa_mask])[::-1]

        # Read time info from the listener
        initialTime = main_reader.readAttribute("initialTime")
        time = (main_reader.readColumn("time") - initialTime) / 60

        list_time.append(time)
        list_normalizedCounts.append(normalizedCounts)
        list_aa_mask.append(aa_mask)
        list_sorted_idx.append(sorted_idx)
        list_metaboliteNames.append(sorted_idx)
        listx.extend(time)
        listy.extend(normalizedCounts[:, aa_mask][:, sorted_idx])

        # Plot everything but amino acids
        ax = plt.subplot(n+1,1,1)
        #ax.set_prop_cycle('color')

        ## Plot and label metabolites that are different from the mean
        mask = ~aa_mask & highlighted
        if np.any(mask):
            plt.plot(time, normalizedCounts[:, mask])

        ## Plot the rest of the metabolites that are not amino acids
        mask = ~aa_mask & ~highlighted
        if np.any(mask):
            plt.plot(time, normalizedCounts[:, mask])

        ## Formatting
        plt.xlabel("Time (min)")
        plt.ylabel("Metabolite fold change")
        plt.title('All metabolites (excluding amino acids)')

    for idx in range(n):
        # Plot only amino acids
        ax = plt.subplot(n+1, 1, idx + 2)
        colors = COLORS_LARGE
        ax.set_prop_cycle('color', colors)
        plt.plot(list_time[idx], list_normalizedCounts[idx][:, list_aa_mask[idx]][:, list_sorted_idx[idx]])
        plt.legend(metaboliteNames[aa_mask][sorted_idx], ncol=2)
        plt.xlim(min(listx)-3, max(listx)+3)
        plt.ylim(np.array(list(map(lambda x1:x1.tolist(), listy))).min()-0.5, np.array(list(map(lambda x1:x1.tolist(), listy))).max()+0.5)
        plt.xlabel("Time (min)")
        plt.ylabel("Metabolite fold change - amino acids only")
        plt.title('Only amino acids')

        plt.tight_layout()

    plt.savefig(exportpath + 'metabolites.pdf')


