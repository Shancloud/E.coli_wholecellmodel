from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt
import six
from six.moves import cPickle, zip

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import units
import math

THRESHOLD = 1e-13 # roughly, the mass of an electron

FG_PER_DALTON = 1.6605402e-9

# TODO: get these from the KB
REPRESENTATIVE_MASSES = {
	"proton":1.007 * FG_PER_DALTON,
	"amino acid":109 * FG_PER_DALTON,
	"ATP":551 * FG_PER_DALTON,
	"protein":40e3 * FG_PER_DALTON,
	"ribosome":2700e3 * FG_PER_DALTON
	}

def mergefigure6(path, exportpath):
    n = len(path)
    plt.figure(figsize=(8.5*n, 11))
    list_index = []
    list_avgProcessMassDifferences = []
    list_mass = []
    listy = []
    for pathx in path:
        with open('/'.join(pathx.split('/')[:-6])+'/kb'+'/simData.cPickle', 'rb') as f:
            sim_data = cPickle.load(f)

            # Listeners used
        main_reader = TableReader(os.path.join(pathx, "Main"))
        mass = TableReader(os.path.join(pathx, "Mass"))
        fba_results = TableReader(os.path.join(pathx, "FBAResults"))

        # Sim time
        initialTime = main_reader.readAttribute("initialTime")
        time = main_reader.readColumn("time") - initialTime

        # Mass differences
        processNames = mass.readAttribute("processNames")
        processMassDifferences = mass.readColumn("processMassDifferences")

        # Adjust metabolism for exchange fluxes
        ## Exchange fluxes will not capture rounding to single molecule level
        ## so need to use change in metabolites as the mass coming into the cell
        metabolites = fba_results.readAttribute('outputMoleculeIDs')
        delta_metabolites = fba_results.readColumn('deltaMetabolites')

        conversion = 1e15 / sim_data.constants.n_avogadro.asNumber(1 / units.mol)
        metabolism_mass_imported = np.zeros(delta_metabolites.shape[0])
        for mol, flux in zip(metabolites, delta_metabolites.T):
            mol_mass = sim_data.getter.get_mass(mol).asNumber(units.g / units.mol)
            metabolism_mass_imported += mol_mass * conversion * flux

        metabolism_mass_difference = processMassDifferences[:, processNames.index('Metabolism')]
        adjusted_metabolism = (metabolism_mass_difference - metabolism_mass_imported).reshape(-1, 1)
        processMassDifferences = np.hstack((adjusted_metabolism, processMassDifferences))

        # Average differences over cell cycle
        processNames = ['Metabolism\n(exchange adjusted)'] + processNames
        avgProcessMassDifferences = np.abs(processMassDifferences).sum(axis=0) / len(time)
        index = np.arange(len(processNames))
        width = 1

        list_index.append(index)
        list_avgProcessMassDifferences.append(avgProcessMassDifferences)
        list_mass.append(mass)
        listy.extend(avgProcessMassDifferences*avgProcessMassDifferences)

    for idx in range(n):
        axes = plt.subplot(1,n, idx + 1)

        r1 = axes.barh(list_index[idx], list_avgProcessMassDifferences[idx] * (list_avgProcessMassDifferences[idx] > THRESHOLD), width, log = True, color = (0.9, 0.2, 0.2))
        r2 = axes.barh(list_index[idx], list_avgProcessMassDifferences[idx] * (list_avgProcessMassDifferences[idx] <= THRESHOLD), width, log = True, color = (0.2, 0.2, 0.9))

        axes.set_yticks(list_index[idx]+width/2)
        axes.set_yticklabels(processNames) #, rotation = -45)

        axes.plot([THRESHOLD, THRESHOLD], [list_index[idx][0], list_index[idx][-1]+width], 'k--', linewidth=3)

        plt.text(THRESHOLD, list_index[idx][-1], "electron", rotation = "vertical", va = "center", ha = "right")

        for name, list_mass[idx] in REPRESENTATIVE_MASSES.items():
            plt.axvline(list_mass[idx])
            plt.text(list_mass[idx], list_index[idx][-1], name, rotation = "vertical", va = "center", ha = "right")

        minx = np.array([i for i in list(map(lambda x1:x1.tolist(), listy)) if i >0]).min()
        maxx = np.array(list(map(lambda x1:x1.tolist(), listy))).max()

        plt.xlim(10**(math.log10(minx)),10**(math.log10(maxx)))

        plt.xlabel("Mass difference (fg)")

        plt.title("Average absolute change in mass by individual processes")

        plt.tight_layout()
        plt.grid(True, which = "major")

    plt.savefig(exportpath + 'processMassBalance.pdf')
