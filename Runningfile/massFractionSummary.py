from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot
from six.moves import zip


COLORS_256 = [ # From colorbrewer2.org, qualitative 8-class set 1
    [228,26,28],
    [55,126,184],
    [77,175,74],
    [152,78,163],
    [255,127,0],
    [255,255,51],
    [166,86,40],
    [247,129,191]
    ]

COLORS = [
    [colorValue/255. for colorValue in color]
    for color in COLORS_256
    ]


def mergefigure3(path, exportpath):

    n = len(path)

    listx = []
    listy1 = []
    listy2 = []
    listfraction = []
    alllistx = []
    alllisty = []

    plt.figure(figsize=(8.5, 11*n))

    for pathx in path:
        mass = TableReader(os.path.join(pathx, "Mass"))
        main_reader = TableReader(os.path.join(pathx, "Main"))

        cell = mass.readColumn("dryMass")
        protein = mass.readColumn("proteinMass")
        tRna = mass.readColumn("tRnaMass")
        rRna = mass.readColumn("rRnaMass")
        mRna = mass.readColumn("mRnaMass")
        dna = mass.readColumn("dnaMass")
        smallMolecules = mass.readColumn("smallMoleculeMass")

        initialTime = main_reader.readAttribute("initialTime")
        t = (main_reader.readColumn("time") - initialTime) / 60.

        masses = np.vstack([
            protein,
            rRna,
            tRna,
            mRna,
            dna,
            smallMolecules,
            ]).T
        fractions = (masses / cell[:, None]).mean(axis=0)

        plt.gca().set_prop_cycle('color', COLORS)

        y1 = list(masses / masses[0, :])
        y2 = list(cell / cell[0])

        listx.append(list(t))
        listy1.append(y1)
        listy2.append(y2)
        listfraction.append(list(fractions))
        alllistx.extend(list(t))
        alllisty.extend(y1[-1])
        alllisty.extend(y2)

    for idx in range(n):
        plt.subplot(n, 1, idx + 1)
        mass_labels = ["Protein", "rRNA", "tRNA", "mRNA", "DNA", "Small Mol.s"]
        legend = [
                     '{} ({:.3f})'.format(label, fraction)
                     for label, fraction in zip(mass_labels, fractions)
                 ] + ['Total dry mass']

        plt.plot(listx[idx], listy1[idx], linewidth=2)
        plt.plot(listx[idx], listy2[idx], color='k', linestyle=':')

        minx = min(alllistx)
        maxx = max(alllistx)
        #miny = min(alllisty)
        maxy = max(alllisty)
        plt.xlim(minx - 2, maxx + 2)
        plt.ylim(1, maxy+0.5)

        plt.title("Biomass components (average fraction of total dry mass in parentheses)")
        plt.xlabel("Time (min)")
        plt.ylabel("Mass (normalized by t = 0 min)")
        plt.legend(legend, loc="best")

        plt.tight_layout()

    plt.savefig(exportpath + 'massFractionSummary.pdf')
