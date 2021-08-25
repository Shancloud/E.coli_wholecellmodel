from aaCounts import mergefigure1
from dntpCounts import mergefigure2
from compartment_mass_fraction_summary import mergefigure15
from massFractionSummary import mergefigure3
from mrnaCounts import mergefigure4
from ntpCounts import mergefigure5
from processMassBalance import mergefigure6
from processMassBalanceDynamics import mergefigure7
from proteinCounts import mergefigure8
from replication import mergefigure9
from ribosomeCapacity import mergefigure10
from ribosomeCounts import mergefigure11
from rnapCapacity import mergefigure12
from rnapCounts import mergefigure13
from rnaseCounts import mergefigure14
from evalutionTime import mergefigure16
from external_exchange_fluxes import mergefigure17
from metabolites import mergefigure18
import sys

path = sys.argv[1:-1]
exportpath = sys.argv[-1]
print(path)
print('exportpath = ' + exportpath)

mergefigure1(path, exportpath)
mergefigure2(path, exportpath)
mergefigure3(path, exportpath)
mergefigure4(path, exportpath)
mergefigure5(path, exportpath)
mergefigure6(path, exportpath)
mergefigure7(path, exportpath)
mergefigure8(path, exportpath)
mergefigure9(path, exportpath)
mergefigure10(path, exportpath)
mergefigure11(path, exportpath)
mergefigure12(path, exportpath)
mergefigure13(path, exportpath)
mergefigure14(path, exportpath)
mergefigure15(path, exportpath)
mergefigure16(path, exportpath)
mergefigure17(path,exportpath)
mergefigure18(path,exportpath)
