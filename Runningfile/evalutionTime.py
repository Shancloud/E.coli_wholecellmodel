from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot
from six.moves import zip
import sys

path = sys.argv[1:-1]
list_update_queries = []
list_partition = []
list_merge = []
list_calculate_mass = []
list_calculate_request = []
list_evolve_state = []
list_update = []
list_append = []
list_state_names = []
list_process_names = []
list_listener_names = []
list_logger_names = []
list_time = []
listx = []
listy1 = []
listy2 = []
listy3 = []
listy4 = []
listy5 = []
listy6 = []
listy7 = []
listy8 = []

for pathx in path:
    evaluationTime = TableReader(os.path.join(pathx, "EvaluationTime"))
    mainReader = TableReader(os.path.join(pathx, "Main"))

    state_names = evaluationTime.readAttribute("state_names")
    process_names = evaluationTime.readAttribute("process_names")
    listener_names = evaluationTime.readAttribute("listener_names")
    logger_names = evaluationTime.readAttribute("logger_names")

    clock_times = evaluationTime.readColumn("clock_time")
    update_queries = evaluationTime.readColumn("update_queries_times")
    partition = evaluationTime.readColumn("partition_times")
    merge = evaluationTime.readColumn("merge_times")
    calculate_mass = evaluationTime.readColumn("calculate_mass_times")
    calculate_request = evaluationTime.readColumn("calculate_request_times")
    evolve_state = evaluationTime.readColumn("evolve_state_times")
    update = evaluationTime.readColumn("update_times")
    append = evaluationTime.readColumn("append_times")

    initialTime = mainReader.readAttribute("initialTime")
    time = (mainReader.readColumn("time") - initialTime) / 60  # min

    list_append.append(append)
    list_merge.append(merge)
    list_update.append(update)
    list_partition.append(partition)
    list_logger_names.append(logger_names)
    list_calculate_mass.append(calculate_mass)
    list_calculate_request.append(calculate_request)
    list_evolve_state.append(evolve_state)
    list_process_names.append(process_names)
    list_listener_names.append(listener_names)
    list_state_names.append(state_names)
    list_update_queries.append(update_queries)
    list_time.append(time)
    listx.extend(time)
    listy1.extend(update_queries)
    listy2.extend(partition)
    listy3.extend(merge)
    listy4.extend(calculate_mass)
    listy5.extend(calculate_request)
    listy6.extend(evolve_state)
    listy7.extend(update)
    listy8.extend(append)

def subplot(gs, x, y,yrange, title, labels, sort=False):
    assert y.ndim == 2, 'y.ndim={}, title={}, labels={}'.format(y.ndim, title, labels)
    ax = plt.subplot(gs)
    if sort:
        idx = np.argsort(y[-1, :])[::-1]
    else:
        idx = np.arange(y.shape[1])

    # Determine legends from total time
    total_time = (y / 60).sum(axis=0)
    legend_labels = np.array(['{} ({:.2f})'.format(name, t)
        for name, t in zip(labels, total_time)])

    # Plot
    ax.semilogy(x, 1000 * y[:, idx])

    # Formatting
    ax.grid(True, which='major')
    ax.set_xlabel('Simulation time (min)')
    ax.set_ylabel('Evaluation time (ms)')
    ax.set_title(title + ' (total {:.2f} mins)'.format(total_time.sum()))
    ax.set_xlim(min(listx), max(listx)+3)
    ax.set_ylim((np.array(list(map(lambda x1:x1.tolist(), yrange))).min())*1000, (np.array(list(map(lambda x1:x1.tolist(), yrange))).max()+0.01)*1000)
    ax.legend(legend_labels[idx], bbox_to_anchor=(1,1), prop={'size':6},
        loc='upper left')

def mergefigure16(path2, exportpath):
    n = len(path2)
    fig = plt.figure(figsize=(12, 15*n))
    gs = GridSpec((4*n), 2)
    for a in range(n):
        subplot(gs=gs[(0+4*a), 0], x=list_time[a], y=list_update_queries[a] ,yrange= listy1, title ='State.updateQueries', labels=list_state_names[a])
        subplot(gs=gs[(1+4*a), 0], x=list_time[a], y=list_partition[a],yrange= listy2,title='State.partition', labels=list_state_names[a])
        subplot(gs=gs[(2+4*a), 0], x=list_time[a], y=list_merge[a], yrange= listy3,title='State.merge', labels=list_state_names[a])
        subplot(gs=gs[(3+4*a), 0], x=list_time[a], y=list_calculate_mass[a], yrange= listy4,title='State.calculateMass', labels=list_state_names[a])
        subplot(gs=gs[(0+4*a), 1], x=list_time[a], y=list_calculate_request[a], yrange= listy5,title='Process.calculateRequest', labels=list_process_names[a], sort=True)
        subplot(gs=gs[(1+4*a), 1], x=list_time[a], y=list_evolve_state[a], yrange= listy6,title='Process.evolveState', labels=list_process_names[a], sort=True)
        subplot(gs=gs[(2+4*a), 1], x=list_time[a], y=list_update[a], yrange= listy7,title='Listener.update', labels=list_listener_names[a], sort=True)
        subplot(gs=gs[(3+4*a), 1], x=list_time[a], y=list_append[a], yrange= listy8,title='Logger.append', labels=list_logger_names[a])

        gs.tight_layout(fig)
        gs.update(top=0.95)
    plt.savefig(exportpath + 'evalutionTime.pdf')
