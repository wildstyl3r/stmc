import matplotlib.pyplot as plt
import pandas as pd
import argparse
import os
import sys
import csv

parser = argparse.ArgumentParser(description='Make plot from all csv files in given directory.')
parser.add_argument('dirname', help="directory with data to plot")
parser.add_argument('-c', '--compare', metavar='dirname2', help="compare two data sources")
parser.add_argument('-x', '--exclude', metavar='filterx', help="filter out by substring in filename")
parser.add_argument('-f', '--filter', metavar='filter', help="filter by substring in filename")
parser.add_argument('-t', '--title', metavar='name', help="plot title")
parser.add_argument('-s', '--stack', help="stack the cumulative datapoints", action='store_true')
parser.add_argument('-p', '--points', help="seek for single-point data", action='store_true')
parser.add_argument('-n', '--norm', help="normalize data by maximum if comparing", action='store_true')
parser.add_argument('-m', '--measures', metavar='measures_arg', help="compare two data sources")
args = parser.parse_args()

exclude = args.exclude.split() if args.exclude is not None else []
include = args.filter.split() if args.filter is not None else []

def read_datafile(dirname):
    datafiles, measures = dict(), []
    max_y = 0
    first = True
    for filename in sorted(os.listdir(dirname)):
        if args.exclude is None or all(e not in filename for e in exclude):
            # if args.filter is None or args.filter in filename:
            if filename.endswith(".csv") or filename.endswith(".txt"):
                datafiles[filename] = pd.read_csv(dirname+'/'+filename, header=(1 if args.measures is None else 0), index_col=0)
                if first:
                    first = False
                    if args.measures is None:
                        with open(dirname+'/'+filename, newline='') as csvfile:
                            reader = csv.reader(csvfile)
                            for row in reader:
                                measures = row
                                break
                    else:
                        measures = args.measures.split()
    return datafiles, measures, max_y

datafiles, measures, max_y = read_datafile(args.dirname)
# print(datafiles)
if args.compare is not None:
    datafiles_c, _, max_y_c = read_datafile(args.compare)

plt.rcParams['text.usetex'] = True
fig, axs = plt.subplots(1)
if args.stack:
    for filename in datafiles.keys():
        for column in datafiles[filename]:
            print(filename, datafiles[filename][column])
        if len (labels) > i:
            name += ": " + labels[i]
        try:
            axs.stackplot(xs, yrows, label=list(name + l for l in labels) if len(labels) > 0 else name)
        except ValueError:
            pass
else:
    for filename in datafiles.keys():
        for column in datafiles[filename]:
            label = f"{filename}: {datafiles[filename][column].name}"
            if args.filter is None or all(f in label for f in include):
                if args.exclude is None or all(e not in label for e in exclude):
                    axs.plot(datafiles[filename][column], "o-", label=label)
    if args.compare:
        for filename in datafiles_c:
            if args.norm:
                for (xs, yrows, name) in data_c:
                    axs.plot(xs, [[y / max_y_c for y in yrow] for yrow in yrows], "o-", label=list(name + l for l in labels) if len(labels) > 0 else name)
            else:
                for column in datafiles_c[filename]:
                    label = f"{filename}: {datafiles_c[filename][column].name}"
                    if args.filter is None or all(f in label for f in include):
                        if args.exclude is None or all(e not in label for e in exclude):
                            axs.plot(datafiles_c[filename][column], "o-", label=label)

axs.legend()
if len(measures) > 0:
    axs.set_xlabel(f"${measures[0]}$")
if len(measures) > 1:
    axs.set_ylabel(f"${measures[1]}$")

fig.canvas.manager.set_window_title(' '.join(sys.argv[1:]))
if args.norm:
    plt.title((args.title or "") + "peak ratio: " + str(max_y_c / max_y))
else:
    plt.title(args.title)
plt.show()