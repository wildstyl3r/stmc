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
args = parser.parse_args()

def read_datafile(dirname):
    datafiles, measures = dict(), []
    max_y = 0
    first = True
    for filename in sorted(os.listdir(dirname)):
        if args.exclude is None or args.exclude not in filename:
            # if args.filter is None or args.filter in filename:
            if filename.endswith(".csv") or filename.endswith(".txt"):
                datafiles[filename] = pd.read_csv(dirname+'/'+filename, header=1, index_col=0)
                # print(datafiles[filename])
                # with open(dirname+'/'+filename, newline='') as csvfile:
                #     reader = csv.reader(csvfile)
                #     xs, yrows = [], []
                #     for row in reader:
                #         if len(row) > 1:
                #             try:
                #                 x = float(row[0])
                #                 ys = []
                #                 for y in row[1:]:
                #                     try:
                #                         if max_y < float(y):
                #                             max_y = float(y)
                #                         ys.append(float(y))
                #                     except:
                #                         pass
                #                 xs.append(x)
                #                 yrows.append(ys)
                #             except ValueError:
                #                 if len(measures) == 0 and first:
                #                     measures = row
                #                 elif len(labels) == 0 and first:
                #                     labels = row
                #                 pass
                #     data.append((xs, yrows, filename))

            # first = False
    return datafiles, measures, max_y

datafiles, measures, max_y = read_datafile(args.dirname)
# print(datafiles)
if args.compare is not None:
    datafiles_c, _, max_y_c = read_datafile(args.compare)

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
            print(column, datafiles[filename][column])
        yrows = [[row[i] for row in yrows] for i in range(len(yrows[0]))]
        if len (labels) > i:
            name += ": " + labels[i]
        try:
            axs.stackplot(xs, yrows, label=list(name + l for l in labels) if len(labels) > 0 else name)
        except ValueError:
            pass
    for (i, (xs, yrows, name)) in enumerate(datafiles):
        if args.norm:
            yrows = [[y / max_y for y in yrow] for yrow in yrows]
        else:
            yrows = [[y / 1 for y in yrow] for yrow in yrows]
        axs.plot(xs, yrows, label=list(name + l for l in labels) if len(labels) > 0 else name)#, marker='.', linewidths=0)
    if args.compare:
        if args.norm:
            for (xs, yrows, name) in data_c:
                axs.plot(xs, [[y / max_y_c for y in yrow] for yrow in yrows], label=list(name + l for l in labels) if len(labels) > 0 else name)
        else:
            for (xs, yrows, name) in data_c:
                axs.plot(xs, [[y / 1 for y in yrow] for yrow in yrows], label=list(name + l for l in labels) if len(labels) > 0 else name)

axs.legend()
if len(measures) > 0:
    axs.set_xlabel(measures[0])
if len(measures) > 1:
    axs.set_ylabel(measures[1])

fig.canvas.manager.set_window_title(' '.join(sys.argv[1:]))
if args.norm:
    plt.title((args.title or "") + "peak ratio: " + str(max_y_c / max_y))
else:
    plt.title(args.title)
plt.show()