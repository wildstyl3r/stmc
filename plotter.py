import matplotlib.pyplot as plt
import argparse
import os
import sys
import csv

parser = argparse.ArgumentParser(description='Make plot from all csv files in given directory.')
parser.add_argument('dirname', help="directory with data to plot")
parser.add_argument('-c', '--compare', metavar='dirname2', help="compare two data sources")
parser.add_argument('-t', '--title', metavar='name', help="plot title")
parser.add_argument('-s', '--stack', help="plot title", action='store_true')
parser.add_argument('-n', '--norm', help="normalize data by maximum if comparing", action='store_true')
args = parser.parse_args()

def read_datafile(dirname):
    data, measures, labels = [], [], []
    max_y = 0
    for filename in sorted(os.listdir(dirname)):
        if filename.endswith(".csv") or filename.endswith(".txt"):
            with open(dirname+'/'+filename, newline='') as csvfile:
                reader = csv.reader(csvfile)
                xs, yrows = [], []
                for row in reader:
                    if len(row) > 1:
                        try:
                            x = float(row[0])
                            ys = []
                            for y in row[1:]:
                                try:
                                    if max_y < float(y):
                                        max_y = float(y)
                                    ys.append(float(y))
                                except:
                                    if len(measures) == 0:
                                        measures = row
                                    elif len(labels) == 0:
                                        labels = row
                            xs.append(x)
                            yrows.append(ys)
                        except ValueError:
                            pass
                data.append((xs, yrows, filename))
    return data, measures, labels, max_y

data, measures, labels, max_y = read_datafile(args.dirname)
if args.compare is not None:
    data_c, _, _, max_y_c = read_datafile(args.compare)

fig, axs = plt.subplots(1)
if args.stack:
    for (i, (xs, yrows, name)) in enumerate(data):
        yrows = [[row[i] for row in yrows] for i in range(len(yrows[0]))]
        if len (labels) > i:
            name += labels[i]
        try:
            axs.stackplot(xs, yrows, label=name)
        except ValueError:
            pass
else:
    for (i, (xs, yrows, name)) in enumerate(data):
        if args.norm:
            yrows = [[y / max_y for y in yrow] for yrow in yrows]
        else:
            yrows = [[y / 1 for y in yrow] for yrow in yrows]
        if len (labels) > i:
            name += labels[i]
        axs.plot(xs, yrows, label=name)
    if args.compare:
        if args.norm:
            for (xs, yrows, name) in data_c:
                axs.plot(xs, [[y / max_y_c for y in yrow] for yrow in yrows], label=name)
        else:
            for (xs, yrows, name) in data_c:
                axs.plot(xs, [[y / 1 for y in yrow] for yrow in yrows], label=name)

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