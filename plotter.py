import matplotlib.pyplot as plt
import argparse
import os
import sys
import csv

parser = argparse.ArgumentParser(description='Make plot from all csv files in given directory.')
parser.add_argument('dirname', help="directory with data to plot")
parser.add_argument('-cn', '--compare-normalized', metavar='dirname2', help="normalize to respective maximums both original data and given by this option, then make plot")
parser.add_argument('-t', '--title', metavar='name', help="plot title")
parser.add_argument('-s', '--stack', help="plot title", action='store_true')
args = parser.parse_args()


data = []
comp_data = []
max_norm = False
max_y = 0
if args.compare_normalized is not None:
    max_norm = True
    comp_max_y = 0
    for filename in os.listdir(args.compare_normalized):
        if filename.endswith(".csv") or filename.endswith(".txt"):
            with open(args.compare_normalized+"/"+filename, newline='') as csvfile:
                reader = csv.reader(csvfile)
                xs = []
                yrows = []
                for row in reader:
                    try:
                        x = float(row[0])
                        ys = []
                        for y in row[1:]:
                            if comp_max_y < float(y):
                                comp_max_y = float(y)
                            ys.append(float(y))
                        xs.append(x)
                        yrows.append(ys)
                    except ValueError:
                        pass
                comp_data.append((xs, yrows, filename))


for filename in os.listdir(args.dirname):
    if filename.endswith(".csv") or filename.endswith(".txt"):
        with open(args.dirname+"/"+filename, newline='') as csvfile:
            reader = csv.reader(csvfile)
            xs = []
            yrows = []
            for row in reader:
                try:
                    x = float(row[0])
                    ys = []
                    for y in row[1:]:
                        if max_y < float(y):
                            max_y = float(y)
                        ys.append(float(y))
                    xs.append(x)
                    yrows.append(ys)
                except ValueError:
                    pass
            data.append((xs, yrows, filename))

fig, axs = plt.subplots(1)
if args.stack:
    yrows = [[row[i] for row in yrows] for i in range(len(yrows[0]))]
    try:
        axs.stackplot(xs, yrows)
    except ValueError:
        pass
else:
    for (xs, yrows, name) in data:
        if max_norm:
            yrows = [[y / max_y for y in yrow] for yrow in yrows]
        axs.plot(xs, yrows, label=name)
    if args.compare_normalized is not None:
        for (xs, yrows, name) in comp_data:
            axs.plot(xs, [[y / max_y for y in yrow] for yrow in yrows], label=name)

axs.legend()

fig.canvas.manager.set_window_title(' '.join(sys.argv[1:]))
plt.title(args.title)
plt.show()