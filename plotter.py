import matplotlib.pyplot as plt
import argparse
import os
import sys
import csv

parser = argparse.ArgumentParser(description='Make plot from all csv files in given directory.')
parser.add_argument('dirname', help="directory with data to plot")
parser.add_argument('-cn', '--compare-normalized', metavar='dirname2', help="normalize to respective maximums both original data and given by this option, then make plot")
parser.add_argument('-t', '--title', metavar='name', help="plot title")
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
                ys = []
                for first,second in reader:
                    try:
                        x,y = float(first), float(second)
                        if comp_max_y < y:
                            comp_max_y = y
                        xs.append(x)
                        ys.append(y)
                    except ValueError:
                        pass
                comp_data.append((xs, ys, filename))


for filename in os.listdir(args.dirname):
    if filename.endswith(".csv") or filename.endswith(".txt"):
        with open(args.dirname+"/"+filename, newline='') as csvfile:
            reader = csv.reader(csvfile)
            xs = []
            ys = []
            for first,second in reader:
                try:
                    x,y = float(first), float(second)
                    if max_y < y:
                        max_y = y
                    xs.append(x)
                    ys.append(y)
                except ValueError:
                    pass
            data.append((xs, ys, filename))


fig, axs = plt.subplots(1)
for (xs, ys, name) in data:
    if max_norm:
        ys = [y / max_y for y in ys]
    axs.plot(xs,ys, label=name)
if args.compare_normalized is not None:
    for (xs, ys, name) in comp_data:
        axs.plot(xs,[y / comp_max_y for y in ys], label=name)

axs.legend()

fig.canvas.manager.set_window_title(' '.join(sys.argv[1:]))
plt.title(args.title)
plt.show()