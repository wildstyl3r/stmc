import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
from natsort import natsorted
import os
import sys
import csv
import pathlib
# import scienceplots
from matplotlib import rc
from enum import Enum
rc('text', usetex=True)
rc('text.latex',preamble=r'\usepackage{amsmath}')
rc('text.latex',preamble=r'\usepackage[T1,T2A]{fontenc}')
rc('text.latex',preamble=r'\usepackage[utf8]{inputenc}')
rc('text.latex',preamble=r'\usepackage[russian]{babel}')

# plt.style.use('retro')

parser = argparse.ArgumentParser(description='Make plot from all csv files in given directory.')
parser.add_argument('query', nargs='*', help="path to the data to plot, with wildcards or modificators")
parser.add_argument('-e', '--exclude', nargs='*', metavar='filterx', help="filter out by substring in filename")
parser.add_argument('-f', '--filter', nargs='*', metavar='filterx', help="filter out by substring in filename")
parser.add_argument('-c', '--compare', nargs='*', metavar='dirname2', help="compare two data sources")
parser.add_argument('-t', '--title', metavar='name', help="plot title")
parser.add_argument('-n', '--norm', help="normalize data by maximum if comparing", action='store_true')
parser.add_argument('-s', '--switch', help="switch x and y columns for one dimensional data files", action='store_true')
parser.add_argument('-x', '--horizontal', metavar='horizontal_axis', help="use that column as x")
parser.add_argument('-y', '--vertical', metavar='vertical_axis', help="use that column as y")
parser.add_argument('-m', '--margin', metavar='vertical_margin', help="use that column as y's margin")
parser.add_argument('-lm', '--lower-margin', metavar='vertical_lower_margin', help="use that column as y's lower error margin")
parser.add_argument('-um', '--upper-margin', metavar='vertical_upper_margin', help="use that column as y's upper error margin")
parser.add_argument('-im', '--ignore-margin', help="ignore margin", action='store_true')
parser.add_argument('-w', '--eval', metavar='eval', help="add symbolic expression to plot" )
parser.add_argument('-r', '--restartcc', help="restart color cycle at",type=int)
parser.add_argument('-ly' '--log-y', help="make y logscale", action='store_true')
args = parser.parse_args()

exclude = args.exclude if args.exclude is not None else []
filters = args.filter if args.filter is not None else []

files = []
markers = ["o-", "v-", "^-"]
marker = 0

for item in pathlib.Path().glob(args.query[0]):
    skip = False
    for exc in exclude:
        if exc in str(item):
            skip = True
            break
    take = True
    if len(filters) > 0:
        take = False
        for f in filters:
            if f in str(item):
                take = True
                break
    if skip or not take:
        continue
    files.append(str(item))
    print(str(item))

files = natsorted(files, lambda f: f[0])


def lists_intersection(lists):
    if len(lists) > 0:
        intersection = set(lists[0])
        for l in lists:
            intersection = intersection.intersection(set(l))
        return intersection
    return []

axis_names = None

# common_parts_set = lists_intersection([f.split('/') for (f,_) in files])
# title = args.title or ' / '.join([seg if seg in common_parts_set else "***" for seg in files[0][0].split('/')]) if len(files) > 0 else args.title

dataframes = []

def get_occur_index(l, query):
    if query is not None:
        for (i, s) in enumerate(l):
            if query in s:
                return i
        else:
            try:
                t = int(query)
                return t
            except:
                return None

densities = dict()


for f in files:
    model, config, line = f.split('/')[-3:]
    densities[config] = [[] for _ in range(5)]

for f in files:
    model, config, line = f.split('/')[-3:]
    i = int(line[-1])
    p = int(line.split('_')[0])
    with open(f+'/PlasmaDensity.csv', newline='') as csvfile:
        spamreader = csv.reader(csvfile)
        first = True
        for row in spamreader:
            if first:
                first = False
                continue
            if (float(row[0]) - 1.5)**2 < 0.02:
                densities[config][i-1].append((p/3, float(row[1])))
                break
    densities[config][i-1] = sorted(densities[config][i-1], key=lambda x: x[0])

for config in densities:
    data = densities[config]
    for i in range(5):
        print(config+f"n_e_{i+1}mA.csv")
        path = '/'.join(files[0].split('/')[:-4])+'/middle_density/'+config
        os.makedirs(path, exist_ok=True)
        with open(path+f"/n_e_{i+1}mA.csv", 'w', newline='') as csvfile:
            spamwriter = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
            spamwriter.writerow(["X (Pa)","n_e (cm^-3)"])
            for row in densities[config][i]:
                spamwriter.writerow([str(row[0]),str(row[1])])