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
parser.add_argument('-lx' '--log-x', help="make x logscale", action='store_true')
args = parser.parse_args()

exclude = args.exclude if args.exclude is not None else []
filters = args.filter if args.filter is not None else []

files = []
markers = ["o-", "v-", "^-"]
marker = 0

for q in args.query:
    fields = q.split(">>", 2)
    d,scale = "",1
    if len(fields) == 2:
        scale = float(fields[0])
        d = fields[1]
    else:
        d = fields[0]
    print(f"scale is {scale}")
    for item in pathlib.Path().glob(d):
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
        files.append((str(item), scale))
        print(str(item))

files = natsorted(files, lambda f: f[0])

class TableType(Enum):
    ONE_DIM_DISTRIBUTION = 1
    ONE_DIM_DISTRIBUTION_WITH_MARGINS = 2
    AGGREGATION = 3
    ONE_DIM_DISTRIBUTION_WITH_TWO_MARGINS = 4

def get_offsets_and_table_type(filepath):
    i = 0
    data_offset = None
    with open(filepath, newline='') as csvfile:
        r = csv.reader(csvfile)
        prev_row = None
        for row in r:
            prev_row = row
            if len(row) > 1:
                if len(row) == 2:
                    return i, TableType.ONE_DIM_DISTRIBUTION
                if len(row) == 3 and row[2] == "margin":
                    return i, TableType.ONE_DIM_DISTRIBUTION_WITH_MARGINS
                if  len(row) == 4 and row[2] == "marginL"and row[3] == "marginU":
                    return i, TableType.ONE_DIM_DISTRIBUTION_WITH_TWO_MARGINS
                else:
                    return i, TableType.AGGREGATION
            i += 1
        return 0, []

def lists_intersection(lists):
    if len(lists) > 0:
        intersection = set(lists[0])
        for l in lists:
            intersection = intersection.intersection(set(l))
        return intersection
    return []

axis_names = None

common_parts_set = lists_intersection([f.split('/') for (f,_) in files])
title = args.title or ' / '.join([seg if seg in common_parts_set else "***" for seg in files[0][0].split('/')]) if len(files) > 0 else args.title

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



def add_math_to_str(s: str) -> str:
    measure_start_index = s.rfind('(')
    measure_end_index = s.rfind(')')
    if measure_end_index != -1 and measure_start_index != -1:
        s = s[:measure_end_index] + " $ " + s[measure_end_index:]
        s = s[:measure_start_index+1] + " $ " + s[measure_start_index+1:]
    print(s)
    return s

fig, axs = plt.subplots(1)
for (i,(filename, scale)) in enumerate(files):
    if args.restartcc is not None and i > 0 and i % args.restartcc == 0:
        axs.set_prop_cycle(None)
        marker = (marker + 1)  % len(markers)
    label = '/'.join([seg for seg in filename.split('/') if seg not in common_parts_set]).removesuffix('.txt').removesuffix('.csv')
    if scale != 1:
        label = "{:.2e} x ".format(scale) + label
    
    offset, table_type = get_offsets_and_table_type(filename)
    dataframe = pd.read_csv(filename, skiprows=offset)

    if table_type == TableType.AGGREGATION:
        x = get_occur_index(dataframe.columns, args.horizontal) or 0
        y = get_occur_index(dataframe.columns, args.vertical) or 1
        m = get_occur_index(dataframe.columns, args.margin) or None
        margin_name = dataframe.columns[m] if m is not None else None
        if label == "":
            label = add_math_to_str(dataframe.columns[y])
        if axis_names is None:
            axis_names = [dataframe.columns[x], dataframe.columns[y]]
        elif axis_names != [dataframe.columns[x], dataframe.columns[y]]:
            print(f"AGGREGATION unexpected column header: wanted {axis_names}, got {[dataframe.columns[x], dataframe.columns[y]]}")
        axs.errorbar(dataframe[dataframe.columns[x]], dataframe[dataframe.columns[y]], yerr=dataframe[margin_name] if margin_name is not None else (dataframe[dataframe.columns[y]+"Margin"] if dataframe.columns[y]+"Margin" in dataframe.columns else None),label=label,fmt=markers[marker], ms=3)
    elif table_type==TableType.ONE_DIM_DISTRIBUTION_WITH_TWO_MARGINS:
        x = get_occur_index(dataframe.columns, args.horizontal) or 0
        y = get_occur_index(dataframe.columns, args.vertical) or 1
        lm = get_occur_index(dataframe.columns, "marginL")
        um = get_occur_index(dataframe.columns, "marginU")
        if x == y:
            y = 1-x
        if label == "":
            label = filename.split('/')[-1].removesuffix('.txt').removesuffix('.csv')
        if axis_names is None:
            axis_names = [dataframe.columns[x], dataframe.columns[y]]
        elif axis_names != [dataframe.columns[x], dataframe.columns[y]]:
            print(f"SCALAR DISTRIBUTION unexpected column header: wanted {axis_names}, got {[dataframe.columns[x], dataframe.columns[y]]}, there are {dataframe.columns} columns")
        # (x, y) = (dataframe[dataframe.columns[1]], dataframe[dataframe.columns[0]]) if args.switch else (dataframe[dataframe.columns[0]], dataframe[dataframe.columns[1]])
        # axs.errorbar(x, y, yerr=dataframe["margin"] if "margin" in dataframe.columns else None,label=label,fmt="o-", ms=3)
        margin = []
        
        if lm is not None and um is not None:
            margin = np.array([dataframe[dataframe.columns[lm]], dataframe[dataframe.columns[um]]])

        axs.errorbar(dataframe[dataframe.columns[x]], dataframe[dataframe.columns[y]], yerr=margin,label=label,fmt=markers[marker], ms=3)
    else:
        x = get_occur_index(dataframe.columns, args.horizontal) or 0
        y = get_occur_index(dataframe.columns, args.vertical) or 1
        m = get_occur_index(dataframe.columns, "margin")
        if x == y:
            y = 1-x
        if label == "":
            label = filename.split('/')[-1].removesuffix('.txt').removesuffix('.csv')
        if axis_names is None:
            axis_names = [dataframe.columns[x], dataframe.columns[y]]
        elif axis_names != [dataframe.columns[x], dataframe.columns[y]]:
            print(f"SCALAR DISTRIBUTION unexpected column header: wanted {axis_names}, got {[dataframe.columns[x], dataframe.columns[y]]}, there are {dataframe.columns} columns")
        # (x, y) = (dataframe[dataframe.columns[1]], dataframe[dataframe.columns[0]]) if args.switch else (dataframe[dataframe.columns[0]], dataframe[dataframe.columns[1]])
        # axs.errorbar(x, y, yerr=dataframe["margin"] if "margin" in dataframe.columns else None,label=label,fmt="o-", ms=3)
        margin = []
        if m is not None:
            margin = dataframe[dataframe.columns[m]]
        else:
            if dataframe.columns[y]+"margin" in dataframe.columns:
                margin = dataframe[dataframe.columns[y]+"Margin"]

        axs.errorbar(dataframe[dataframe.columns[x]], scale*dataframe[dataframe.columns[y]], yerr=margin if args.ignore_margin is False and m is not None else None, fmt=markers[marker], label=label, ms=3)

if args.eval is not None:
    low,high,lmd = args.eval.split('#')
    x = [i for i in range(int(low),int(high), (int(high)-int(low))//100)]
    f = eval(lmd)
    y = [f(i) for i in x]
    
    axs.errorbar(x, y, label = lmd, fmt='--',ms=3)

# if args.compare:
#     datafiles_c, _, max_y_c, lineLabels_c = read_datafiles(args.compare)
#     plt.gca().set_prop_cycle(None)
#     if args.norm:
#         max_y_c = [datafiles_c[filename].max() for filename in datafiles_c]
#         print('dir-wide max is [',max_y_c, ']')
#         for filename in datafiles_c:
#             for column in datafiles_c[filename]:
#                 label = f"{filename}: {datafiles_c[filename][column].name}"
#                 if args.filter is None or all(f in label for f in include):
#                     if args.exclude is None or all(e not in label for e in exclude):
#                         axs.errorbar(datafiles_c[filename].index, column,
#                                     fmt="-", capsize=3,
#                                     yerr=None,
#                                     label=label,
#                                     data=datafiles_c[filename])
#     else:
#         for filename in datafiles_c:
#             for column in datafiles_c[filename]:
#                 if not column.endswith("_conf_interval"):
# if args.compare:
#     datafiles_c, _, max_y_c, lineLabels_c = read_datafiles(args.compare)
#     plt.gca().set_prop_cycle(None)
#     if args.norm:
#         max_y_c = [datafiles_c[filename].max() for filename in datafiles_c]
#         print('dir-wide max is [',max_y_c, ']')
#         for filename in datafiles_c:
#             for column in datafiles_c[filename]:
#                 label = f"{filename}: {datafiles_c[filename][column].name}"
#                 if args.filter is None or all(f in label for f in include):
#                     if args.exclude is None or all(e not in label for e in exclude):
#                         axs.errorbar(datafiles_c[filename].index, column,
#                                     fmt="-", capsize=3,
#                                     yerr=None,
#                                     label=label,
#                                     data=datafiles_c[filename])
#     else:
#         for filename in datafiles_c:
#             for column in datafiles_c[filename]:
#                 if not column.endswith("_conf_interval"):
#                     label = f"{filename}: {datafiles_c[filename][column].name}"
#                     if args.filter is None or all(f in label for f in include):
#                         
#                     label = f"{filename}: {datafiles_c[filename][column].name}"
#                     if args.filter is None or all(f in label for f in include):
#                         if args.exclude is None or all(e not in label for e in exclude):
#                             axs.errorbar(
#                                 datafiles_c[filename].index, column,
#                                 fmt="-", capsize=3,
#                                 yerr=column + "_conf_interval" if column + "_conf_interval" in datafiles_c[filename] else None,
#                                 label=label,data=datafiles_c[filename])

# axs.legend()
# if len(measures) > 0:
#     axs.set_xlabel(f"{measures[0]}")
# if len(measures) > 1:
#     axs.set_ylabel(f"{measures[1]}")

axs.legend()
if axis_names is not None:

    axs.set_xlabel(f"{add_math_to_str(axis_names[0])}")
    axs.set_ylabel(f"{add_math_to_str(axis_names[1])}")

fcm = fig.canvas.manager
if fcm is not None:
    fcm.set_window_title(' '.join(sys.argv[1:]))
if args.norm:
    title +="peak ratio: " + str(max_y_c / max_y)
plt.title(title)
print(args)
if args.ly__log_y is not None and args.ly__log_y:
    axs.set_yscale('log')
if args.lx__log_x is not None and args.lx__log_x:
    axs.set_xscale('log')
plt.show()