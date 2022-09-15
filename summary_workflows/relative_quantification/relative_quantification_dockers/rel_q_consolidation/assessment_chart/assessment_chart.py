import pandas
import io
import json
import os
import numpy as np
import matplotlib
matplotlib.use("SVG")

import matplotlib.pyplot as plt
plt.ioff()

"""
    INFO:
    This module contains functions that generate an assessment 2D chart from the aggregation datasets generated in any of the OpenEBench benchmarking workflows.
    For more details about the OpenEbench scientific results visualization, please visit: https://github.com/inab/OpenEBench_scientific_visualizer

    @Javier Garrayo Ventas. Barcelona Supercomputing Center. Spain. 2019"
"""
def pareto_frontier(Xs, Ys, maxX=True, maxY=True):
    # Sort the list in either ascending or descending order of X
    myList = sorted([[Xs[i], Ys[i]] for i, val in enumerate(Xs, 0)], reverse=maxX)
    # Start the Pareto frontier with the first value in the sorted list
    p_front = [myList[0]]
    # Loop through the sorted list
    for pair in myList[1:]:
        if maxY:
            if pair[1] >= p_front[-1][1]:  # Look for higher values of Y
                p_front.append(pair)  # and add them to the Pareto frontier
        else:
            if pair[1] <= p_front[-1][1]:  # look for lower values of Y
                p_front.append(pair)  # and add them to the pareto frontier
    # Turn resulting pairs back into a list of Xs and Ys
    p_frontX = [pair[0] for pair in p_front]
    p_frontY = [pair[1] for pair in p_front]
    return p_frontX, p_frontY

# funtion that gets quartiles for x and y values
def plot_square_quartiles(x_values, means, tools, better, ax, percentile=50):
    x_percentile, y_percentile = (np.nanpercentile(x_values, percentile), np.nanpercentile(means, percentile))
    plt.axvline(x=x_percentile, linestyle='--', color='#0A58A2', linewidth=1.5)
    plt.axhline(y=y_percentile, linestyle='--', color='#0A58A2', linewidth=1.5)

    # create a dictionary with tools and their corresponding quartile
    tools_quartiles = {}
    if better == "bottom-right":

        # add quartile numbers to plot
        plt.text(0.99, 0.15, '1', verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes, fontsize=25, alpha=0.2)
        plt.text(0.01, 0.15, '2', verticalalignment='bottom', horizontalalignment='left', transform=ax.transAxes, fontsize=25, alpha=0.2)
        plt.text(0.99, 0.85, '3', verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, fontsize=25, alpha=0.2)
        plt.text(0.01, 0.85, '4', verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, fontsize=25, alpha=0.2)

        for i, val in enumerate(tools, 0):
            if x_values[i] >= x_percentile and means[i] <= y_percentile:
                tools_quartiles[tools[i]] = 1
            elif x_values[i] >= x_percentile and means[i] > y_percentile:
                tools_quartiles[tools[i]] = 3
            elif x_values[i] < x_percentile and means[i] > y_percentile:
                tools_quartiles[tools[i]] = 4
            elif x_values[i] < x_percentile and means[i] <= y_percentile:
                tools_quartiles[tools[i]] = 2

    elif better == "top-right":
        
        # add quartile numbers to plot
        plt.text(0.99, 0.85, '1', verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, fontsize=25, alpha=0.2)
        plt.text(0.01, 0.85, '2', verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, fontsize=25, alpha=0.2)
        plt.text(0.99, 0.15, '3', verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes, fontsize=25, alpha=0.2)
        plt.text(0.01, 0.15, '4', verticalalignment='bottom', horizontalalignment='left', transform=ax.transAxes, fontsize=25, alpha=0.2)

        for i, val in enumerate(tools, 0):
            if x_values[i] >= x_percentile and means[i] < y_percentile:
                tools_quartiles[tools[i]] = 3
            elif x_values[i] >= x_percentile and means[i] >= y_percentile:
                tools_quartiles[tools[i]] = 1
            elif x_values[i] < x_percentile and means[i] >= y_percentile:
                tools_quartiles[tools[i]] = 2
            elif x_values[i] < x_percentile and means[i] < y_percentile:
                tools_quartiles[tools[i]] = 4

    elif better == "top-left":
        # add quartile numbers to plot
        plt.text(0.99, 0.85, '2', verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, fontsize=25, alpha=0.2)
        plt.text(0.01, 0.85, '1', verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, fontsize=25, alpha=0.2)
        plt.text(0.99, 0.15, '4', verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes, fontsize=25, alpha=0.2)
        plt.text(0.01, 0.15, '3', verticalalignment='bottom', horizontalalignment='left', transform=ax.transAxes, fontsize=25, alpha=0.2)

        for i, val in enumerate(tools, 0):
            if x_values[i] >= x_percentile and means[i] < y_percentile:
                tools_quartiles[tools[i]] = 4
            elif x_values[i] >= x_percentile and means[i] >= y_percentile:
                tools_quartiles[tools[i]] = 2
            elif x_values[i] < x_percentile and means[i] >= y_percentile:
                tools_quartiles[tools[i]] = 1
            elif x_values[i] < x_percentile and means[i] < y_percentile:
                tools_quartiles[tools[i]] = 3

    elif better == "bottom-left":
        # add quartile numbers to plot
        plt.text(0.99, 0.85, '4', verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, fontsize=25, alpha=0.2)
        plt.text(0.01, 0.85, '3', verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, fontsize=25, alpha=0.2)
        plt.text(0.99, 0.15, '2', verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes, fontsize=25, alpha=0.2)
        plt.text(0.01, 0.15, '1', verticalalignment='bottom', horizontalalignment='left', transform=ax.transAxes, fontsize=25, alpha=0.2)

        for i, val in enumerate(tools, 0):
            if x_values[i] >= x_percentile and means[i] < y_percentile:
                tools_quartiles[tools[i]] = 2
            elif x_values[i] >= x_percentile and means[i] >= y_percentile:
                tools_quartiles[tools[i]] = 4
            elif x_values[i] < x_percentile and means[i] >= y_percentile:
                tools_quartiles[tools[i]] = 3
            elif x_values[i] < x_percentile and means[i] < y_percentile:
                tools_quartiles[tools[i]] = 1

    return (tools_quartiles)


# function to normalize the x and y axis to 0-1 range
def normalize_data(x_values, means):
    maxX = max(x_values)
    minX = min(x_values)
    maxY = max(means)
    minY = min(means)
    # maxX = ax.get_xlim()[1]
    # minX = ax.get_xlim()[0]
    # maxY = ax.get_ylim()[1]
    # minY = ax.get_ylim()[0]
    # x_norm = [(x - minX) / (maxX - minX) for x in x_values]
    # means_norm = [(y - minY) / (maxY - minY) for y in means]
    if maxX == 0.0:
        x_norm = x_values
    else:
        x_norm = [x / maxX for x in x_values]
    if maxY == 0.0:
        means_norm = means
    else:
        means_norm = [y / maxY for y in means]
    return x_norm, means_norm


# funtion that plots a diagonal line separating the values by the given quartile
def draw_diagonal_line(scores_and_values, quartile, better, max_x, max_y):
    for i, val in enumerate(scores_and_values, 0):
        # find out which are the two points that contain the percentile value
        if scores_and_values[i][0] <= quartile:
            target = [(scores_and_values[i - 1][1], scores_and_values[i - 1][2]),
                      (scores_and_values[i][1], scores_and_values[i][2])]
            break
    # get the the mid point between the two, where the quartile line will pass
    half_point = (target[0][0] + target[1][0]) / 2, (target[0][1] + target[1][1]) / 2
    # plt.plot(half_point[0], half_point[1], '*')
    # draw the line depending on which is the optimal corner
    if better == "bottom-right":
        x_coords = (half_point[0] - max_x, half_point[0] + max_x)
        y_coords = (half_point[1] - max_y, half_point[1] + max_y)
    elif better == "top-right":
        x_coords = (half_point[0] + max_x, half_point[0] - max_x)
        y_coords = (half_point[1] - max_y, half_point[1] + max_y)
    elif better == "top-left":
        x_coords = (half_point[0] + max_x, half_point[0] - max_x)
        y_coords = (half_point[1] + max_y, half_point[1] - max_y)

    elif better == "bottom-left":
        x_coords = (half_point[0] - max_x, half_point[0] + max_x)
        y_coords = (half_point[1] + max_y, half_point[1] - max_y)


    plt.plot(x_coords, y_coords, linestyle='--', color='#0A58A2', linewidth=1.5)


# funtion that splits the analysed tools into four quartiles, according to the asigned score
def get_quartile_points(scores_and_values, first_quartile, second_quartile, third_quartile):
    tools_quartiles = {}
    for i, val in enumerate(scores_and_values, 0):
        if scores_and_values[i][0] > third_quartile:
            tools_quartiles[scores_and_values[i][3]] = 1
        elif second_quartile < scores_and_values[i][0] <= third_quartile:
            tools_quartiles[scores_and_values[i][3]] = 2
        elif first_quartile < scores_and_values[i][0] <= second_quartile:
            tools_quartiles[scores_and_values[i][3]] = 3
        elif scores_and_values[i][0] <= first_quartile:
            tools_quartiles[scores_and_values[i][3]] = 4
    return (tools_quartiles)


# funtion that separate the points through diagonal quartiles based on the distance to the 'best corner'
def plot_diagonal_quartiles(x_values, means, tools, better):
    # get distance to lowest score corner

    # normalize data to 0-1 range
    x_norm, means_norm = normalize_data(x_values, means)
    max_x = max(x_values)
    max_y = max(means)
    # compute the scores for each of the tool. based on their distance to the x and y axis
    scores = []
    for i, val in enumerate(x_norm, 0):
        if better == "bottom-right":
            scores.append(x_norm[i] + (1 - means_norm[i]))
        elif better == "top-right":
            scores.append(x_norm[i] + means_norm[i])
        elif better == "top-left":
            scores.append(1-x_norm[i] + means_norm[i])
        elif better == "bottom-left":
            scores.append(1-x_norm[i] + (1 - means_norm[i]))

    # add plot annotation boxes with info about scores and tool names
    for counter, scr in enumerate(scores):
        plt.annotate(
            tools[counter] + "\n" +
            # str(round(x_norm[counter], 6)) + " * " + str(round(1 - means_norm[counter], 6)) + " = " + str(round(scr, 8)),
            "score = " + str(round(scr, 3)),
            xy=(x_values[counter], means[counter]), xytext=(0, 20),
            textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.15),
            size=7,
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

    # region sort the list in descending order
    scores_and_values = sorted([[scores[i], x_values[i], means[i], tools[i]] for i, val in enumerate(scores, 0)],
                               reverse=True)
    scores = sorted(scores, reverse=True)
    # print (scores_and_values)
    # print (scores)
    # endregion
    first_quartile, second_quartile, third_quartile = (
        np.nanpercentile(scores, 25), np.nanpercentile(scores, 50), np.nanpercentile(scores, 75))
    # print (first_quartile, second_quartile, third_quartile)
    draw_diagonal_line(scores_and_values, first_quartile, better, max_x, max_y)
    draw_diagonal_line(scores_and_values, second_quartile, better, max_x, max_y)
    draw_diagonal_line(scores_and_values, third_quartile, better, max_x, max_y)

    # split in quartiles
    tools_quartiles = get_quartile_points(scores_and_values, first_quartile, second_quartile, third_quartile)
    return (tools_quartiles)


# function that prints a table with the list of tools and the corresponding quartiles
def print_quartiles_table(tools_quartiles):
    row_names = tools_quartiles.keys()
    quartiles_1 = tools_quartiles.values()

    colnames = ["TOOL", "Quartile"] 
    celltxt = zip(row_names, quartiles_1) 
    df = pandas.DataFrame(celltxt)
    vals = df.values

    # set cell colors depending on the quartile
    # green color scale
    colors = df.applymap(lambda x: '#238b45' if x == 1 else '#74c476' if x == 2 else '#bae4b3' if x == 3
    else '#edf8e9' if x == 4 else '#ffffff')

    colors = colors.values

    the_table = plt.table(cellText=vals,
                          colLabels=colnames,
                          cellLoc='center',
                          loc='right',
                          bbox=[1.1, 0.15, 0.5, 0.8],
                          colWidths=[1.2, 0.5],
                          cellColours=colors,
                          colColours=['#ffffff'] * 2)
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12)
    plt.subplots_adjust(right=0.65, bottom=0.2)



def print_chart(challenge_dir, summary_dir, challenge_type, classification_type):

    tools = []
    x_values = []
    y_values = []
    # with io.open(summary_dir, mode='r', encoding="utf-8") as f:
    #     aggregation_file = json.load(f)

    aggregation_file = summary_dir

    # participants and their metrics
    for participant_data in aggregation_file["datalink"]["inline_data"]["challenge_participants"]:

        tools.append(participant_data['participant_id'])
        x_values.append(participant_data['metric_x'])
        y_values.append(participant_data['metric_y'])
    # metrics names for axes
    x_metric = aggregation_file["datalink"]["inline_data"]["visualization"]["x_axis"]
    y_metric = aggregation_file["datalink"]["inline_data"]["visualization"]["y_axis"]
    
    ax = plt.subplot()
    for i, val in enumerate(tools, 0):
        markers = [".", "o", "v", "^", "<", ">", "1", "2", "3", "4", "8", "s", "p", "P", "*", "h", "H", "+",
                   "x", "X",
                   "D",
                   "d", "|", "_", ","]
        colors = ['#5b2a49', '#a91310', '#9693b0', '#e7afd7', '#fb7f6a', '#0566e5', '#00bdc8', '#cf4119', '#8b123f',
                  '#b35ccc', '#dbf6a6', '#c0b596', '#516e85', '#1343c3', '#7b88be']

        ax.errorbar(x_values[i], y_values[i], linestyle='None', marker=markers[i],
                    markersize='15', markerfacecolor=colors[i], markeredgecolor=colors[i], capsize=6,
                    ecolor=colors[i], label=tools[i])

    # change plot style
    # set plot title

    plt.title(challenge_type, fontsize=18, fontweight='bold')

    ax.set_xlabel(x_metric, fontsize=12)
    ax.set_ylabel(y_metric, fontsize=12)

    # Shrink current axis's height  on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.25,
                     box.width, box.height * 0.75])

    # Put a legend below current axis
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), markerscale=0.7,
               fancybox=True, shadow=True, ncol=5, prop={'size': 12})


    # set the axis limits
    x_lims = ax.get_xlim()
    plt.xlim(x_lims)
    y_lims = ax.get_ylim()
    plt.ylim(y_lims)
    if x_lims[0] >= 1000:
        ax.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
    if y_lims[0] >= 1000:
        ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda y, loc: "{:,}".format(int(y))))

    # set parameters for optimization
    better = aggregation_file["datalink"]["inline_data"]["visualization"]["optimization"]
    max_x = True
    max_y = True

    # get pareto frontier and plot
    p_frontX, p_frontY = pareto_frontier(x_values, y_values, maxX=max_x, maxY=max_y)
    plt.plot(p_frontX, p_frontY, linestyle='--', color='grey', linewidth=1)
    # append edges to pareto frontier
    if better == 'bottom-right':
        left_edge = [[x_lims[0], p_frontX[-1]], [p_frontY[-1], p_frontY[-1]]]
        right_edge = [[p_frontX[0], p_frontX[0]], [p_frontY[0], y_lims[1]]]
        plt.plot(left_edge[0], left_edge[1], right_edge[0], right_edge[1], linestyle='--', color='red',
                 linewidth=1)

    elif better == 'top-right':
        left_edge = [[x_lims[0], p_frontX[-1]], [p_frontY[-1], p_frontY[-1]]]
        right_edge = [[p_frontX[0], p_frontX[0]], [p_frontY[0], y_lims[0]]]
        plt.plot(left_edge[0], left_edge[1], right_edge[0], right_edge[1], linestyle='--', color='red',
                 linewidth=1)

    elif better == 'top-left':
        left_edge = [[x_lims[1], p_frontX[-1]], [p_frontY[-1], p_frontY[-1]]]
        right_edge = [[p_frontX[0], p_frontX[0]], [p_frontY[0], y_lims[0]]]
        plt.plot(left_edge[0], left_edge[1], right_edge[0], right_edge[1], linestyle='--', color='red',
                 linewidth=1)

    elif better == 'bottom-left':
        left_edge = [[x_lims[1], p_frontX[-1]], [p_frontY[-1], p_frontY[-1]]]
        right_edge = [[p_frontX[0], p_frontX[0]], [p_frontY[0], y_lims[1]]]
        plt.plot(left_edge[0], left_edge[1], right_edge[0], right_edge[1], linestyle='--', color='red',
                 linewidth=1)

    # add 'better' annotation and quartile numbers to plot
    if better == 'bottom-right':
        plt.annotate('better', xy=(0.98, 0.04), xycoords='axes fraction',
                     xytext=(-30, 30), textcoords='offset points',
                     ha="right", va="bottom",
                     arrowprops=dict(facecolor='black', shrink=0.05, width=0.9))

    elif better == 'top-right':
        plt.annotate('better', xy=(0.98, 0.95), xycoords='axes fraction',
                     xytext=(-30, -30), textcoords='offset points',
                     ha="right", va="top",
                     arrowprops=dict(facecolor='black', shrink=0.05, width=0.9))

    elif better == 'top-left':
        plt.annotate('better', xy=(0.04, 0.95), xycoords='axes fraction',
                     xytext=(30, -30), textcoords='offset points',
                     ha="left", va="top",
                     arrowprops=dict(facecolor='black', shrink=0.05, width=0.9))

    elif better == 'bottom-left':
        plt.annotate('better', xy=(0.04, 0.04), xycoords='axes fraction',
                     xytext=(30, 30), textcoords='offset points',
                     ha="left", va="bottom",
                     arrowprops=dict(facecolor='black', shrink=0.05, width=0.9))

    # add chart grid
    plt.grid(b=None, which='major', axis='both', linewidth=0.5)


    if classification_type == "SQR":
        tools_quartiles = plot_square_quartiles(x_values, y_values, tools, better, ax)
        print_quartiles_table(tools_quartiles)

    elif classification_type == "DIAG":
        tools_quartiles = plot_diagonal_quartiles(x_values, y_values, tools, better)
        print_quartiles_table(tools_quartiles)

    out_id = aggregation_file["_id"].split("_")
    del out_id[0]
    out_id = "_".join(out_id)

    outname = os.path.join(challenge_dir, out_id + "_benchmark_" + classification_type + ".svg")
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    fig.savefig(outname, dpi=100)

    plt.close("all")

def print_barplot(challenge_dir, aggregation, challenge_acronym):
    """
    Print bar plots when there is only a single metric in the aggregation
    """

    tools = []
    values = []
    orange_hex = "#f47c21" #Bar plots in orange like on the OEB website

    # participants and their metrics
    for participant_data in aggregation["datalink"]["inline_data"]["challenge_participants"]:
        tools.append(participant_data['participant_id'])
        values.append(participant_data['metric_value'])

    # metrics names for axes
    metric_name = aggregation["datalink"]["inline_data"]["visualization"][
        "metric"]

    ax = plt.subplot()
    ax.bar(tools, values, color = orange_hex, edgecolor = orange_hex)

    ax.set_xlabel("Tools", fontsize=12)
    ax.set_ylabel(f"{metric_name.capitalize()}", fontsize=12)
    ax.set_title(f"{metric_name.capitalize()} bar-plot in challenge {challenge_acronym}")
    ax.legend()

    # Extract the challenge name and "Aggregation" from the id for the file name
    out_id = aggregation["_id"].split("_")
    del out_id[0]
    out_id = "_".join(out_id)

    outname = os.path.join(challenge_dir,
                           out_id + "_benchmark_" + metric_name + "_barplot.svg")
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    fig.savefig(outname, dpi=100)

    plt.close("all")
