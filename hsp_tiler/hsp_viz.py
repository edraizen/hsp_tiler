#!/usr/bin/python2.7
#Author: Eli Draizen 
#File: hsp_viz.py

#Standard libraries
import sys
import argparse
from itertools import izip
import collections
import math

#Import other libraries
import prettyplotlib as ppl
import numpy as np
from prettyplotlib import plt
import matplotlib as mpl
#import matplotlib.pyplot as plt
from prettyplotlib import brewer2mpl
from matplotlib.backends.backend_pdf import PdfPages

#Import custum libraries
from hsp_score import read_scores

"""View the improvements of HSP-Tiler in scatterplots
and histograms"""

def analyze(original, updated, scatterName=None, histName=None, log=False):
    """plot both scatterplot and histogram"""
    scatter(original, updated, save=scatterName)
    histogram(original, updated, save=histName, log=log)


   
def scatter(original, updated, main="", save=None):
    """Plot a scatterplot of updated bitscores vs. original bitscores

    """ 
    #Remove hits with no improvement and calcate the number of hits with no
    #improvement(udated == original), positive imporvent (updated > original), 
    #and negative improvment (updated < original)
    print len(original)
    positiveImprovement = []
    negativeImprovement = []
    noImprovement = 0
    for o, u in izip(original, updated):
        if int(o) == int(u):
            noImprovement +=1
        elif u > o:
            positiveImprovement.append((o,u))
        elif u < o:
            negativeImprovement.append((o,u))
        else:
            noImprovement +=1

    if not positiveImprovement:
        positiveImprovement = [()]
    if not negativeImprovement:
        negativeImprovement = [()]

    print positiveImprovement
    print negativeImprovement
    print noImprovement

    #Set deimensions
    x, y = zip(*positiveImprovement+negativeImprovement)
    xMax = int(round(sorted(x)[-1]/500.0)*500.0)
    yMax = int(round(sorted(y)[-1]/500.0)*500.0)
    sep = 500
    xticks = range(0, xMax, sep)
    yticks = range(0,yMax,sep)
    color_cycle = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors

    fig, ax = plt.subplots()
    ax.set_title(main)
    ax.set_xlabel("Original Bitscores")
    ax.set_ylabel("Updated Bitscores")

    
    #Plot postive improvement (green, automatically by prettyplotlib)
    if positiveImprovement:
        ppl.scatter(ax, *zip(*positiveImprovement), 
                    label="Positive Improvement ({} seqs)".format(len(positiveImprovement)),
                    color=color_cycle[0])

    #Draw no improvement line
    ppl.plot(ax, (0,xMax), (0,xMax), color='k', linestyle='-', linewidth=2,
             label="No Improvement ({} seqs)".format(noImprovement))

    #Plot negative improvement (red, automatically by prettyplotlib)
    if negativeImprovement:
        ppl.scatter(ax, *zip(*negativeImprovement),
                    label="Negative Improvement ({} seqs)".format(len(negativeImprovement)),
                    color=color_cycle[1])

    #Draw labels
    ppl.legend(ax)

    #Set axis
    ax.set_ylim([0,yMax])
    ax.set_xlim([0,xMax])

    if save is None:
        plt.show()
    else:
        pp = PdfPages(save)
        pp.savefig(fig)
        pp.close()

def histogram(original, updated, bins=None, main="", save=None, log=False):
    """Plot a histogram of score improvements (updated-origianl)

    Input:
    original - list of original scores
    updated - list of updates scores in same order as original
    bins - number of bins to represent improvements
    """
    #Lengths of score lists must be identical, assume in same order
    assert len(original) == len(original)

    #Set up bins:
    if bins is not None and bins > 0:
        imoprovements = {(-1,-1):0}
        for i in xrange(0, len(original), bins):
            improvements[(0,i+bins)] = 0
    else:
        improvements = {(-1,-1):0, (-5,0):0, (0,1):0, (1,25):0, (25,50):0, (50,75):0, (75,100):0, (100,125):0, (125,150):0, (150,200):0, (200,300):0, (300,400):0, (500,10000):0} #defaultdict(int)
    
    #Calcualte improvements
    for o, u in izip(original, updated):
        if o>u: 
            improvements[(-1,-1)] += 1
            continue
        for lower, upper in improvements:
            if lower <= int(u-o) < upper:
                improvements[(lower,upper)] += 1
                break
    keys = sorted(improvements.keys(), key=lambda x:x[0])
    values = [improvements[r] for r in keys]

    fig, ax = plt.subplots()
    ax.set_title(main)
    ax.set_xlabel("Improvement (updated-original) bitscores")
    ax.set_ylabel("log(Frequency)")
    #ax.set_yscale('log')

    width = 1.0
    #ax.set_xticks(np.arange(len(improvements)))
    #ax.set_xticklabels([l for l, u in keys])
    bar(ax, np.arange(len(improvements)), values, log=log,
        annotate=True, grid='y', xticklabels=[l for l, u in keys])

    if save is None:
        plt.show()
    else:
        plt.savefig(save)

def bar(*args, **kwargs):
    """
    Creates a bar plot, with white outlines and a fill color that defaults to
     the first teal-ish green in ColorBrewer's Set2. Optionally accepts
     grid='y' or grid='x' to draw a white grid over the bars,
     to show the scale. Almost like "erasing" some of the plot,
     but it adds more information!

    Can also add an annotation of the height of the barplots directly onto
    the bars with the `annotate` parameter, which can either be True,
    which will annotate the values, or a list of strings, which will annotate
    with the supplied strings.

    Can support stacked bars with the value of each stack shown on the stack
    (Added by Salil Banerjee)

    This fixes the bug with log scale in prettyplotlib

    @param ax: matplotlib.axes instance
    @param left: Vector of values of where to put the left side of the bar
    @param height: Vector of values of the bar heights
    @param kwargs: Besides xticklabels, which is a prettyplotlib-specific
    argument, any additional arguments to matplotlib.bar(): http://matplotlib
    .org/api/axes_api.html#matplotlib.axes.Axes.bar is accepted.
    """
    ax, args, kwargs = maybe_get_ax(*args, **kwargs)
    color_cycle = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
    almost_black = '#262626'
    kwargs.setdefault('color', color_cycle[0])
    kwargs.setdefault('edgecolor', 'white')
    middle = 0.4 if 'width' not in kwargs else kwargs['width']/2.0

    # Check if data contains stacks
    stacked = kwargs.pop('stacked',False)
    # Check if stack text should be included
    stack_text = kwargs.pop('stack_text',False)
    # Get legend if available
    legend = kwargs.pop('legend',False)

    left = args[0]
    height = np.array(args[1])

    # Label each individual bar, if xticklabels is provided
    xtickabels = kwargs.pop('xticklabels', None)
    # left+0.4 is the center of the bar
    xticks = np.array(left) + middle

    # Whether or not to annotate each bar with the height value
    annotate = kwargs.pop('annotate', False)

    show_ticks = kwargs.pop('show_ticks', False)

    # If no grid specified, don't draw one.
    grid = kwargs.pop('grid', None)

    # Check if stacked and plot data accordingly
    if stacked:
        num_stacks, num_data = height.shape
        bottom = np.zeros(num_data)
        for i in np.arange(num_stacks):
            lst = list(args)
            lst[1] = height[i]
            args = tuple(lst)
            kwargs['color'] = set2[i]
            kwargs['bottom'] = bottom
            rectangles = ax.bar(*args, **kwargs)
            bottom += height[i]
    else:
        rectangles = ax.bar(*args, **kwargs)

    # add legend
    if isinstance(legend, collections.Iterable):
        ax.legend(legend,loc='upper center',bbox_to_anchor=(0.5,1.11), ncol=5)

    # add whitespace padding on left
    xmin, xmax = ax.get_xlim()
    xmin -= 0.2
    if stacked:
        xmax = num_data
    ax.set_xlim(xmin, xmax)

    # If the user is only plotting one bar, make it an iterable
    if not isinstance(height, collections.Iterable):
        height = [height]


    # If there are negative counts, remove the bottom axes
    # and add a line at y=0
    if any(h < 0 for h in height.tolist()):
        axes_to_remove = ['top', 'right', 'bottom']
        ax.hlines(y=0, xmin=xmin, xmax=xmax,
                      linewidths=0.75)
    else:
        axes_to_remove = ['top', 'right']

    # Remove excess axes
    remove_chartjunk(ax, axes_to_remove, grid=grid, show_ticks=show_ticks)

    if stacked:
        data = height
        height = height.sum(axis=0)

    # Add the xticklabels if they are there
    if xtickabels is not None:
        ax.set_xticks(xticks)
        ax.set_xticklabels(xtickabels)

    if annotate or isinstance(annotate, collections.Iterable):
        annotate_yrange_factor = 0.025
        ymin, ymax = ax.get_ylim()
        yrange = ymax - ymin

        # Reset ymax and ymin so there's enough room to see the annotation of
        # the top-most
        if ymax > 0:
            ymax += yrange * 0.1
        if ymin < 0:
            ymin -= yrange * 0.1
        ax.set_ylim(ymin, ymax)
        yrange = ymax - ymin

        offset_ = math.log(yrange) + math.log(annotate_yrange_factor+1)
        print offset_
        print yrange * annotate_yrange_factor
        print math.log(yrange) + math.log(annotate_yrange_factor)
        if isinstance(annotate, collections.Iterable):
            annotations = map(str, annotate)
        else:
            annotations = ['%.3f' % h if type(h) is np.float_ else str(h)
                               for h in height]

        for x, h, annotation in zip(xticks, height, annotations):
            # Adjust the offset to account for negative bars
            offset = offset_ if h >= 0 else -1 * offset_
            verticalalignment = 'bottom' if h >= 0 else 'top'

            # Finally, add the text to the axes
            ax.annotate(annotation, (x, h + annotate_yrange_factor), 
                        verticalalignment=verticalalignment,
                        horizontalalignment='center',
                        color=almost_black)

    # Text for each block of stack
    # This was partially inspired by the following article by Tableau software
    # http://www.tableausoftware.com/about/blog/2014/1/new-whitepaper-survey-data-less-ugly-more-understandable-27812
    if stack_text:
        bottom = np.zeros(num_data)
        max_h = max(height)
        for i in np.arange(num_stacks):
            for x, d, b in zip(xticks, data[i], bottom):
                if (d*100.0/max_h) > 4.0:
                    ax.text(x,b+d/2.0,d, ha='center', va='center', color=almost_black)
            bottom += data[i]
    return rectangles

def maybe_get_ax(*args, **kwargs):
    """
    It used to be that the first argument of prettyplotlib had to be the 'ax'
    object, but that's not the case anymore.

    @param args:
    @type args:
    @param kwargs:
    @type kwargs:
    @return:
    @rtype:
    """

    if 'ax' in kwargs:
        ax = kwargs.pop('ax')
    elif len(args) == 0:
        fig = plt.gcf()
        ax = plt.gca()
    elif isinstance(args[0], mpl.axes.Axes):
        ax = args[0]
        args = args[1:]
    else:
        ax = plt.gca()
    return ax, args, dict(kwargs)

def remove_chartjunk(ax, spines, grid=None, ticklabels=None, show_ticks=False):
    '''
    Removes "chartjunk", such as extra lines of axes and tick marks.

    If grid="y" or "x", will add a white grid at the "y" or "x" axes,
    respectively

    If ticklabels="y" or "x", or ['x', 'y'] will remove ticklabels from that
    axis
    '''
    all_spines = ['top', 'bottom', 'right', 'left', 'polar']
    for spine in spines:
        # The try/except is for polar coordinates, which only have a 'polar'
        # spine and none of the others
        try:
            ax.spines[spine].set_visible(False)
        except KeyError:
            pass

    # For the remaining spines, make their line thinner and a slightly
    # off-black dark grey
    for spine in all_spines:
        if spine not in spines:
            # The try/except is for polar coordinates, which only have a 'polar'
            # spine and none of the others
            try:
                ax.spines[spine].set_linewidth(0.5)
            except KeyError:
                pass
                # ax.spines[spine].set_color(almost_black)
            #            ax.spines[spine].set_tick_params(color=almost_black)
            # Check that the axes are not log-scale. If they are, leave the ticks
            # because otherwise people assume a linear scale.
    x_pos = set(['top', 'bottom'])
    y_pos = set(['left', 'right'])
    xy_pos = [x_pos, y_pos]
    xy_ax_names = ['xaxis', 'yaxis']

    for ax_name, pos in zip(xy_ax_names, xy_pos):
        axis = ax.__dict__[ax_name]
        # axis.set_tick_params(color=almost_black)
        #print 'axis.get_scale()', axis.get_scale()
        if show_ticks or axis.get_scale() == 'log':
            # if this spine is not in the list of spines to remove
            for p in pos.difference(spines):
                #print 'p', p
                axis.set_tick_params(direction='out')
                axis.set_ticks_position(p)
                #                axis.set_tick_params(which='both', p)
        else:
            axis.set_ticks_position('none')

    if grid is not None:
        for g in grid:
            assert g in ('x', 'y')
            ax.grid(axis=grid, color='white', linestyle='-', linewidth=0.5)

    if ticklabels is not None:
        if type(ticklabels) is str:
            assert ticklabels in set(('x', 'y'))
            if ticklabels == 'x':
                ax.set_xticklabels([])
            if ticklabels == 'y':
                ax.set_yticklabels([])
        else:
            assert set(ticklabels) | set(('x', 'y')) > 0
            if 'x' in ticklabels:
                ax.set_xticklabels([])
            elif 'y' in ticklabels:
                ax.set_yticklabels([])

def parse_args():
    parser = argparse.ArgumentParser(description="Visualize HSP-Tiler output")
    parser.add_argument("-s", "--scores",
                        type=argparse.FileType('r'),
                        help="BLAST output file with updated scores")
    parser.add_argument("-l", "--log",
                        default=False,
                        action="store_true",
                        help="Make histogram on log scale")
    #Define output
    parser.add_argument("--scatterplot",
                        required=False,
                        default=None,
                        help="File to save scatterplot sequences")
    parser.add_argument("--histogram",
                        required=False,
                        default=None,
                        help="File to save scatterplot sequences")
    return parser.parse_args()

if __name__ == "__main__":
    #Parse arguments
    args = parse_args()

    print args

    original_scores, updated_scores, scores = read_scores(args.scores)

    analyze(original_scores, updated_scores, scatterName=args.scatterplot, histName=args.histogram)



