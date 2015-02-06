#!/usr/bin/env python

""" Plotting for the chemical tagging project. """

from __future__ import absolute_import, print_function, with_statement

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Standard library.
import cPickle as pickle
import logging

# Third-party.
import cubehelix
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FixedLocator, FixedFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Create logger.
logger = logging.getLogger("code")

def cmap_discretize(cmap, N):
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in xrange(N+1) ]
    # Return colormap object.
    return LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)


def histogram2d(results_table, column=None, bins=None, discretise=False, **kwargs):
    """
    Produce a 2D histogram showing the true number of clusters `N_true` on the
    x-axis, and the difference between the inferred and true cluster count
    (`N_inferred - N_true`) on the y-axis.

    :param results_table:
        The results table. This is expected to have the column 'N_true', and the
        third column is the number of inferred clusters (unless overwritten by
        the `column` parameter).

    :type results_table:
        :class:`astropy.table.Table`

    :param column: [optional]
        The name of the column that refers to the inferred number of clusters.
        By default the third column in `results_table` will be given.

    :type column:
        str

    :param bins: [optional]
        The number of bins to use in the histogram. This can be a positive
        integer to specify the total number of bins, or a two-length tuple of
        positive integers to specify the x- and y-axis bin counts. By default
        the number of bins will be `(max(N_true), 2 * max(N_true))`.

    :type bins:
        int or two-length tuple of ints
    """

    # Results table is expected to contain columns:
    # unique_hash, N_true, N_X
    # where N_x is the number of inferred clusters from some method.
    # If the column is not provided, just grab the third column

    if column is None:
        column = results_table.dtype.names[2]
    else:
        if column not in results_table.dtype.names:
            raise ValueError("column '{}' is not in the results table".format(
                column))

    # [TODO] allow for ax to be provided so that many axes can be on the same
    # color scale
    fig, ax = plt.subplots()

    if bins is None:
        _ = max(results_table["N_true"])
        bins = (_ - 1, 2*_)

    x = results_table["N_true"]
    y = results_table[column] - x

    # Don't use results where we could not infer the number of clusters.
    ok = results_table[column] > 0
    max_true = max(results_table["N_true"])
    heatmap, xedges, yedges = np.histogram2d(x[ok], y[ok], bins=bins,
        range=([1, max_true], [-max_true, max_true]))
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    # Create the colormap
    if discretise:
        max_colorbar_segments = kwargs.pop("__max_colorbar_segments", 10)
        colorbar_labels = np.sort(list(set(map(int, 
            np.linspace(heatmap.min(), heatmap.max(), max_colorbar_segments)))))
        cmap = cmap_discretize(cubehelix.cmap(reverse=True, start=0.3, rot=-0.5),
            len(colorbar_labels))

    else:
        cmap = cubehelix.cmap(reverse=True, start=0.3, rot=-0.5)
        
    imshow_default_kwargs = {
        # Since we will save to PDF, we cannot use 'none' for interpolation:
        # http://matplotlib.org/examples/images_contours_and_fields/interpolation_none_vs_nearest.html
        
        # Actually it's not clear which will work in the end. Might be affected
        # by the matplotlibrc file in code/
        # TODO
        "interpolation": "none",
        "origin": "lower",
        "aspect": "equal",
        "cmap": cmap
    }
    imshow_default_kwargs.update(kwargs)
    image = ax.imshow(heatmap.T, extent=extent, **imshow_default_kwargs)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "10%", pad="5%")

    colorbar = plt.colorbar(image, cax=cax)

    # Minor ticks and grid.
    minortick_frequency = 1
    ax.set_xticks(np.arange(max_true + 1)[::minortick_frequency], minor=True)
    ax.set_yticks((np.arange(2*max_true + 1) - max_true)[::minortick_frequency],
        minor=True)
    ax.grid(True, which="minor", linestyle="-", c="#cccccc")
    ax.tick_params(which="minor", width=0)
    
    # Major ticks and ticklabels.
    ax.set_xticks(ax.get_xticks() + 0.5)
    ax.set_xticklabels(["{:.0f}".format(round(_) - 1) for _ in ax.get_xticks()])
    ax.set_yticks(ax.get_yticks() + 0.5)
    yticklabels = []
    for ytick in ax.get_yticks():
        if ytick < 0:
            val = "$-${:.0f}".format(round(abs(ytick)))
        else:
            val = "{:.0f}".format(round(abs(ytick)) - 1)
        yticklabels.append(val)
    ax.set_yticklabels(yticklabels)

    # Colorbar (has at most 10 labels)
    colorbar.set_label("$N_{\\rm realisations}$")
    if discretise:
        colorbar.locator = FixedLocator(np.arange(1, 
            len(colorbar_labels) + 1) * heatmap.max()/len(colorbar_labels) \
                - heatmap.max()/(2.*len(colorbar_labels)))
        colorbar.formatter = FixedFormatter(["{:.0f}".format(np.floor(_)) \
            for _ in colorbar_labels])
    
    # Remove ticks from colorbar
    colorbar.ax.get_children()[1].set_tick_params(width=0)
    colorbar.update_ticks()

    # Axis labels and limits.
    ax.set_xlabel("$N_{\\rm clusters}$")
    ax.set_ylabel("$\Delta{}N_{\\rm clusters}$")
    ax.set_xlim(extent[0], extent[1] + 1)
    ax.set_ylim(extent[2], extent[3] + 1)

    fig.tight_layout()
    plt.draw()

    # Remember: save this with fig.savefig(filename, bbox_inches="tight")
    return fig
