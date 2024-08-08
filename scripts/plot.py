#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
Plotting functions.
'''


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.collections import LineCollection
from scipy.stats import gaussian_kde


RC = {'axes.labelsize': 14, 'axes.titlesize': 14,
      'xtick.labelsize': 14, 'ytick.labelsize': 14,
      'lines.color': 'darkgrey', 'lines.linestyle': 'solid',
      'savefig.bbox': 'tight', 'savefig.dpi': 600,
      'figure.figsize': (4, 4)}


def round_value(value, figures):
    return '{:g}'.format(float('{:.{p}g}'.format(value, p=figures)))


def scatter(data, x, y, kde=False, hue=None, hue_order=None,
            color='#023eff', palette=None, cmap='seismic', edgecolor=None,
            size=50, figsize=None, xlim=None, ylim=None, xticks=None,
            yticks=None, vlines=None, hlines=None, diagonal=False,
            regline=False, xlabel=None, ylabel=None, outfile=None):
    plt.style.use('default')
    plt.rcParams.update(RC)

    if hue and hue_order:
        data.sort_values(
            hue, key=np.vectorize(hue_order.index), inplace=True)
    if hue and not palette:
        palette = sns.color_palette('bright', data[hue].nunique())
    if figsize:
        plt.figure(figsize=figsize)

    if kde:
        values = np.vstack([data[x], data[y]])
        kde = gaussian_kde(values)(values)
        plt.scatter(
            data[x], data[y], c=kde, cmap=cmap, edgecolor=edgecolor, s=size)
    else:
        sns.scatterplot(
            data, x=x, y=y, hue=hue, hue_order=hue_order, color=color,
            palette=palette, edgecolor=edgecolor, s=size, legend=False)

    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xticks(np.arange(xticks[0], xticks[1]+xticks[2], xticks[2]))
    plt.yticks(np.arange(yticks[0], yticks[1]+yticks[2], yticks[2]))
    if vlines:
        for line in vlines:
            plt.gca().axvline(x=line)
    if hlines:
        for line in hlines:
            plt.gca().axhline(y=line)

    if diagonal:
        xlim = plt.gca().get_xlim()
        ylim = plt.gca().get_ylim()
        xy1 = min([xlim[0], ylim[0]])
        xy2 = max([xlim[1], ylim[1]])
        plt.gca().add_collection(LineCollection([[(xy1,)*2, (xy2,)*2]]))
    if regline:
        sns.regplot(
            data, x=x, y=y, scatter=False, truncate=False, ci=None,
            color='black', line_kws={'linewidth': 1})

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if outfile:
        plt.savefig(outfile)
    plt.show()
    plt.clf()


def cumulative_hist(data, x, hue=None, hue_order=None, bins=1000, legend=False,
                    color='#023eff', palette=None, figsize=None, xlim=None,
                    xticks=None, xlabel=None, outfile=None):
    plt.style.use('default')
    plt.rcParams.update(RC)

    if hue and not palette:
        palette = sns.color_palette('bright', data[hue].nunique())
    if figsize:
        plt.figure(figsize=figsize)

    sns.histplot(
        data, x=x, hue=hue, hue_order=hue_order, color=color, palette=palette,
        legend=legend, bins=bins, cumulative=True, element='step', fill=False,
        stat='density', common_norm=False, lw=2)

    plt.xlim(xlim)
    plt.ylim(0, 1)
    plt.xticks(np.arange(xticks[0], xticks[1]+xticks[2], xticks[2]))
    plt.yticks(np.arange(0, 1.2, 0.2))

    plt.xlabel(xlabel)
    plt.ylabel('Cumulative fraction')
    if outfile:
        plt.savefig(outfile)
    plt.show()
    plt.clf()


def box(data, x, y, order=None, palette=None, figsize=None, xlim=None,
        ylim=None, xticks=None, yticks=None, xticklabels=None, hlines=None,
        xlabel=None, ylabel=None, outfile=None):
    plt.style.use('default')
    plt.rcParams.update(RC)

    if not palette:
        palette = sns.color_palette('bright', data[x].nunique())
    if figsize:
        plt.figure(figsize=figsize)

    sns.boxplot(
        data=data, x=x, y=y, order=order, hue=x, palette=palette,
        notch=True, showfliers=True, legend=False)

    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xticks(np.arange(xticks[0], xticks[1]+xticks[2], xticks[2]))
    plt.yticks(np.arange(yticks[0], yticks[1]+yticks[2], yticks[2]))
    plt.gca().set_xticklabels(xticklabels)
    if hlines:
        for line in hlines:
            plt.gca().axhline(y=line)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if outfile:
        plt.savefig(outfile)
    plt.show()
    plt.clf()


def heatmap(matrix, vmin=None, vmax=None, cmap='coolwarm', cbar=True,
            cbar_label='', figsize=None, xticklabels=None, yticklabels=None,
            xlabel=None, ylabel=None, outfile=None):
    plt.style.use('default')
    plt.rcParams.update(RC)

    if figsize:
        plt.figure(figsize=figsize)

    ax = sns.heatmap(matrix, vmin=vmin, vmax=vmax, cmap=cmap, cbar=cbar,
                     linewidths=1, clip_on=False)

    if cbar:
        ax.collections[0].colorbar.set_label(cbar_label, fontsize=14)
        ax.collections[0].colorbar.ax.tick_params(labelsize=14)

    plt.gca().set_xticklabels(xticklabels)
    plt.gca().set_yticklabels(yticklabels, rotation=360)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if outfile:
        plt.savefig(outfile)
    plt.show()
    plt.clf()
