#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
Plotting functions.
'''


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.colors import to_hex
from matplotlib.collections import LineCollection
from matplotlib_venn import venn2, venn3, venn2_circles
from scipy.stats import gaussian_kde, linregress, f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd


RC = {'axes.labelsize': 14, 'axes.titlesize': 14,
      'xtick.labelsize': 14, 'ytick.labelsize': 14,
      'lines.color': 'black', 'lines.linestyle': 'dashed',#'solid',
      'lines.linewidth': 1,
      'savefig.bbox': 'tight', 'savefig.dpi': 600,
      'figure.figsize': (4, 4)}


def round_value(value, figures):
    return '{:g}'.format(float('{:.{p}g}'.format(value, p=figures)))


def scatter(data, x, y, kde=False, hue=None, hue_order=None, hue_sort=False,
            color='#023eff', palette=None, cmap='seismic', edgecolor=None,
            alpha=1, size=50, figsize=None, xlim=None, ylim=None, xticks=None,
            yticks=None, vlines=None, hlines=None, diagonal=False,
            regline=False, regline_hue=True,
            xlabel=None, ylabel=None, outfile=None):
    plt.style.use('default')
    plt.rcParams.update(RC)

    if hue and hue_order and hue_sort:
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
            data[x], data[y], c=kde, cmap=cmap, edgecolor=edgecolor,
            alpha=alpha, s=size)
    else:
        sns.scatterplot(
            data, x=x, y=y, hue=hue, hue_order=hue_order, color=color,
            palette=palette, edgecolor=edgecolor, alpha=alpha,
            s=size, legend=False)

    plt.xlim(xlim)
    plt.ylim(ylim)
    if xticks:
        plt.xticks(np.arange(xticks[0], xticks[1]+xticks[2], xticks[2]))
    if yticks:
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
        plt.gca().add_collection(
            LineCollection([[(xy1,)*2, (xy2,)*2]], linestyle='dashed'))

    if regline:
        if hue and regline_hue:
            for value, color in zip(data[hue].unique(), palette):
                sns.regplot(
                    data[data[hue] == value], x=x, y=y, scatter=False,
                    truncate=False, ci=None, color=color,
                    line_kws={'linewidth': 2, 'linestyle': 'solid'})
        else:
            sns.regplot(
                data, x=x, y=y, scatter=False, truncate=False, ci=None,
                color=color, line_kws={'linewidth': 2, 'linestyle': 'solid'})

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if outfile:
        plt.savefig(outfile)
    plt.show()
    plt.clf()


def kde(data, x, hue=None, hue_order=None, legend=False, color='#023eff',
        palette=None, alpha=1, fill=False, figsize=None, xlim=None, ylim=None,
        xticks=None, yticks=None, xlabel=None, ylabel=None, outfile=None):
    plt.style.use('default')
    plt.rcParams.update(RC)

    if hue and not palette:
        palette = sns.color_palette('bright', data[hue].nunique())
    if figsize:
        plt.figure(figsize=figsize)

    sns.kdeplot(data=data, x=x, hue=hue, hue_order=hue_order, legend=legend,
                color=color, palette=palette, alpha=alpha, fill=fill)

    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    if xticks:
        plt.xticks(np.arange(xticks[0], xticks[1]+xticks[2], xticks[2]))
    if yticks:
        plt.yticks(np.arange(yticks[0], yticks[1]+yticks[2], yticks[2]))

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if outfile:
        plt.savefig(outfile)
    plt.show()
    plt.clf()


def hist(data, x, hue=None, hue_order=None, bins=1000, legend=False,
         color='#023eff', palette=None, edgecolor=None, alpha=1,
         fill=True, figsize=None, xlim=None, ylim=None, xticks=None,
         yticks=None, xlabel=None, ylabel=None, outfile=None):
    plt.style.use('default')
    plt.rcParams.update(RC)

    if hue and not palette:
        palette = sns.color_palette('bright', data[hue].nunique())
    if figsize:
        plt.figure(figsize=figsize)

    sns.histplot(data=data, x=x, hue=hue, hue_order=hue_order, bins=bins,
                 legend=legend, color=color, palette=palette, alpha=alpha,
                 element='step', fill=fill, ls='solid')

    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    if xticks:
        plt.xticks(np.arange(xticks[0], xticks[1]+xticks[2], xticks[2]))
    if yticks:
        plt.yticks(np.arange(yticks[0], yticks[1]+yticks[2], yticks[2]))

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
        stat='density', common_norm=False, lw=2, ls='solid')

    plt.xlim(xlim)
    plt.ylim(0, 1)
    if xticks:
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
    if xticks:
        plt.xticks(np.arange(xticks[0], xticks[1]+xticks[2], xticks[2]))
    if yticks:
        plt.yticks(np.arange(yticks[0], yticks[1]+yticks[2], yticks[2]))
    if xticklabels:
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


def venn(sets, labels=None, colors=None, alpha=0.7, figsize=None,
         outfile=None):
    plt.style.use('default')
    plt.rcParams.update(RC)

    if figsize:
        plt.figure(figsize=figsize)
    if not colors:
        colors = sns.color_palette('bright', len(sets)).as_hex()
    if len(sets) == 2:
        vd2 = venn2(sets, labels, set_colors=colors, alpha=alpha)
        venn2_circles(sets, color='#000000', linewidth=1)
        vd2.get_label_by_id('100').set_text('')
        vd2.get_label_by_id('110').set_text('')
        vd2.get_label_by_id('010').set_text('')
    elif len(sets) == 3:
        venn3(sets, labels, set_colors=colors, alpha=alpha)

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


def anova(data, hue, groups, var, hsd=False):
    model = f_oneway(*[data[data[hue] == group][var] for group in groups])
    print('p =', round_value(model.pvalue, 2), '(ANOVA)\n')
    if hsd:
        hsd = pairwise_tukeyhsd(
            endog=data[var], groups=data[hue], alpha=0.05)
        print(*str(hsd.summary()).split('\n')[2:-1], sep='\n')
        print('\n')


def print_linregress(data, x_var, y_var):
    model = linregress(data[x_var], data[y_var])
    print('R2 =', round_value(model.rvalue**2, 2), '/',
          'Slope =', round_value(model.slope, 2), '/',
          'Count:', data.shape[0], '\n')


def print_counts(data, hue, groups):
    counts = data[hue].value_counts().reset_index()
    counts.sort_values(hue, key=np.vectorize(groups.index), inplace=True)
    print(counts, '\n')
