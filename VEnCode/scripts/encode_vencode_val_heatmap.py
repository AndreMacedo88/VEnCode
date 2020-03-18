import os

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from VEnCode.scripts import encode_vencode_val_variables as var


def get_heatmap(data, path_out=None, order_by=None, one_at_time=False, label_size=7):

    def heatmap(df, lbl_size=label_size):
        fig_size = (15, 15)
        plt.figure(figsize=fig_size, dpi=600)

        colormap = plt.cm.get_cmap("plasma_r")
        # colormap = truncate_colormap(colormap, 0.03, 1)  # if you need to use only part of the colormap
        colormap.set_under('w')

        plot = plt.pcolormesh(df, cmap=colormap, vmin=0.0001)
        plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
        plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
        plt.tick_params(
            axis='x',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            labelsize=lbl_size  # label number size is reduced to fit
        )
        plt.tick_params(
            axis="y",
            which="both",
            labelsize=lbl_size  # label number size is reduced to fit
        )
        cb = plt.colorbar(plot, shrink=0.7)
        cb.ax.tick_params(labelsize=20)
        plt.axes().set_aspect("auto")
        plt.tight_layout()
        if path_out:
            plt.savefig(path_out, transparent=True, dpi="figure")
        else:
            plt.savefig(transparent=True, dpi="figure")
        plt.close()

    def heatmap_loop(df):
        fig_size = (10, 6)
        fig = plt.figure(figsize=fig_size, dpi=600)
        ax1 = fig.add_subplot(111)
        colormap = plt.cm.get_cmap("viridis_r")
        colormap.set_under('w')
        plot = ax1.pcolormesh(df, cmap=colormap, vmin=0.0001)
        plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
        plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
        plt.tick_params(
            axis='x',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            labelsize=7  # label number size is reduced to fit
        )
        plt.tick_params(
            axis="y",
            which="both",
            labelsize=7  # label number size is reduced to fit
        )
        plt.axes().set_aspect("equal")
        plt.tight_layout()
        cbaxes = fig.add_axes([0.3, 0.8, 0.6, 0.03])
        plt.colorbar(plot, cax=cbaxes, orientation="horizontal")
        if path_out:
            plt.savefig(path_out, transparent=True, dpi="figure")
        else:
            plt.savefig(transparent=True, dpi="figure")
        plt.close()

    if order_by is not None:
        data = data.reindex(order_by[0])
        data = data[order_by[1]]
    if one_at_time:
        for idx in data.index:
            data_temp = data.loc[idx]
            data_temp = data_temp.to_frame()
            heatmap_loop(data_temp.T)
    else:
        heatmap(data, lbl_size=label_size)


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


data_path_parent = "D:/Utilizador HDD/OneDrive - Nova Medical School Faculdade de Ciências Médicas da UNL/1-Research/" \
                   "3-Vencode/Fantom5/Files/Validation_files/ENCODE/"
data_name = "CAGE VEn in ENCODE DNase matrix enhancers_200bp_random-normalized.csv"
data_path = os.path.join(data_path_parent, data_name)
figure_path = os.path.join(data_path_parent, "ENCODE enh VEn validation heatmap_matching_200bp_random.png")

data_file = pd.read_csv(data_path, sep=";", engine="python", index_col=0)

y_axis = var.index_enhancer_matching[::-1]  # we reverse the list because plot axis get reversed during procedure
x_axis = var.columns_enhancer_matching
get_heatmap(data_file, figure_path, order_by=(y_axis, x_axis),
            one_at_time=False, label_size=17)
