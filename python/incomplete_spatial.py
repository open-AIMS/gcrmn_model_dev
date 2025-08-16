
# + tags=["parameters"]
# declare a list tasks whose products you want to use as inputs
upstream = None
# -

import multiprocessing

## ---- python incomplete spatial libraries
import pandas as pd
import numpy as np
import xarray as xr
import pymc as pm
import arviz as az
import preliz as pz
import pymc_bart as pmb
#from pymc_extras.model_builder import ModelBuilder
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pickle
import cloudpickle
from matplotlib.ticker import FuncFormatter
## ----end

def incomplete_spatial_global_parameters_(product):
    ## ---- python incomplete spatial global parameters
    # Assign a variable to the global namespace
    data_path = "../data/"
    output_path = "../output/"
    # Create a dictionary of paths
    paths = {
        "data_path": data_path,
        "synthetic_path": f"{data_path}synthetic/",
        "fig_path": f"{output_path}figures/"
    }
    # Ensure the directories exist
    for path in paths.values():
        if not os.path.exists(path):
            os.makedirs(path)
    ## ----end
    pickle.dump(paths, open(product, "wb"))

def incomplete_spatial_read_all_reefs_data_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    ## ---- python incomplete spatial read all reefs data
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ----end
    pickle.dump(benthos_reefs_sf, open(product, "wb"))

def incomplete_spatial_read_all_reefs_data_show_(product, upstream):
    benthos_reefs_sf = pickle.load(open(upstream["incomplete_spatial_read_all_reefs_data_"], "rb"))
    ## ---- python incomplete spatial read all reefs data show
    pd.set_option("display.max_columns", None)  # Show all columns
    benthos_reefs_sf 
    pd.reset_option("display.max_columns")
    ## ----end
    pickle.dump("done", open(product, "wb"))

def incomplete_spatial_read_all_temporal_summary_(product, upstream):
    benthos_reefs_sf = pickle.load(open(upstream["incomplete_spatial_read_all_reefs_data_"], "rb"))
    ## ---- python incomplete spatial all reefs temporal summary
    benthos_reefs_temporal_summary = (
        benthos_reefs_sf.groupby("Year").agg(
            Mean=("HCC", "mean"),
            Median=("HCC", "median"),
            SD=("HCC", "std"),
            Lower=("HCC", lambda x: x.quantile(0.025)),
            Upper=("HCC", lambda x: x.quantile(0.975))
        ).reset_index()
    )
    ## ----end
    pickle.dump(benthos_reefs_temporal_summary, open(product, "wb"))

def incomplete_spatial_read_all_temporal_summary_plot_(product, upstream):
    benthos_reefs_temporal_summary = pickle.load(open(upstream["incomplete_spatial_read_all_temporal_summary_"], "rb"))
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    ## ---- python incomplete spatial all reefs temporal summary plot
    # Create the plot
    plt.figure(figsize=(10, 6))
    # Add ribbon (shaded area)
    plt.fill_between(
        benthos_reefs_temporal_summary["Year"],
        benthos_reefs_temporal_summary["Lower"],
        benthos_reefs_temporal_summary["Upper"],
        alpha=0.2,
        label="Confidence Interval"
    )
    # Add mean line
    plt.plot(
        benthos_reefs_temporal_summary["Year"],
        benthos_reefs_temporal_summary["Mean"],
        label="Mean",
        color="blue"
    )
    # Add median line
    plt.plot(
        benthos_reefs_temporal_summary["Year"],
        benthos_reefs_temporal_summary["Median"],
        label="Median",
        color="orange"
    )
    # Add labels and theme
    plt.xlabel("Year")
    plt.ylabel("Value")
    plt.title("Temporal Summary")
    plt.legend()
    plt.grid(True)
    # plt.show()
    # Save the plot as a PNG file
    plt.savefig(f"{paths['fig_path']}Python_all_temporal_summary_plot2.png", dpi=300, bbox_inches="tight")
    ## ----end
    plt.savefig(product, dpi=72, bbox_inches="tight")

## All Northern reefs --------------------------------------------------

def incomplete_spatial_read_sampled_reefs_data_10_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    ## ---- python read sampled reefs data 10
    benthos_fixed_locs_obs_10 = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_pts_northern.csv")
    # Back transform HCC
    benthos_fixed_locs_obs_10['HCC'] = 1 / (1 + np.exp(-benthos_fixed_locs_obs_10['HCC']))
    benthos_fixed_locs_obs_10[['Longitude', 'Latitude']] = (
      benthos_fixed_locs_obs_10['geometry']
        .str.strip('c()')
        .str.split(',', expand = True)
    .    astype(float)
    )
    benthos_fixed_locs_obs_10
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_10, open(product, "wb"))

def read_sampled_reefs_data_10_plot_1_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_10 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_10_"], "rb"))
    ## ---- python read sampled reefs data 10 plot 1
    # Create a FacetGrid for faceting by Year
    g = sns.FacetGrid(benthos_fixed_locs_obs_10, col="Year", col_wrap=4, height=4, sharex=False, sharey=False)

    # Map a scatterplot to the grid (use `fill` and `color` aesthetics)
    g.map_dataframe(
    sns.scatterplot,
      x="Longitude", y="Latitude", hue="HCC", palette="viridis",
      edgecolor="none"
    )
    # Ensure the x and y axes are on the same scale
    for ax in g.axes.flat:
      ax.set_aspect('equal', adjustable='datalim')
    # Add a color bar for the hue (HCC)
    #g.add_legend()
    norm = mpl.colors.Normalize(vmin=benthos_fixed_locs_obs_10["HCC"].min(), vmax=benthos_fixed_locs_obs_10["HCC"].max())
    sm = mpl.cm.ScalarMappable(cmap="viridis", norm=norm)
    sm.set_array([])  # Required for ScalarMappable
    # Add the color bar to the figure
    g.fig.colorbar(sm, ax=g.axes, orientation="vertical", label="HCC")
    # Apply a clean theme (similar to theme_bw in ggplot2)
    sns.set_theme(style="whitegrid")
    plt.savefig(f"{paths['fig_path']}Python_sampled_reefs_10_plot_1.png", dpi=300, bbox_inches="tight")
    ## ----end
    pickle.dump("done", open(product, "wb"))

def read_sampled_reefs_data_10_plot_2_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_10 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_10_"], "rb"))
    ## ---- python read sampled reefs data 10 plot 2
    benthos_fixed_locs_obs_10 = benthos_fixed_locs_obs_10.assign(
      fYear = benthos_fixed_locs_obs_10['Year'].astype('category'),
      Reef = benthos_fixed_locs_obs_10['Reef'].astype('category'),
      cover = benthos_fixed_locs_obs_10['HCC']
    )
    benthos_reefs_temporal_summary_10 = (
        benthos_fixed_locs_obs_10.groupby("Year").agg(
            mean = ("cover", "mean"),
            median = ("cover", "median"),
            sd = ("cover", "std"),
            lower = ("cover", lambda x: x.quantile(0.025)),
            upper = ("cover", lambda x: x.quantile(0.975))
        ).reset_index()
    )
    # Initialize the plot
    plt.figure(figsize=(10, 6))
    # Plot the ribbon (shaded area between lower and upper)
    plt.fill_between(benthos_reefs_temporal_summary_10["Year"], benthos_reefs_temporal_summary_10["lower"], benthos_reefs_temporal_summary_10["upper"], color="gray", alpha=0.3, label="Quantile range (lower-upper)")

    # Plot the mean line
    sns.lineplot(x="Year", y="mean", data=benthos_reefs_temporal_summary_10, label="Mean", color="blue")

    # Plot the median line
    sns.lineplot(x="Year", y="median", data=benthos_reefs_temporal_summary_10, label="Median", color="orange")

    # Add labels and legend
    plt.xlabel("Year")
    plt.ylabel("Value")
    plt.title("Mean, Median, and Quantile range Over Years")
    plt.legend()
    plt.savefig(f"{paths['fig_path']}Python_sampled_reefs_10_plot_2.png", dpi=300, bbox_inches="tight")
    ## ----end
    pickle.dump(benthos_reefs_temporal_summary_10, open(product, "wb"))

## Northern sampled reefs --------------------------------------------------

def incomplete_spatial_read_sampled_reefs_data_12_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    ## ---- python read sampled reefs data 12
    benthos_fixed_locs_obs_12 = pd.read_csv(f"{paths['data_path']}synthetic/benthos_fixed_locs_northern_obs_disturb.csv")
    # Back transform HCC
    benthos_fixed_locs_obs_12['HCC'] = benthos_fixed_locs_obs_12['HCC']/100
    benthos_fixed_locs_obs_12
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_12, open(product, "wb"))

def read_sampled_reefs_data_12_plot_1_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_12 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_12_"], "rb"))
    ## ---- python read sampled reefs data 12 plot 1
    # Create a FacetGrid for faceting by Year
    g = sns.FacetGrid(benthos_fixed_locs_obs_12, col="Year", col_wrap=4, height=4, sharex=False, sharey=False)

    # Map a scatterplot to the grid (use `fill` and `color` aesthetics)
    g.map_dataframe(
    sns.scatterplot,
      x="Longitude", y="Latitude", hue="HCC", palette="viridis",
      edgecolor="none"
    )
    # Ensure the x and y axes are on the same scale
    for ax in g.axes.flat:
      ax.set_aspect('equal', adjustable='datalim')
    # Add a color bar for the hue (HCC)
    #g.add_legend()
    norm = mpl.colors.Normalize(vmin=benthos_fixed_locs_obs_12["HCC"].min(), vmax=benthos_fixed_locs_obs_12["HCC"].max())
    sm = mpl.cm.ScalarMappable(cmap="viridis", norm=norm)
    sm.set_array([])  # Required for ScalarMappable
    # Add the color bar to the figure
    g.fig.colorbar(sm, ax=g.axes, orientation="vertical", label="HCC")
    # Apply a clean theme (similar to theme_bw in ggplot2)
    sns.set_theme(style="whitegrid")
    plt.savefig(f"{paths['fig_path']}Python_sampled_reefs_12_plot_1.png", dpi=300, bbox_inches="tight")
    ## ----end
    pickle.dump("done", open(product, "wb"))

def read_sampled_reefs_data_12_plot_2_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_12 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_12_"], "rb"))
    ## ---- python read sampled reefs data 12 plot 2
    benthos_fixed_locs_obs_12 = benthos_fixed_locs_obs_12.assign(
      fYear = benthos_fixed_locs_obs_12['Year'].astype('category'),
      Reef = benthos_fixed_locs_obs_12['Reef'].astype('category'),
      cover = benthos_fixed_locs_obs_12['HCC']
    )
    benthos_reefs_temporal_summary_12 = (
        benthos_fixed_locs_obs_12.groupby("Year").agg(
            mean = ("cover", "mean"),
            median = ("cover", "median"),
            sd = ("cover", "std"),
            lower = ("cover", lambda x: x.quantile(0.025)),
            upper = ("cover", lambda x: x.quantile(0.975))
        ).reset_index()
    )
    # Initialize the plot
    plt.figure(figsize=(10, 6))
    # Plot the ribbon (shaded area between lower and upper)
    plt.fill_between(benthos_reefs_temporal_summary_12["Year"], benthos_reefs_temporal_summary_12["lower"], benthos_reefs_temporal_summary_12["upper"], color="gray", alpha=0.3, label="Quantile range (lower-upper)")

    # Plot the mean line
    sns.lineplot(x="Year", y="mean", data=benthos_reefs_temporal_summary_12, label="Mean", color="blue")

    # Plot the median line
    sns.lineplot(x="Year", y="median", data=benthos_reefs_temporal_summary_12, label="Median", color="orange")

    # Add labels and legend
    plt.xlabel("Year")
    plt.ylabel("Value")
    plt.title("Mean, Median, and Quantile range Over Years")
    plt.legend()
    plt.savefig(f"{paths['fig_path']}Python_sampled_reefs_12_plot_2.png", dpi=300, bbox_inches="tight")
    ## ----end
    pickle.dump(benthos_reefs_temporal_summary_12, open(product, "wb"))

## All Southern reefs --------------------------------------------------

def incomplete_spatial_read_sampled_reefs_data_11_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    ## ---- python read sampled reefs data 11
    benthos_fixed_locs_obs_11 = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_pts_southern.csv")
    # Back transform HCC
    benthos_fixed_locs_obs_11['HCC'] = 1 / (1 + np.exp(-benthos_fixed_locs_obs_11['HCC']))
    benthos_fixed_locs_obs_11[['Longitude', 'Latitude']] = (
      benthos_fixed_locs_obs_11['geometry']
        .str.strip('c()')
        .str.split(',', expand = True)
    .    astype(float)
    )
    benthos_fixed_locs_obs_11
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_11, open(product, "wb"))

def read_sampled_reefs_data_11_plot_1_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_11 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_11_"], "rb"))
    ## ---- python read sampled reefs data 11 plot 1
    # Create a FacetGrid for faceting by Year
    g = sns.FacetGrid(benthos_fixed_locs_obs_11, col="Year", col_wrap=4, height=4, sharex=False, sharey=False)

    # Map a scatterplot to the grid (use `fill` and `color` aesthetics)
    g.map_dataframe(
    sns.scatterplot,
      x="Longitude", y="Latitude", hue="HCC", palette="viridis",
      edgecolor="none"
    )
    # Ensure the x and y axes are on the same scale
    for ax in g.axes.flat:
      ax.set_aspect('equal', adjustable='datalim')
    # Add a color bar for the hue (HCC)
    #g.add_legend()
    norm = mpl.colors.Normalize(vmin=benthos_fixed_locs_obs_11["HCC"].min(), vmax=benthos_fixed_locs_obs_11["HCC"].max())
    sm = mpl.cm.ScalarMappable(cmap="viridis", norm=norm)
    sm.set_array([])  # Required for ScalarMappable
    # Add the color bar to the figure
    g.fig.colorbar(sm, ax=g.axes, orientation="vertical", label="HCC")
    # Apply a clean theme (similar to theme_bw in ggplot2)
    sns.set_theme(style="whitegrid")
    plt.savefig(f"{paths['fig_path']}Python_sampled_reefs_11_plot_1.png", dpi=300, bbox_inches="tight")
    ## ----end
    pickle.dump("done", open(product, "wb"))

def read_sampled_reefs_data_11_plot_2_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_11 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_11_"], "rb"))
    ## ---- python read sampled reefs data 11 plot 2
    benthos_fixed_locs_obs_11 = benthos_fixed_locs_obs_11.assign(
      fYear = benthos_fixed_locs_obs_11['Year'].astype('category'),
      Reef = benthos_fixed_locs_obs_11['Reef'].astype('category'),
      cover = benthos_fixed_locs_obs_11['HCC']
    )
    benthos_reefs_temporal_summary_11 = (
        benthos_fixed_locs_obs_11.groupby("Year").agg(
            mean = ("cover", "mean"),
            median = ("cover", "median"),
            sd = ("cover", "std"),
            lower = ("cover", lambda x: x.quantile(0.025)),
            upper = ("cover", lambda x: x.quantile(0.975))
        ).reset_index()
    )
    # Initialize the plot
    plt.figure(figsize=(10, 6))
    # Plot the ribbon (shaded area between lower and upper)
    plt.fill_between(benthos_reefs_temporal_summary_11["Year"], benthos_reefs_temporal_summary_11["lower"], benthos_reefs_temporal_summary_11["upper"], color="gray", alpha=0.3, label="Quantile range (lower-upper)")

    # Plot the mean line
    sns.lineplot(x="Year", y="mean", data=benthos_reefs_temporal_summary_11, label="Mean", color="blue")

    # Plot the median line
    sns.lineplot(x="Year", y="median", data=benthos_reefs_temporal_summary_11, label="Median", color="orange")

    # Add labels and legend
    plt.xlabel("Year")
    plt.ylabel("Value")
    plt.title("Mean, Median, and Quantile range Over Years")
    plt.legend()
    plt.savefig(f"{paths['fig_path']}Python_sampled_reefs_11_plot_2.png", dpi=300, bbox_inches="tight")
    ## ----end
    pickle.dump(benthos_reefs_temporal_summary_11, open(product, "wb"))

## Southern sampled reefs --------------------------------------------------

def incomplete_spatial_read_sampled_reefs_data_13_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    ## ---- python read sampled reefs data 13
    benthos_fixed_locs_obs_13 = pd.read_csv(f"{paths['data_path']}synthetic/benthos_fixed_locs_southern_obs_disturb.csv")
    # Back transform HCC
    benthos_fixed_locs_obs_13['HCC'] = benthos_fixed_locs_obs_13['HCC']/100
    benthos_fixed_locs_obs_13
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_13, open(product, "wb"))

def read_sampled_reefs_data_13_plot_1_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_13 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_13_"], "rb"))
    ## ---- python read sampled reefs data 13 plot 1
    # Create a FacetGrid for faceting by Year
    g = sns.FacetGrid(benthos_fixed_locs_obs_13, col="Year", col_wrap=4, height=4, sharex=False, sharey=False)

    # Map a scatterplot to the grid (use `fill` and `color` aesthetics)
    g.map_dataframe(
    sns.scatterplot,
      x="Longitude", y="Latitude", hue="HCC", palette="viridis",
      edgecolor="none"
    )
    # Ensure the x and y axes are on the same scale
    for ax in g.axes.flat:
      ax.set_aspect('equal', adjustable='datalim')
    # Add a color bar for the hue (HCC)
    #g.add_legend()
    norm = mpl.colors.Normalize(vmin=benthos_fixed_locs_obs_13["HCC"].min(), vmax=benthos_fixed_locs_obs_13["HCC"].max())
    sm = mpl.cm.ScalarMappable(cmap="viridis", norm=norm)
    sm.set_array([])  # Required for ScalarMappable
    # Add the color bar to the figure
    g.fig.colorbar(sm, ax=g.axes, orientation="vertical", label="HCC")
    # Apply a clean theme (similar to theme_bw in ggplot2)
    sns.set_theme(style="whitegrid")
    plt.savefig(f"{paths['fig_path']}Python_sampled_reefs_13_plot_1.png", dpi=300, bbox_inches="tight")
    ## ----end
    pickle.dump("done", open(product, "wb"))

def read_sampled_reefs_data_13_plot_2_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_13 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_13_"], "rb"))
    ## ---- python read sampled reefs data 13 plot 2
    benthos_fixed_locs_obs_13 = benthos_fixed_locs_obs_13.assign(
      fYear = benthos_fixed_locs_obs_13['Year'].astype('category'),
      Reef = benthos_fixed_locs_obs_13['Reef'].astype('category'),
      cover = benthos_fixed_locs_obs_13['HCC']
    )
    benthos_reefs_temporal_summary_13 = (
        benthos_fixed_locs_obs_13.groupby("Year").agg(
            mean = ("cover", "mean"),
            median = ("cover", "median"),
            sd = ("cover", "std"),
            lower = ("cover", lambda x: x.quantile(0.025)),
            upper = ("cover", lambda x: x.quantile(0.975))
        ).reset_index()
    )
    # Initialize the plot
    plt.figure(figsize=(10, 6))
    # Plot the ribbon (shaded area between lower and upper)
    plt.fill_between(benthos_reefs_temporal_summary_13["Year"], benthos_reefs_temporal_summary_13["lower"], benthos_reefs_temporal_summary_13["upper"], color="gray", alpha=0.3, label="Quantile range (lower-upper)")

    # Plot the mean line
    sns.lineplot(x="Year", y="mean", data=benthos_reefs_temporal_summary_13, label="Mean", color="blue")

    # Plot the median line
    sns.lineplot(x="Year", y="median", data=benthos_reefs_temporal_summary_13, label="Median", color="orange")

    # Add labels and legend
    plt.xlabel("Year")
    plt.ylabel("Value")
    plt.title("Mean, Median, and Quantile range Over Years")
    plt.legend()
    plt.savefig(f"{paths['fig_path']}Python_sampled_reefs_13_plot_2.png", dpi=300, bbox_inches="tight")
    ## ----end
    pickle.dump(benthos_reefs_temporal_summary_13, open(product, "wb"))

## All Western reefs --------------------------------------------------

def incomplete_spatial_read_sampled_reefs_data_14_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    ## ---- python read sampled reefs data 14
    benthos_fixed_locs_obs_14 = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_pts_western.csv")
    # Back transform HCC
    benthos_fixed_locs_obs_14['HCC'] = 1 / (1 + np.exp(-benthos_fixed_locs_obs_14['HCC']))
    benthos_fixed_locs_obs_14[['Longitude', 'Latitude']] = (
      benthos_fixed_locs_obs_14['geometry']
        .str.strip('c()')
        .str.split(',', expand = True)
    .    astype(float)
    )
    benthos_fixed_locs_obs_14
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_14, open(product, "wb"))

def read_sampled_reefs_data_14_plot_1_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_14 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_14_"], "rb"))
    ## ---- python read sampled reefs data 14 plot 1
    # Create a FacetGrid for faceting by Year
    g = sns.FacetGrid(benthos_fixed_locs_obs_14, col="Year", col_wrap=4, height=4, sharex=False, sharey=False)

    # Map a scatterplot to the grid (use `fill` and `color` aesthetics)
    g.map_dataframe(
    sns.scatterplot,
      x="Longitude", y="Latitude", hue="HCC", palette="viridis",
      edgecolor="none"
    )
    # Ensure the x and y axes are on the same scale
    for ax in g.axes.flat:
      ax.set_aspect('equal', adjustable='datalim')
    # Add a color bar for the hue (HCC)
    #g.add_legend()
    norm = mpl.colors.Normalize(vmin=benthos_fixed_locs_obs_14["HCC"].min(), vmax=benthos_fixed_locs_obs_14["HCC"].max())
    sm = mpl.cm.ScalarMappable(cmap="viridis", norm=norm)
    sm.set_array([])  # Required for ScalarMappable
    # Add the color bar to the figure
    g.fig.colorbar(sm, ax=g.axes, orientation="vertical", label="HCC")
    # Apply a clean theme (similar to theme_bw in ggplot2)
    sns.set_theme(style="whitegrid")
    plt.savefig(f"{paths['fig_path']}Python_sampled_reefs_14_plot_1.png", dpi=300, bbox_inches="tight")
    ## ----end
    pickle.dump("done", open(product, "wb"))

def read_sampled_reefs_data_14_plot_2_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_14 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_14_"], "rb"))
    ## ---- python read sampled reefs data 14 plot 2
    benthos_fixed_locs_obs_14 = benthos_fixed_locs_obs_14.assign(
      fYear = benthos_fixed_locs_obs_14['Year'].astype('category'),
      Reef = benthos_fixed_locs_obs_14['Reef'].astype('category'),
      cover = benthos_fixed_locs_obs_14['HCC']
    )
    benthos_reefs_temporal_summary_14 = (
        benthos_fixed_locs_obs_14.groupby("Year").agg(
            mean = ("cover", "mean"),
            median = ("cover", "median"),
            sd = ("cover", "std"),
            lower = ("cover", lambda x: x.quantile(0.025)),
            upper = ("cover", lambda x: x.quantile(0.975))
        ).reset_index()
    )
    # Initialize the plot
    plt.figure(figsize=(10, 6))
    # Plot the ribbon (shaded area between lower and upper)
    plt.fill_between(benthos_reefs_temporal_summary_14["Year"], benthos_reefs_temporal_summary_14["lower"], benthos_reefs_temporal_summary_14["upper"], color="gray", alpha=0.3, label="Quantile range (lower-upper)")

    # Plot the mean line
    sns.lineplot(x="Year", y="mean", data=benthos_reefs_temporal_summary_14, label="Mean", color="blue")

    # Plot the median line
    sns.lineplot(x="Year", y="median", data=benthos_reefs_temporal_summary_14, label="Median", color="orange")

    # Add labels and legend
    plt.xlabel("Year")
    plt.ylabel("Value")
    plt.title("Mean, Median, and Quantile range Over Years")
    plt.legend()
    plt.savefig(f"{paths['fig_path']}Python_sampled_reefs_14_plot_2.png", dpi=300, bbox_inches="tight")
    ## ----end
    pickle.dump(benthos_reefs_temporal_summary_14, open(product, "wb"))

## Western sampled reefs --------------------------------------------------

def incomplete_spatial_read_sampled_reefs_data_16_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    ## ---- python read sampled reefs data 16
    benthos_fixed_locs_obs_16 = pd.read_csv(f"{paths['data_path']}synthetic/benthos_fixed_locs_western_obs_disturb.csv")
    # Back transform HCC
    benthos_fixed_locs_obs_16['HCC'] = benthos_fixed_locs_obs_16['HCC']/100
    benthos_fixed_locs_obs_16
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_16, open(product, "wb"))

def read_sampled_reefs_data_16_plot_1_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_16 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_16_"], "rb"))
    ## ---- python read sampled reefs data 16 plot 1
    # Create a FacetGrid for faceting by Year
    g = sns.FacetGrid(benthos_fixed_locs_obs_16, col="Year", col_wrap=4, height=4, sharex=False, sharey=False)

    # Map a scatterplot to the grid (use `fill` and `color` aesthetics)
    g.map_dataframe(
    sns.scatterplot,
      x="Longitude", y="Latitude", hue="HCC", palette="viridis",
      edgecolor="none"
    )
    # Ensure the x and y axes are on the same scale
    for ax in g.axes.flat:
      ax.set_aspect('equal', adjustable='datalim')
    # Add a color bar for the hue (HCC)
    #g.add_legend()
    norm = mpl.colors.Normalize(vmin=benthos_fixed_locs_obs_16["HCC"].min(), vmax=benthos_fixed_locs_obs_16["HCC"].max())
    sm = mpl.cm.ScalarMappable(cmap="viridis", norm=norm)
    sm.set_array([])  # Required for ScalarMappable
    # Add the color bar to the figure
    g.fig.colorbar(sm, ax=g.axes, orientation="vertical", label="HCC")
    # Apply a clean theme (similar to theme_bw in ggplot2)
    sns.set_theme(style="whitegrid")
    plt.savefig(f"{paths['fig_path']}Python_sampled_reefs_16_plot_1.png", dpi=300, bbox_inches="tight")
    ## ----end
    pickle.dump("done", open(product, "wb"))

def read_sampled_reefs_data_16_plot_2_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_16 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_16_"], "rb"))
    ## ---- python read sampled reefs data 16 plot 2
    benthos_fixed_locs_obs_16 = benthos_fixed_locs_obs_16.assign(
      fYear = benthos_fixed_locs_obs_16['Year'].astype('category'),
      Reef = benthos_fixed_locs_obs_16['Reef'].astype('category'),
      cover = benthos_fixed_locs_obs_16['HCC']
    )
    benthos_reefs_temporal_summary_16 = (
        benthos_fixed_locs_obs_16.groupby("Year").agg(
            mean = ("cover", "mean"),
            median = ("cover", "median"),
            sd = ("cover", "std"),
            lower = ("cover", lambda x: x.quantile(0.025)),
            upper = ("cover", lambda x: x.quantile(0.975))
        ).reset_index()
    )
    # Initialize the plot
    plt.figure(figsize=(10, 6))
    # Plot the ribbon (shaded area between lower and upper)
    plt.fill_between(benthos_reefs_temporal_summary_16["Year"], benthos_reefs_temporal_summary_16["lower"], benthos_reefs_temporal_summary_16["upper"], color="gray", alpha=0.3, label="Quantile range (lower-upper)")

    # Plot the mean line
    sns.lineplot(x="Year", y="mean", data=benthos_reefs_temporal_summary_16, label="Mean", color="blue")

    # Plot the median line
    sns.lineplot(x="Year", y="median", data=benthos_reefs_temporal_summary_16, label="Median", color="orange")

    # Add labels and legend
    plt.xlabel("Year")
    plt.ylabel("Value")
    plt.title("Mean, Median, and Quantile range Over Years")
    plt.legend()
    plt.savefig(f"{paths['fig_path']}Python_sampled_reefs_16_plot_2.png", dpi=300, bbox_inches="tight")
    ## ----end
    pickle.dump(benthos_reefs_temporal_summary_16, open(product, "wb"))

## All Eastern reefs --------------------------------------------------

def incomplete_spatial_read_sampled_reefs_data_15_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    ## ---- python read sampled reefs data 15
    benthos_fixed_locs_obs_15 = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_pts_eastern.csv")
    # Back transform HCC
    benthos_fixed_locs_obs_15['HCC'] = 1 / (1 + np.exp(-benthos_fixed_locs_obs_15['HCC']))
    benthos_fixed_locs_obs_15[['Longitude', 'Latitude']] = (
      benthos_fixed_locs_obs_15['geometry']
        .str.strip('c()')
        .str.split(',', expand = True)
    .    astype(float)
    )
    benthos_fixed_locs_obs_15
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_15, open(product, "wb"))

def read_sampled_reefs_data_15_plot_1_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_15 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_15_"], "rb"))
    ## ---- python read sampled reefs data 15 plot 1
    # Create a FacetGrid for faceting by Year
    g = sns.FacetGrid(benthos_fixed_locs_obs_15, col="Year", col_wrap=4, height=4, sharex=False, sharey=False)

    # Map a scatterplot to the grid (use `fill` and `color` aesthetics)
    g.map_dataframe(
    sns.scatterplot,
      x="Longitude", y="Latitude", hue="HCC", palette="viridis",
      edgecolor="none"
    )
    # Ensure the x and y axes are on the same scale
    for ax in g.axes.flat:
      ax.set_aspect('equal', adjustable='datalim')
    # Add a color bar for the hue (HCC)
    #g.add_legend()
    norm = mpl.colors.Normalize(vmin=benthos_fixed_locs_obs_15["HCC"].min(), vmax=benthos_fixed_locs_obs_15["HCC"].max())
    sm = mpl.cm.ScalarMappable(cmap="viridis", norm=norm)
    sm.set_array([])  # Required for ScalarMappable
    # Add the color bar to the figure
    g.fig.colorbar(sm, ax=g.axes, orientation="vertical", label="HCC")
    # Apply a clean theme (similar to theme_bw in ggplot2)
    sns.set_theme(style="whitegrid")
    plt.savefig(f"{paths['fig_path']}Python_sampled_reefs_15_plot_1.png", dpi=300, bbox_inches="tight")
    ## ----end
    pickle.dump("done", open(product, "wb"))

def read_sampled_reefs_data_15_plot_2_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_15 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_15_"], "rb"))
    ## ---- python read sampled reefs data 15 plot 2
    benthos_fixed_locs_obs_15 = benthos_fixed_locs_obs_15.assign(
      fYear = benthos_fixed_locs_obs_15['Year'].astype('category'),
      Reef = benthos_fixed_locs_obs_15['Reef'].astype('category'),
      cover = benthos_fixed_locs_obs_15['HCC']
    )
    benthos_reefs_temporal_summary_15 = (
        benthos_fixed_locs_obs_15.groupby("Year").agg(
            mean = ("cover", "mean"),
            median = ("cover", "median"),
            sd = ("cover", "std"),
            lower = ("cover", lambda x: x.quantile(0.025)),
            upper = ("cover", lambda x: x.quantile(0.975))
        ).reset_index()
    )
    # Initialize the plot
    plt.figure(figsize=(10, 6))
    # Plot the ribbon (shaded area between lower and upper)
    plt.fill_between(benthos_reefs_temporal_summary_15["Year"], benthos_reefs_temporal_summary_15["lower"], benthos_reefs_temporal_summary_15["upper"], color="gray", alpha=0.3, label="Quantile range (lower-upper)")

    # Plot the mean line
    sns.lineplot(x="Year", y="mean", data=benthos_reefs_temporal_summary_15, label="Mean", color="blue")

    # Plot the median line
    sns.lineplot(x="Year", y="median", data=benthos_reefs_temporal_summary_15, label="Median", color="orange")

    # Add labels and legend
    plt.xlabel("Year")
    plt.ylabel("Value")
    plt.title("Mean, Median, and Quantile range Over Years")
    plt.legend()
    plt.savefig(f"{paths['fig_path']}Python_sampled_reefs_15_plot_2.png", dpi=300, bbox_inches="tight")
    ## ----end
    pickle.dump(benthos_reefs_temporal_summary_15, open(product, "wb"))

## Eastern sampled reefs --------------------------------------------------

def incomplete_spatial_read_sampled_reefs_data_17_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    ## ---- python read sampled reefs data 17
    benthos_fixed_locs_obs_17 = pd.read_csv(f"{paths['data_path']}synthetic/benthos_fixed_locs_eastern_obs_disturb.csv")
    # Back transform HCC
    benthos_fixed_locs_obs_17['HCC'] = benthos_fixed_locs_obs_17['HCC']/100
    benthos_fixed_locs_obs_17
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_17, open(product, "wb"))

def read_sampled_reefs_data_17_plot_1_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_17 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_17_"], "rb"))
    ## ---- python read sampled reefs data 17 plot 1
    # Create a FacetGrid for faceting by Year
    g = sns.FacetGrid(benthos_fixed_locs_obs_17, col="Year", col_wrap=4, height=4, sharex=False, sharey=False)

    # Map a scatterplot to the grid (use `fill` and `color` aesthetics)
    g.map_dataframe(
    sns.scatterplot,
      x="Longitude", y="Latitude", hue="HCC", palette="viridis",
      edgecolor="none"
    )
    # Ensure the x and y axes are on the same scale
    for ax in g.axes.flat:
      ax.set_aspect('equal', adjustable='datalim')
    # Add a color bar for the hue (HCC)
    #g.add_legend()
    norm = mpl.colors.Normalize(vmin=benthos_fixed_locs_obs_17["HCC"].min(), vmax=benthos_fixed_locs_obs_17["HCC"].max())
    sm = mpl.cm.ScalarMappable(cmap="viridis", norm=norm)
    sm.set_array([])  # Required for ScalarMappable
    # Add the color bar to the figure
    g.fig.colorbar(sm, ax=g.axes, orientation="vertical", label="HCC")
    # Apply a clean theme (similar to theme_bw in ggplot2)
    sns.set_theme(style="whitegrid")
    plt.savefig(f"{paths['fig_path']}Python_sampled_reefs_17_plot_1.png", dpi=300, bbox_inches="tight")
    ## ----end
    pickle.dump("done", open(product, "wb"))

def read_sampled_reefs_data_17_plot_2_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_17 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_17_"], "rb"))
    ## ---- python read sampled reefs data 17 plot 2
    benthos_fixed_locs_obs_17 = benthos_fixed_locs_obs_17.assign(
      fYear = benthos_fixed_locs_obs_17['Year'].astype('category'),
      Reef = benthos_fixed_locs_obs_17['Reef'].astype('category'),
      cover = benthos_fixed_locs_obs_17['HCC']
    )
    benthos_reefs_temporal_summary_17 = (
        benthos_fixed_locs_obs_17.groupby("Year").agg(
            mean = ("cover", "mean"),
            median = ("cover", "median"),
            sd = ("cover", "std"),
            lower = ("cover", lambda x: x.quantile(0.025)),
            upper = ("cover", lambda x: x.quantile(0.975))
        ).reset_index()
    )
    # Initialize the plot
    plt.figure(figsize=(10, 6))
    # Plot the ribbon (shaded area between lower and upper)
    plt.fill_between(benthos_reefs_temporal_summary_17["Year"], benthos_reefs_temporal_summary_17["lower"], benthos_reefs_temporal_summary_17["upper"], color="gray", alpha=0.3, label="Quantile range (lower-upper)")

    # Plot the mean line
    sns.lineplot(x="Year", y="mean", data=benthos_reefs_temporal_summary_17, label="Mean", color="blue")

    # Plot the median line
    sns.lineplot(x="Year", y="median", data=benthos_reefs_temporal_summary_17, label="Median", color="orange")

    # Add labels and legend
    plt.xlabel("Year")
    plt.ylabel("Value")
    plt.title("Mean, Median, and Quantile range Over Years")
    plt.legend()
    plt.savefig(f"{paths['fig_path']}Python_sampled_reefs_17_plot_2.png", dpi=300, bbox_inches="tight")
    ## ----end
    pickle.dump(benthos_reefs_temporal_summary_17, open(product, "wb"))

## Data transformations ============================================

## All Northern reefs --------------------------------------------------

def incomplete_spatial_data_prep_10_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_10 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_10_"], "rb"))
    ## ---- python sampled data prep 10
    # Back transform HCC
    benthos_fixed_locs_obs_10 = benthos_fixed_locs_obs_10.assign(
      fYear = benthos_fixed_locs_obs_10['Year'].astype('category'),
      Reef = benthos_fixed_locs_obs_10['Reef'].astype('category'),
      cover = benthos_fixed_locs_obs_10['HCC']
    )
    benthos_fixed_locs_obs_10
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_10, open(product, "wb"))

def incomplete_spatial_newdata_10_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_10 = pickle.load(open(upstream["incomplete_spatial_data_prep_10_"], "rb"))
    ## ---- python sampled newdata 10
    newdata_10 = pd.MultiIndex.from_product(
      [benthos_fixed_locs_obs_10["fYear"].unique()],
      names=["fYear"]
    ).to_frame(index=False)
    newdata_10
    ## ----end
    pickle.dump(newdata_10, open(product, "wb"))

def incomplete_spatial_newdata_10b_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_10 = pickle.load(open(upstream["incomplete_spatial_data_prep_10_"], "rb"))
    ## ---- python sampled newdata 10b
    benthos_fixed_locs_obs_10[['Longitude', 'Latitude']] = (
      benthos_fixed_locs_obs_10['geometry']
      .str.strip('c()')
      .str.split(',', expand = True)
      .astype(float)
    )
    newdata_10b = (
      benthos_fixed_locs_obs_10.groupby(["Year", "Reef"], as_index=False, observed=True).agg(
        HCC = ("HCC", "mean"),
        CYC = ("CYC", "mean"),
        DHW = ("DHW", "mean"),
        OTHER = ("OTHER", "mean"),
        Longitude = ("Longitude", "mean"),
        Latitude = ("Latitude", "mean")
      ).reset_index()
    )
    newdata_10b
    ## ----end
    pickle.dump(newdata_10b, open(product, "wb"))


## All Southern reefs --------------------------------------------------

def incomplete_spatial_data_prep_11_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_11 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_11_"], "rb"))
    ## ---- python sampled data prep 11
    # Back transform HCC
    benthos_fixed_locs_obs_11 = benthos_fixed_locs_obs_11.assign(
      fYear = benthos_fixed_locs_obs_11['Year'].astype('category'),
      Reef = benthos_fixed_locs_obs_11['Reef'].astype('category'),
      cover = benthos_fixed_locs_obs_11['HCC']
    )
    benthos_fixed_locs_obs_11
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_11, open(product, "wb"))

def incomplete_spatial_newdata_11_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_11 = pickle.load(open(upstream["incomplete_spatial_data_prep_11_"], "rb"))
    ## ---- python sampled newdata 11
    newdata_11 = pd.MultiIndex.from_product(
      [benthos_fixed_locs_obs_11["fYear"].unique()],
      names=["fYear"]
    ).to_frame(index=False)
    newdata_11
    ## ----end
    pickle.dump(newdata_11, open(product, "wb"))

def incomplete_spatial_newdata_11b_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_11 = pickle.load(open(upstream["incomplete_spatial_data_prep_11_"], "rb"))
    ## ---- python sampled newdata 11b
    benthos_fixed_locs_obs_11[['Longitude', 'Latitude']] = (
      benthos_fixed_locs_obs_11['geometry']
      .str.strip('c()')
      .str.split(',', expand = True)
      .astype(float)
    )
    newdata_11b = (
      benthos_fixed_locs_obs_11.groupby(["Year", "Reef"], as_index=False, observed=True).agg(
        HCC = ("HCC", "mean"),
        CYC = ("CYC", "mean"),
        DHW = ("DHW", "mean"),
        OTHER = ("OTHER", "mean"),
        Longitude = ("Longitude", "mean"),
        Latitude = ("Latitude", "mean")
      ).reset_index()
    )
    newdata_11b
    ## ----end
    pickle.dump(newdata_11b, open(product, "wb"))

## Northern sampled reefs --------------------------------------------------

def incomplete_spatial_data_prep_12_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_12 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_12_"], "rb"))
    ## ---- python sampled data prep 12
    # Back transform HCC
    benthos_fixed_locs_obs_12 = benthos_fixed_locs_obs_12.assign(
      fYear = benthos_fixed_locs_obs_12['Year'].astype('category'),
      Reef = benthos_fixed_locs_obs_12['Reef'].astype('category'),
      Site = benthos_fixed_locs_obs_12['Reef'].astype(str) + "_" +
            benthos_fixed_locs_obs_12['Site'].astype(str),
      Transect = benthos_fixed_locs_obs_12['Site'] + "_" +
            benthos_fixed_locs_obs_12['Transect'].astype(str),
      cover = benthos_fixed_locs_obs_12['HCC']
    )
    benthos_fixed_locs_obs_12
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_12, open(product, "wb"))

def incomplete_spatial_newdata_12_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_12 = pickle.load(open(upstream["incomplete_spatial_data_prep_12_"], "rb"))
    ## ---- python sampled newdata 12
    newdata_12 = pd.MultiIndex.from_product(
      [benthos_fixed_locs_obs_12["fYear"].unique()],
      names=["fYear"]
    ).to_frame(index=False)
    newdata_12
    ## ----end
    pickle.dump(newdata_12, open(product, "wb"))

def incomplete_spatial_newdata_12b_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_12 = pickle.load(open(upstream["incomplete_spatial_data_prep_12_"], "rb"))
    ## ---- python sampled newdata 12b
    benthos_fixed_locs_obs_12[['Longitude', 'Latitude']] = (
      benthos_fixed_locs_obs_12['geometry']
      .str.strip('c()')
      .str.split(',', expand = True)
      .astype(float)
    )
    newdata_12b = (
      benthos_fixed_locs_obs_12.groupby(["Year", "Reef"], as_index=False, observed=True).agg(
        HCC = ("HCC", "mean"),
        CYC = ("CYC", "mean"),
        DHW = ("DHW", "mean"),
        OTHER = ("OTHER", "mean"),
        Longitude = ("Longitude", "mean"),
        Latitude = ("Latitude", "mean")
      ).reset_index()
    )
    newdata_12b
    ## ----end
    pickle.dump(newdata_12b, open(product, "wb"))

## Southern sampled reefs --------------------------------------------------

def incomplete_spatial_data_prep_13_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_13 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_13_"], "rb"))
    ## ---- python sampled data prep 13
    # Back transform HCC
    benthos_fixed_locs_obs_13 = benthos_fixed_locs_obs_13.assign(
      fYear = benthos_fixed_locs_obs_13['Year'].astype('category'),
      Reef = benthos_fixed_locs_obs_13['Reef'].astype('category'),
      Site = benthos_fixed_locs_obs_13['Reef'].astype(str) + "_" +
            benthos_fixed_locs_obs_13['Site'].astype(str),
      Transect = benthos_fixed_locs_obs_13['Site'] + "_" +
            benthos_fixed_locs_obs_13['Transect'].astype(str),
      cover = benthos_fixed_locs_obs_13['HCC']
    )
    benthos_fixed_locs_obs_13
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_13, open(product, "wb"))

def incomplete_spatial_newdata_13_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_13 = pickle.load(open(upstream["incomplete_spatial_data_prep_13_"], "rb"))
    ## ---- python sampled newdata 13
    newdata_13 = pd.MultiIndex.from_product(
      [benthos_fixed_locs_obs_13["fYear"].unique()],
      names=["fYear"]
    ).to_frame(index=False)
    newdata_13
    ## ----end
    pickle.dump(newdata_13, open(product, "wb"))

def incomplete_spatial_newdata_13b_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_13 = pickle.load(open(upstream["incomplete_spatial_data_prep_13_"], "rb"))
    ## ---- python sampled newdata 13b
    benthos_fixed_locs_obs_13[['Longitude', 'Latitude']] = (
      benthos_fixed_locs_obs_13['geometry']
      .str.strip('c()')
      .str.split(',', expand = True)
      .astype(float)
    )
    newdata_13b = (
      benthos_fixed_locs_obs_13.groupby(["Year", "Reef"], as_index=False, observed=True).agg(
        HCC = ("HCC", "mean"),
        CYC = ("CYC", "mean"),
        DHW = ("DHW", "mean"),
        OTHER = ("OTHER", "mean"),
        Longitude = ("Longitude", "mean"),
        Latitude = ("Latitude", "mean")
      ).reset_index()
    )
    newdata_13b
    ## ----end
    pickle.dump(newdata_13b, open(product, "wb"))

## All Western reefs --------------------------------------------------

def incomplete_spatial_data_prep_14_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_14 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_14_"], "rb"))
    ## ---- python sampled data prep 14
    # Back transform HCC
    benthos_fixed_locs_obs_14 = benthos_fixed_locs_obs_14.assign(
      fYear = benthos_fixed_locs_obs_14['Year'].astype('category'),
      Reef = benthos_fixed_locs_obs_14['Reef'].astype('category'),
      cover = benthos_fixed_locs_obs_14['HCC']
    )
    benthos_fixed_locs_obs_14
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_14, open(product, "wb"))

def incomplete_spatial_newdata_14_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_14 = pickle.load(open(upstream["incomplete_spatial_data_prep_14_"], "rb"))
    ## ---- python sampled newdata 14
    newdata_14 = pd.MultiIndex.from_product(
      [benthos_fixed_locs_obs_14["fYear"].unique()],
      names=["fYear"]
    ).to_frame(index=False)
    newdata_14
    ## ----end
    pickle.dump(newdata_14, open(product, "wb"))

def incomplete_spatial_newdata_14b_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_14 = pickle.load(open(upstream["incomplete_spatial_data_prep_14_"], "rb"))
    ## ---- python sampled newdata 14b
    benthos_fixed_locs_obs_14[['Longitude', 'Latitude']] = (
      benthos_fixed_locs_obs_14['geometry']
      .str.strip('c()')
      .str.split(',', expand = True)
      .astype(float)
    )
    newdata_14b = (
      benthos_fixed_locs_obs_14.groupby(["Year", "Reef"], as_index=False, observed=True).agg(
        HCC = ("HCC", "mean"),
        CYC = ("CYC", "mean"),
        DHW = ("DHW", "mean"),
        OTHER = ("OTHER", "mean"),
        Longitude = ("Longitude", "mean"),
        Latitude = ("Latitude", "mean")
      ).reset_index()
    )
    newdata_14b
    ## ----end
    pickle.dump(newdata_14b, open(product, "wb"))

## Western sampled reefs --------------------------------------------------

def incomplete_spatial_data_prep_16_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_16 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_16_"], "rb"))
    ## ---- python sampled data prep 16
    # Back transform HCC
    benthos_fixed_locs_obs_16 = benthos_fixed_locs_obs_16.assign(
      fYear = benthos_fixed_locs_obs_16['Year'].astype('category'),
      Reef = benthos_fixed_locs_obs_16['Reef'].astype('category'),
      Site = benthos_fixed_locs_obs_16['Reef'].astype(str) + "_" +
            benthos_fixed_locs_obs_16['Site'].astype(str),
      Transect = benthos_fixed_locs_obs_16['Site'] + "_" +
            benthos_fixed_locs_obs_16['Transect'].astype(str),
      cover = benthos_fixed_locs_obs_16['HCC']
    )
    benthos_fixed_locs_obs_16
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_16, open(product, "wb"))

def incomplete_spatial_newdata_16_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_16 = pickle.load(open(upstream["incomplete_spatial_data_prep_16_"], "rb"))
    ## ---- python sampled newdata 16
    newdata_16 = pd.MultiIndex.from_product(
      [benthos_fixed_locs_obs_16["fYear"].unique()],
      names=["fYear"]
    ).to_frame(index=False)
    newdata_16
    ## ----end
    pickle.dump(newdata_16, open(product, "wb"))

def incomplete_spatial_newdata_16b_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_16 = pickle.load(open(upstream["incomplete_spatial_data_prep_16_"], "rb"))
    ## ---- python sampled newdata 16b
    benthos_fixed_locs_obs_16[['Longitude', 'Latitude']] = (
      benthos_fixed_locs_obs_16['geometry']
      .str.strip('c()')
      .str.split(',', expand = True)
      .astype(float)
    )
    newdata_16b = (
      benthos_fixed_locs_obs_16.groupby(["Year", "Reef"], as_index=False, observed=True).agg(
        HCC = ("HCC", "mean"),
        CYC = ("CYC", "mean"),
        DHW = ("DHW", "mean"),
        OTHER = ("OTHER", "mean"),
        Longitude = ("Longitude", "mean"),
        Latitude = ("Latitude", "mean")
      ).reset_index()
    )
    newdata_16b
    ## ----end
    pickle.dump(newdata_16b, open(product, "wb"))

## All Eastern reefs --------------------------------------------------

def incomplete_spatial_data_prep_15_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_15 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_15_"], "rb"))
    ## ---- python sampled data prep 15
    # Back transform HCC
    benthos_fixed_locs_obs_15 = benthos_fixed_locs_obs_15.assign(
      fYear = benthos_fixed_locs_obs_15['Year'].astype('category'),
      Reef = benthos_fixed_locs_obs_15['Reef'].astype('category'),
      cover = benthos_fixed_locs_obs_15['HCC']
    )
    benthos_fixed_locs_obs_15
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_15, open(product, "wb"))

def incomplete_spatial_newdata_15_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_15 = pickle.load(open(upstream["incomplete_spatial_data_prep_15_"], "rb"))
    ## ---- python sampled newdata 15
    newdata_15 = pd.MultiIndex.from_product(
      [benthos_fixed_locs_obs_15["fYear"].unique()],
      names=["fYear"]
    ).to_frame(index=False)
    newdata_15
    ## ----end
    pickle.dump(newdata_15, open(product, "wb"))

def incomplete_spatial_newdata_15b_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_15 = pickle.load(open(upstream["incomplete_spatial_data_prep_15_"], "rb"))
    ## ---- python sampled newdata 15b
    benthos_fixed_locs_obs_15[['Longitude', 'Latitude']] = (
      benthos_fixed_locs_obs_15['geometry']
      .str.strip('c()')
      .str.split(',', expand = True)
      .astype(float)
    )
    newdata_15b = (
      benthos_fixed_locs_obs_15.groupby(["Year", "Reef"], as_index=False, observed=True).agg(
        HCC = ("HCC", "mean"),
        CYC = ("CYC", "mean"),
        DHW = ("DHW", "mean"),
        OTHER = ("OTHER", "mean"),
        Longitude = ("Longitude", "mean"),
        Latitude = ("Latitude", "mean")
      ).reset_index()
    )
    newdata_15b
    ## ----end
    pickle.dump(newdata_15b, open(product, "wb"))

## Eastern sampled reefs --------------------------------------------------

def incomplete_spatial_data_prep_17_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_17 = pickle.load(open(upstream["incomplete_spatial_read_sampled_reefs_data_17_"], "rb"))
    ## ---- python sampled data prep 17
    # Back transform HCC
    benthos_fixed_locs_obs_17 = benthos_fixed_locs_obs_17.assign(
      fYear = benthos_fixed_locs_obs_17['Year'].astype('category'),
      Reef = benthos_fixed_locs_obs_17['Reef'].astype('category'),
      Site = benthos_fixed_locs_obs_17['Reef'].astype(str) + "_" +
            benthos_fixed_locs_obs_17['Site'].astype(str),
      Transect = benthos_fixed_locs_obs_17['Site'] + "_" +
            benthos_fixed_locs_obs_17['Transect'].astype(str),
      cover = benthos_fixed_locs_obs_17['HCC']
    )
    benthos_fixed_locs_obs_17
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_17, open(product, "wb"))

def incomplete_spatial_newdata_17_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_17 = pickle.load(open(upstream["incomplete_spatial_data_prep_17_"], "rb"))
    ## ---- python sampled newdata 17
    newdata_17 = pd.MultiIndex.from_product(
      [benthos_fixed_locs_obs_17["fYear"].unique()],
      names=["fYear"]
    ).to_frame(index=False)
    newdata_17
    ## ----end
    pickle.dump(newdata_17, open(product, "wb"))

def incomplete_spatial_newdata_17b_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_17 = pickle.load(open(upstream["incomplete_spatial_data_prep_17_"], "rb"))
    ## ---- python sampled newdata 17b
    benthos_fixed_locs_obs_17[['Longitude', 'Latitude']] = (
      benthos_fixed_locs_obs_17['geometry']
      .str.strip('c()')
      .str.split(',', expand = True)
      .astype(float)
    )
    newdata_17b = (
      benthos_fixed_locs_obs_17.groupby(["Year", "Reef"], as_index=False, observed=True).agg(
        HCC = ("HCC", "mean"),
        CYC = ("CYC", "mean"),
        DHW = ("DHW", "mean"),
        OTHER = ("OTHER", "mean"),
        Longitude = ("Longitude", "mean"),
        Latitude = ("Latitude", "mean")
      ).reset_index()
    )
    newdata_17b
    ## ----end
    pickle.dump(newdata_17b, open(product, "wb"))



## Fit models ============================================

## Northern sampled reefs --------------------------------------------------

def sampled_reefs_bart_data_12_(product, upstream):
    ## This function prepares the data for the BART model
    ## It prepares two datasets:
    ## 1. `sample_data`: The sampled data (25 sites (etc), 12 years) that
    ##    will be used to fit the model.
    ## 2. `full_data`: The full dataset (all sites, all years) that will
    ##    be used to make predictions on the scale of the full
    ##    spatio-temporal domain.
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_12 = pickle.load(open(upstream["incomplete_spatial_data_prep_12_"], "rb"))
    newdata_10b = pickle.load(open(upstream["incomplete_spatial_newdata_10b_"], "rb"))
    newdata_11b = pickle.load(open(upstream["incomplete_spatial_newdata_11b_"], "rb"))
    newdata_12b = pickle.load(open(upstream["incomplete_spatial_newdata_12b_"], "rb"))
    newdata_13b = pickle.load(open(upstream["incomplete_spatial_newdata_13b_"], "rb"))
    ## ---- python sampled data barts data 12
    upstream = {"incomplete_spatial_data_prep_12_": "../python/cache/incomplete_spatial_data_prep_12_.pkl",
                "incomplete_spatial_global_parameters_": "../python/cache/incomplete_spatial_global_parameters_.pkl",
                "incomplete_spatial_newdata_10b_": "../python/cache/incomplete_spatial_newdata_10b_.pkl",
                "incomplete_spatial_newdata_11b_": "../python/cache/incomplete_spatial_newdata_11b_.pkl",
                "incomplete_spatial_newdata_12b_": "../python/cache/incomplete_spatial_newdata_12b_.pkl",
                "incomplete_spatial_newdata_13b_": "../python/cache/incomplete_spatial_newdata_13b_.pkl"}


    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_12 = pickle.load(open(upstream["incomplete_spatial_data_prep_12_"], "rb"))
    newdata_12b = pickle.load(open(upstream["incomplete_spatial_newdata_12b_"], "rb"))
    # Prepare the training data 
    benthos_fixed_locs_obs_12_prep = benthos_fixed_locs_obs_12.copy()
    for col in ["Reef", "Site", "Transect"]:
        benthos_fixed_locs_obs_12_prep[col] = benthos_fixed_locs_obs_12_prep[col].astype("category").cat.remove_unused_categories()
    
    sample_data = benthos_fixed_locs_obs_12_prep
    # prepare the northern reefs prediction data (newdata_12b)
    newdata_12b = newdata_12b
    newdata_13b = newdata_13b
    ## ----end
    data_dict = {
        "sample_data": sample_data,
        "newdata_10b": newdata_10b,
        "newdata_11b": newdata_11b,
        "newdata_12b": newdata_12b,
        "newdata_13b": newdata_13b
    }
    pickle.dump(data_dict, open(product, "wb"))

def sampled_reefs_bart_data_12_sample_data_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_12 = pickle.load(open(upstream["incomplete_spatial_data_prep_12_"], "rb"))
    ## ---- python sampled data barts data 12 sample data
    upstream = {"incomplete_spatial_data_prep_12_": "../python/cache/incomplete_spatial_data_prep_12_.pkl",
                "incomplete_spatial_global_parameters_": "../python/cache/incomplete_spatial_global_parameters_.pkl"
                }

    benthos_fixed_locs_obs_12 = pickle.load(open(upstream["incomplete_spatial_data_prep_12_"], "rb"))
    benthos_fixed_locs_obs_12
    benthos_fixed_locs_obs_12.groupby("Year").agg(
        Mean=("HCC", "mean"),
        Median=("HCC", "median")
    ).reset_index()
    ## ----end
    result = 1
    pickle.dump(result, open(product, "wb"))

def sampled_reefs_bart_data_12_newdata_12b_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    newdata_12b = pickle.load(open(upstream["incomplete_spatial_newdata_12b_"], "rb"))
    ## ---- python sampled data barts data 12 newdata 12b
    upstream = {"incomplete_spatial_newdata_12b_": "../python/cache/incomplete_spatial_newdata_12b_.pkl",
                "incomplete_spatial_global_parameters_": "../python/cache/incomplete_spatial_global_parameters_.pkl"
                }
    newdata_12b = pickle.load(open(upstream["incomplete_spatial_newdata_12b_"], "rb"))
    newdata_12b
    newdata_12b.groupby("Year").agg(
        Mean=("HCC", "mean"),
        Median=("HCC", "median")
    ).reset_index()
    ## ----end
    result = 1
    pickle.dump(result, open(product, "wb"))

def sampled_reefs_bart_data_12_newdata_13b_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    newdata_13b = pickle.load(open(upstream["incomplete_spatial_newdata_13b_"], "rb"))
    ## ---- python sampled data barts data 12 newdata 13b
    upstream = {"incomplete_spatial_newdata_13b_": "../python/cache/incomplete_spatial_newdata_13b_.pkl",
                "incomplete_spatial_global_parameters_": "../python/cache/incomplete_spatial_global_parameters_.pkl"
                }
    newdata_13b = pickle.load(open(upstream["incomplete_spatial_newdata_13b_"], "rb"))
    newdata_13b
    newdata_13b.groupby("Year").agg(
        Mean=("HCC", "mean"),
        Median=("HCC", "median")
    ).reset_index()
    ## ----end
    result = 1
    pickle.dump(result, open(product, "wb"))

def sampled_reefs_bart_fit_12_(product, upstream):
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_12_"], "rb"))
    sample_data = data_dict["sample_data"]
    newdata_10b = data_dict["newdata_10b"]
    newdata_11b = data_dict["newdata_11b"]
    newdata_12b = data_dict["newdata_12b"]
    newdata_13b = data_dict["newdata_13b"]
    ## ---- python sampled data barts fit 12
    X = sample_data[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = sample_data["cover"]
    with pm.Model() as model_12:
        data_X = pm.Data("data_X", X)
        data_Y = pm.Data("data_Y", Y)
        μ = pmb.BART("μ", data_X, Y, m=100)
        alpha = pm.HalfNormal("alpha", sigma=2)
        beta = pm.Deterministic("beta", alpha * (1 - μ) / μ)
        y = pm.Beta("y", alpha=alpha, beta=beta, observed=data_Y)
    with model_12:
        idata = pm.sample(random_seed=123)
    with model_12:
        ppc1 = pm.sample_posterior_predictive(idata)

    another_X = newdata_12b[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]].copy()
    X = another_X[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = newdata_12b["HCC"]
    with model_12:
            pm.set_data({"data_X": another_X,
                        "data_Y": np.arange(another_X.shape[0])
                        })
            ppc2 = pm.sample_posterior_predictive(idata)

    another_X = newdata_10b[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]].copy()
    X = another_X[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = newdata_10b["HCC"]
    with model_12:
            pm.set_data({"data_X": another_X,
                        "data_Y": np.arange(another_X.shape[0])
                        })
            ppc3 = pm.sample_posterior_predictive(idata)

    another_X = newdata_11b[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]].copy()
    X = another_X[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = newdata_11b["HCC"]
    with model_12:
            pm.set_data({"data_X": another_X,
                        "data_Y": np.arange(another_X.shape[0])
                        })
            ppc4 = pm.sample_posterior_predictive(idata)
    ## ----end
    trace_dict = {
        "idata": idata,
        "ppc1": ppc1,
        "ppc2": ppc2,
        "ppc3": ppc3,
        "ppc4": ppc4
    }
    pickle.dump(trace_dict, open(product, "wb"))

def sampled_reefs_bart_prep_plot_12_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_12_"], "rb"))
    ppc1 = trace_dict["ppc1"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_12_"], "rb"))
    sample_data = data_dict["sample_data"]
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ---- python sampled data pymc_barts prep plot 12 
    posterior_predictive = ppc1.posterior_predictive["y"]  # replace "y_pred" with your variable name
    year_array = sample_data["Year"].values  # must match the shape of predictions

    # Add 'Year' as a coordinate to the xarray
    posterior_predictive = posterior_predictive.assign_coords({"y_dim_0": year_array})
    # Group by Year and compute mean per chain/draw
    avg_preds = posterior_predictive.groupby("y_dim_0").mean("y_dim_0")
    # Rename dimension for clarity
    avg_preds = avg_preds.rename({"y_dim_0": "Year"})
    modelled_trend = az.summary(avg_preds, hdi_prob=0.95)
    modelled_trend["Year"] = avg_preds.coords["Year"].values

    ## Summarise the sampled data
    sample_data_summary = (
        sample_data.groupby("Year").agg(
            Mean=("HCC", "mean"),
            Median=("HCC", "median")
        ).reset_index()
    )

    ## Summarise the full grid
    benthos_reefs_temporal_summary = (
        benthos_reefs_sf.groupby("Year").agg(
            Mean=("HCC", "mean"),
            Median=("HCC", "median"),
            SD=("HCC", "std"),
            Lower=("HCC", lambda x: x.quantile(0.025)),
            Upper=("HCC", lambda x: x.quantile(0.975))
        ).reset_index()
    )
    years = np.unique(year_array)
    cover_mean = avg_preds.mean(dim=["draw", "chain"])
    cover = avg_preds
    ## ----end
    summary_dict = {
        "modelled_trend": modelled_trend,
        "sample_data_summary": sample_data_summary,
        "benthos_reefs_temporal_summary": benthos_reefs_temporal_summary,
        "years": years,
        "cover_mean": cover_mean,
        "cover": cover
    }
    pickle.dump(summary_dict, open(product, "wb"))

def sampled_reefs_bart_plot_12_(product, upstream):
    # paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_12_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    benthos_reefs_temporal_summary = summary_dict["benthos_reefs_temporal_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python sampled data pymc_barts plot 12 
    def multiply_by_100(x, pos):
        return f"{x * 100:.0f}"
    plt.figure(figsize=(8, 6))
    plt.plot(years, cover_mean, "w", lw=3)
    ## Add true mean line
    #plt.plot(
    #    benthos_reefs_temporal_summary["Year"],
    #    benthos_reefs_temporal_summary["Mean"],
    #    label="True mean",
    #    color="blue"
    #)
    ## Add true median line
    #plt.plot(
    #    benthos_reefs_temporal_summary["Year"],
    #    benthos_reefs_temporal_summary["Median"],
    #    label="True median",
    #    color="green"
    #)
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Mean"],
        label="Sample mean",
        color="blue",
        linestyle="--"
    )
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Median"],
        label="Sample median",
        color="green",
        linestyle="--"
    )
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"95% HDI"}, hdi_prob=0.95)
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"50% HDI", "alpha": 1}, hdi_prob=0.50)
    # Add labels and title
    # Apply the custom formatter to the y-axis
    plt.gca().yaxis.set_major_formatter(FuncFormatter(multiply_by_100))
    plt.xlabel("Year")
    plt.ylabel("Mean")
    plt.title("Mean with HDI per Year")
    plt.legend()
    ## ----end
    plt.savefig(f"{paths['fig_path']}python_pdp_1_mod_pymc_bart_12b.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_mse_12b_1_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_12_"], "rb"))
    ppc2 = trace_dict["ppc2"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_12_"], "rb"))
    newdata_12b = data_dict["newdata_12b"]
    ## ---- python sampled data pymc_barts mse 12b 1 
    posterior_predictive = ppc2.posterior_predictive["y"]
    mse = (newdata_12b["HCC"] - posterior_predictive.mean(dim=["draw", "chain"]))**2
    acc = np.abs(np.exp(np.log(newdata_12b["HCC"]) - np.log(posterior_predictive.mean(dim=["draw", "chain"]))) - 1)
    df = pd.DataFrame({
        "mse_mean": [np.mean(mse)],
        "mse_median": [np.median(mse)],
        "acc_mean": [np.mean(acc)],
        "acc_median": [np.median(acc)],
        "model": ["pymc-bart"],
        "lower": [np.quantile(mse, 0.025)],
        "upper": [np.quantile(mse, 0.975)],
        "type": 1,
        "model_type": "covariates"
    })
    df.to_csv(f"{paths['data_path']}/modelled/pymc_bart_12b_mse_1.csv", index=False)
    ## ----end
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_prep_plot_12_2_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_12_"], "rb"))
    ppc3 = trace_dict["ppc3"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_12_"], "rb"))
    ## sample_data = data_dict["sample_data"]
    ## benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    newdata_10b = data_dict["newdata_10b"]
    ## ---- python sampled data pymc_barts prep plot 12 2
    posterior_predictive = ppc3.posterior_predictive["y"]  # replace "y_pred" with your variable name
    year_array = newdata_10b["Year"].values  # must match the shape of predictions

    # Add 'Year' as a coordinate to the xarray
    posterior_predictive = posterior_predictive.assign_coords({"y_dim_0": year_array})
    # Group by Year and compute mean per chain/draw
    avg_preds = posterior_predictive.groupby("y_dim_0").mean("y_dim_0")
    # Rename dimension for clarity
    avg_preds = avg_preds.rename({"y_dim_0": "Year"})
    modelled_trend = az.summary(avg_preds, hdi_prob=0.95)
    modelled_trend["Year"] = avg_preds.coords["Year"].values

    ## Summarise the sampled data
    sample_data_summary = (
        newdata_10b.groupby("Year").agg(
            Mean=("HCC", "mean"),
            Median=("HCC", "median")
        ).reset_index()
    )

    years = np.unique(year_array)
    cover_mean = avg_preds.mean(dim=["draw", "chain"])
    cover = avg_preds
    ## ----end
    summary_dict = {
        "modelled_trend": modelled_trend,
        "sample_data_summary": sample_data_summary,
        "years": years,
        "cover_mean": cover_mean,
        "cover": cover
    }
    pickle.dump(summary_dict, open(product, "wb"))

def sampled_reefs_bart_plot_12_2_(product, upstream):
    # paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_12_2_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python sampled data pymc_barts plot 12 2
    def multiply_by_100(x, pos):
        return f"{x * 100:.0f}"
    plt.figure(figsize=(8, 6))
    plt.plot(years, cover_mean, "w", lw=3)
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Mean"],
        label="Sample mean",
        color="blue",
        linestyle="--"
    )
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Median"],
        label="Sample median",
        color="green",
        linestyle="--"
    )
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"95% HDI"}, hdi_prob=0.95)
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"50% HDI", "alpha": 1}, hdi_prob=0.50)
    # Add labels and title
    # Apply the custom formatter to the y-axis
    plt.gca().yaxis.set_major_formatter(FuncFormatter(multiply_by_100))
    plt.xlabel("Year")
    plt.ylabel("Mean")
    plt.title("Mean with HDI per Year")
    plt.legend()
    ## ----end
    plt.savefig(f"{paths['fig_path']}python_pdp_2_mod_pymc_bart_12b.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_mse_12b_2_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_12_"], "rb"))
    ppc3 = trace_dict["ppc3"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_12_"], "rb"))
    newdata_10b = data_dict["newdata_10b"]
    ## ---- python sampled data pymc_barts mse 12b 2 
    posterior_predictive = ppc3.posterior_predictive["y"]
    mse = (newdata_10b["HCC"] - posterior_predictive.mean(dim=["draw", "chain"]))**2
    acc = np.abs(np.exp(np.log(newdata_10b["HCC"]) - np.log(posterior_predictive.mean(dim=["draw", "chain"]))) - 1)
    df = pd.DataFrame({
        "mse_mean": [np.mean(mse)],
        "mse_median": [np.median(mse)],
        "acc_mean": [np.mean(acc)],
        "acc_median": [np.median(acc)],
        "model": ["pymc-bart"],
        "lower": [np.quantile(mse, 0.025)],
        "upper": [np.quantile(mse, 0.975)],
        "type": 2,
        "model_type": "covariates"
    })
    df.to_csv(f"{paths['data_path']}/modelled/pymc_bart_12b_mse_2.csv", index=False)
    ## ----end
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_prep_plot_12_3_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_12_"], "rb"))
    ppc4 = trace_dict["ppc4"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_12_"], "rb"))
    ## sample_data = data_dict["sample_data"]
    ## benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    newdata_11b = data_dict["newdata_11b"]
    ## ---- python sampled data pymc_barts prep plot 12 3
    posterior_predictive = ppc4.posterior_predictive["y"]  # replace "y_pred" with your variable name
    year_array = newdata_11b["Year"].values  # must match the shape of predictions

    # Add 'Year' as a coordinate to the xarray
    posterior_predictive = posterior_predictive.assign_coords({"y_dim_0": year_array})
    # Group by Year and compute mean per chain/draw
    avg_preds = posterior_predictive.groupby("y_dim_0").mean("y_dim_0")
    # Rename dimension for clarity
    avg_preds = avg_preds.rename({"y_dim_0": "Year"})
    modelled_trend = az.summary(avg_preds, hdi_prob=0.95)
    modelled_trend["Year"] = avg_preds.coords["Year"].values

    ## Summarise the sampled data
    sample_data_summary = (
        newdata_11b.groupby("Year").agg(
            Mean=("HCC", "mean"),
            Median=("HCC", "median")
        ).reset_index()
    )

    years = np.unique(year_array)
    cover_mean = avg_preds.mean(dim=["draw", "chain"])
    cover = avg_preds
    ## ----end
    summary_dict = {
        "modelled_trend": modelled_trend,
        "sample_data_summary": sample_data_summary,
        "years": years,
        "cover_mean": cover_mean,
        "cover": cover
    }
    pickle.dump(summary_dict, open(product, "wb"))

def sampled_reefs_bart_plot_12_3_(product, upstream):
    # paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_12_3_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python sampled data pymc_barts plot 12 3
    def multiply_by_100(x, pos):
        return f"{x * 100:.0f}"
    plt.figure(figsize=(8, 6))
    plt.plot(years, cover_mean, "w", lw=3)
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Mean"],
        label="Sample mean",
        color="blue",
        linestyle="--"
    )
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Median"],
        label="Sample median",
        color="green",
        linestyle="--"
    )
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"95% HDI"}, hdi_prob=0.95)
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"50% HDI", "alpha": 1}, hdi_prob=0.50)
    # Add labels and title
    # Apply the custom formatter to the y-axis
    plt.gca().yaxis.set_major_formatter(FuncFormatter(multiply_by_100))
    plt.xlabel("Year")
    plt.ylabel("Mean")
    plt.title("Mean with HDI per Year")
    plt.legend()
    ## ----end
    plt.savefig(f"{paths['fig_path']}python_pdp_3_mod_pymc_bart_12b.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_mse_12b_3_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_12_"], "rb"))
    ppc4 = trace_dict["ppc4"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_12_"], "rb"))
    newdata_11b = data_dict["newdata_11b"]
    ## ---- python sampled data pymc_barts mse 12b 3 
    posterior_predictive = ppc4.posterior_predictive["y"]
    mse = (newdata_11b["HCC"] - posterior_predictive.mean(dim=["draw", "chain"]))**2
    acc = np.abs(np.exp(np.log(newdata_11b["HCC"]) - np.log(posterior_predictive.mean(dim=["draw", "chain"]))) - 1)
    df = pd.DataFrame({
        "mse_mean": [np.mean(mse)],
        "mse_median": [np.median(mse)],
        "acc_mean": [np.mean(acc)],
        "acc_median": [np.median(acc)],
        "model": ["pymc-bart"],
        "lower": [np.quantile(mse, 0.025)],
        "upper": [np.quantile(mse, 0.975)],
        "type": 3,
        "model_type": "covariates"
    })
    df.to_csv(f"{paths['data_path']}/modelled/pymc_bart_12b_mse_3.csv", index=False)
    ## ----end
    results=1
    pickle.dump(results, open(product, "wb"))


## Southern sampled reefs --------------------------------------------------

def sampled_reefs_bart_data_13_(product, upstream):
    ## This function prepares the data for the BART model
    ## It prepares two datasets:
    ## 1. `sample_data`: The sampled data (25 sites (etc), 12 years) that
    ##    will be used to fit the model.
    ## 2. `full_data`: The full dataset (all sites, all years) that will
    ##    be used to make predictions on the scale of the full
    ##    spatio-temporal domain.
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_13 = pickle.load(open(upstream["incomplete_spatial_data_prep_13_"], "rb"))
    newdata_10b = pickle.load(open(upstream["incomplete_spatial_newdata_10b_"], "rb"))
    newdata_11b = pickle.load(open(upstream["incomplete_spatial_newdata_11b_"], "rb"))
    newdata_12b = pickle.load(open(upstream["incomplete_spatial_newdata_12b_"], "rb"))
    newdata_13b = pickle.load(open(upstream["incomplete_spatial_newdata_13b_"], "rb"))
    ## ---- python sampled data barts data 13
    upstream = {"incomplete_spatial_data_prep_13_": "../python/cache/incomplete_spatial_data_prep_13_.pkl",
                "incomplete_spatial_global_parameters_": "../python/cache/incomplete_spatial_global_parameters_.pkl",
                "incomplete_spatial_newdata_10b_": "../python/cache/incomplete_spatial_newdata_10b_.pkl",
                "incomplete_spatial_newdata_11b_": "../python/cache/incomplete_spatial_newdata_11b_.pkl",
                "incomplete_spatial_newdata_12b_": "../python/cache/incomplete_spatial_newdata_12b_.pkl",
                "incomplete_spatial_newdata_13b_": "../python/cache/incomplete_spatial_newdata_13b_.pkl"}


    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_13 = pickle.load(open(upstream["incomplete_spatial_data_prep_13_"], "rb"))
    newdata_13b = pickle.load(open(upstream["incomplete_spatial_newdata_13b_"], "rb"))
    # Prepare the training data 
    benthos_fixed_locs_obs_13_prep = benthos_fixed_locs_obs_13.copy()
    for col in ["Reef", "Site", "Transect"]:
        benthos_fixed_locs_obs_13_prep[col] = benthos_fixed_locs_obs_13_prep[col].astype("category").cat.remove_unused_categories()
    
    sample_data = benthos_fixed_locs_obs_13_prep
    # prepare the northern reefs prediction data (newdata_13b)
    newdata_12b = newdata_12b
    newdata_13b = newdata_13b
    ## ----end
    data_dict = {
        "sample_data": sample_data,
        "newdata_10b": newdata_10b,
        "newdata_11b": newdata_11b,
        "newdata_12b": newdata_12b,
        "newdata_13b": newdata_13b
    }
    pickle.dump(data_dict, open(product, "wb"))

def sampled_reefs_bart_data_13_sample_data_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_13 = pickle.load(open(upstream["incomplete_spatial_data_prep_13_"], "rb"))
    ## ---- python sampled data barts data 13 sample data
    upstream = {"incomplete_spatial_data_prep_13_": "../python/cache/incomplete_spatial_data_prep_13_.pkl",
                "incomplete_spatial_global_parameters_": "../python/cache/incomplete_spatial_global_parameters_.pkl"
                }

    benthos_fixed_locs_obs_13 = pickle.load(open(upstream["incomplete_spatial_data_prep_13_"], "rb"))
    benthos_fixed_locs_obs_13
    benthos_fixed_locs_obs_13.groupby("Year").agg(
        Mean=("HCC", "mean"),
        Median=("HCC", "median")
    ).reset_index()
    ## ----end
    result = 1
    pickle.dump(result, open(product, "wb"))

def sampled_reefs_bart_fit_13_(product, upstream):
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_13_"], "rb"))
    sample_data = data_dict["sample_data"]
    newdata_10b = data_dict["newdata_10b"]
    newdata_11b = data_dict["newdata_11b"]
    newdata_12b = data_dict["newdata_12b"]
    newdata_13b = data_dict["newdata_13b"]
    ## ---- python sampled data barts fit 13
    X = sample_data[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = sample_data["cover"]
    with pm.Model() as model_13:
        data_X = pm.Data("data_X", X)
        data_Y = pm.Data("data_Y", Y)
        μ = pmb.BART("μ", data_X, Y, m=100)
        alpha = pm.HalfNormal("alpha", sigma=2)
        beta = pm.Deterministic("beta", alpha * (1 - μ) / μ)
        y = pm.Beta("y", alpha=alpha, beta=beta, observed=data_Y)
    with model_13:
        idata = pm.sample(random_seed=123)
    with model_13:
        ppc1 = pm.sample_posterior_predictive(idata)

    another_X = newdata_13b[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]].copy()
    X = another_X[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = newdata_13b["HCC"]
    with model_13:
            pm.set_data({"data_X": another_X,
                        "data_Y": np.arange(another_X.shape[0])
                        })
            ppc2 = pm.sample_posterior_predictive(idata)

    another_X = newdata_11b[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]].copy()
    X = another_X[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = newdata_11b["HCC"]
    with model_13:
            pm.set_data({"data_X": another_X,
                        "data_Y": np.arange(another_X.shape[0])
                        })
            ppc3 = pm.sample_posterior_predictive(idata)

    another_X = newdata_10b[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]].copy()
    X = another_X[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = newdata_10b["HCC"]
    with model_13:
            pm.set_data({"data_X": another_X,
                        "data_Y": np.arange(another_X.shape[0])
                        })
            ppc4 = pm.sample_posterior_predictive(idata)
    ## ----end
    trace_dict = {
        "idata": idata,
        "ppc1": ppc1,
        "ppc2": ppc2,
        "ppc3": ppc3,
        "ppc4": ppc4
    }
    pickle.dump(trace_dict, open(product, "wb"))

def sampled_reefs_bart_prep_plot_13_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_13_"], "rb"))
    ppc1 = trace_dict["ppc1"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_13_"], "rb"))
    sample_data = data_dict["sample_data"]
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ---- python sampled data pymc_barts prep plot 13 
    posterior_predictive = ppc1.posterior_predictive["y"]  # replace "y_pred" with your variable name
    year_array = sample_data["Year"].values  # must match the shape of predictions

    # Add 'Year' as a coordinate to the xarray
    posterior_predictive = posterior_predictive.assign_coords({"y_dim_0": year_array})
    # Group by Year and compute mean per chain/draw
    avg_preds = posterior_predictive.groupby("y_dim_0").mean("y_dim_0")
    # Rename dimension for clarity
    avg_preds = avg_preds.rename({"y_dim_0": "Year"})
    modelled_trend = az.summary(avg_preds, hdi_prob=0.95)
    modelled_trend["Year"] = avg_preds.coords["Year"].values

    ## Summarise the sampled data
    sample_data_summary = (
        sample_data.groupby("Year").agg(
            Mean=("HCC", "mean"),
            Median=("HCC", "median")
        ).reset_index()
    )

    ## Summarise the full grid
    benthos_reefs_temporal_summary = (
        benthos_reefs_sf.groupby("Year").agg(
            Mean=("HCC", "mean"),
            Median=("HCC", "median"),
            SD=("HCC", "std"),
            Lower=("HCC", lambda x: x.quantile(0.025)),
            Upper=("HCC", lambda x: x.quantile(0.975))
        ).reset_index()
    )
    years = np.unique(year_array)
    cover_mean = avg_preds.mean(dim=["draw", "chain"])
    cover = avg_preds
    ## ----end
    summary_dict = {
        "modelled_trend": modelled_trend,
        "sample_data_summary": sample_data_summary,
        "benthos_reefs_temporal_summary": benthos_reefs_temporal_summary,
        "years": years,
        "cover_mean": cover_mean,
        "cover": cover
    }
    pickle.dump(summary_dict, open(product, "wb"))

def sampled_reefs_bart_plot_13_(product, upstream):
    # paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_13_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    benthos_reefs_temporal_summary = summary_dict["benthos_reefs_temporal_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python sampled data pymc_barts plot 13 
    def multiply_by_100(x, pos):
        return f"{x * 100:.0f}"
    plt.figure(figsize=(8, 6))
    plt.plot(years, cover_mean, "w", lw=3)
    ## Add true mean line
    #plt.plot(
    #    benthos_reefs_temporal_summary["Year"],
    #    benthos_reefs_temporal_summary["Mean"],
    #    label="True mean",
    #    color="blue"
    #)
    ## Add true median line
    #plt.plot(
    #    benthos_reefs_temporal_summary["Year"],
    #    benthos_reefs_temporal_summary["Median"],
    #    label="True median",
    #    color="green"
    #)
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Mean"],
        label="Sample mean",
        color="blue",
        linestyle="--"
    )
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Median"],
        label="Sample median",
        color="green",
        linestyle="--"
    )
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"95% HDI"}, hdi_prob=0.95)
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"50% HDI", "alpha": 1}, hdi_prob=0.50)
    # Add labels and title
    # Apply the custom formatter to the y-axis
    plt.gca().yaxis.set_major_formatter(FuncFormatter(multiply_by_100))
    plt.xlabel("Year")
    plt.ylabel("Mean")
    plt.title("Mean with HDI per Year")
    plt.legend()
    ## ----end
    plt.savefig(f"{paths['fig_path']}python_pdp_1_mod_pymc_bart_13b.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_mse_13b_1_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_13_"], "rb"))
    ppc2 = trace_dict["ppc2"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_13_"], "rb"))
    newdata_13b = data_dict["newdata_13b"]
    ## ---- python sampled data pymc_barts mse 13b 1 
    posterior_predictive = ppc2.posterior_predictive["y"]
    mse = (newdata_13b["HCC"] - posterior_predictive.mean(dim=["draw", "chain"]))**2
    acc = np.abs(np.exp(np.log(newdata_13b["HCC"]) - np.log(posterior_predictive.mean(dim=["draw", "chain"]))) - 1)
    df = pd.DataFrame({
        "mse_mean": [np.mean(mse)],
        "mse_median": [np.median(mse)],
        "acc_mean": [np.mean(acc)],
        "acc_median": [np.median(acc)],
        "model": ["pymc-bart"],
        "lower": [np.quantile(mse, 0.025)],
        "upper": [np.quantile(mse, 0.975)],
        "type": 1,
        "model_type": "covariates"
    })
    df.to_csv(f"{paths['data_path']}/modelled/pymc_bart_13b_mse_1.csv", index=False)
    ## ----end
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_prep_plot_13_2_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_13_"], "rb"))
    ppc3 = trace_dict["ppc3"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_13_"], "rb"))
    ## sample_data = data_dict["sample_data"]
    ## benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    newdata_11b = data_dict["newdata_11b"]
    ## ---- python sampled data pymc_barts prep plot 13 2
    posterior_predictive = ppc3.posterior_predictive["y"]  # replace "y_pred" with your variable name
    year_array = newdata_11b["Year"].values  # must match the shape of predictions

    # Add 'Year' as a coordinate to the xarray
    posterior_predictive = posterior_predictive.assign_coords({"y_dim_0": year_array})
    # Group by Year and compute mean per chain/draw
    avg_preds = posterior_predictive.groupby("y_dim_0").mean("y_dim_0")
    # Rename dimension for clarity
    avg_preds = avg_preds.rename({"y_dim_0": "Year"})
    modelled_trend = az.summary(avg_preds, hdi_prob=0.95)
    modelled_trend["Year"] = avg_preds.coords["Year"].values

    ## Summarise the sampled data
    sample_data_summary = (
        newdata_11b.groupby("Year").agg(
            Mean=("HCC", "mean"),
            Median=("HCC", "median")
        ).reset_index()
    )

    years = np.unique(year_array)
    cover_mean = avg_preds.mean(dim=["draw", "chain"])
    cover = avg_preds
    ## ----end
    summary_dict = {
        "modelled_trend": modelled_trend,
        "sample_data_summary": sample_data_summary,
        "years": years,
        "cover_mean": cover_mean,
        "cover": cover
    }
    pickle.dump(summary_dict, open(product, "wb"))

def sampled_reefs_bart_plot_13_2_(product, upstream):
    # paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_13_2_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python sampled data pymc_barts plot 13 2
    def multiply_by_100(x, pos):
        return f"{x * 100:.0f}"
    plt.figure(figsize=(8, 6))
    plt.plot(years, cover_mean, "w", lw=3)
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Mean"],
        label="Sample mean",
        color="blue",
        linestyle="--"
    )
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Median"],
        label="Sample median",
        color="green",
        linestyle="--"
    )
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"95% HDI"}, hdi_prob=0.95)
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"50% HDI", "alpha": 1}, hdi_prob=0.50)
    # Add labels and title
    # Apply the custom formatter to the y-axis
    plt.gca().yaxis.set_major_formatter(FuncFormatter(multiply_by_100))
    plt.xlabel("Year")
    plt.ylabel("Mean")
    plt.title("Mean with HDI per Year")
    plt.legend()
    ## ----end
    plt.savefig(f"{paths['fig_path']}python_pdp_2_mod_pymc_bart_13b.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_mse_13b_2_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_13_"], "rb"))
    ppc3 = trace_dict["ppc3"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_13_"], "rb"))
    newdata_11b = data_dict["newdata_11b"]
    ## ---- python sampled data pymc_barts mse 13b 2 
    posterior_predictive = ppc3.posterior_predictive["y"]
    mse = (newdata_11b["HCC"] - posterior_predictive.mean(dim=["draw", "chain"]))**2
    acc = np.abs(np.exp(np.log(newdata_11b["HCC"]) - np.log(posterior_predictive.mean(dim=["draw", "chain"]))) - 1)
    df = pd.DataFrame({
        "mse_mean": [np.mean(mse)],
        "mse_median": [np.median(mse)],
        "acc_mean": [np.mean(acc)],
        "acc_median": [np.median(acc)],
        "model": ["pymc-bart"],
        "lower": [np.quantile(mse, 0.025)],
        "upper": [np.quantile(mse, 0.975)],
        "type": 2,
        "model_type": "covariates"
    })
    df.to_csv(f"{paths['data_path']}/modelled/pymc_bart_13b_mse_2.csv", index=False)
    ## ----end
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_prep_plot_13_3_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_13_"], "rb"))
    ppc4 = trace_dict["ppc4"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_13_"], "rb"))
    ## sample_data = data_dict["sample_data"]
    ## benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    newdata_10b = data_dict["newdata_10b"]
    ## ---- python sampled data pymc_barts prep plot 12 3
    posterior_predictive = ppc4.posterior_predictive["y"]  # replace "y_pred" with your variable name
    year_array = newdata_10b["Year"].values  # must match the shape of predictions

    # Add 'Year' as a coordinate to the xarray
    posterior_predictive = posterior_predictive.assign_coords({"y_dim_0": year_array})
    # Group by Year and compute mean per chain/draw
    avg_preds = posterior_predictive.groupby("y_dim_0").mean("y_dim_0")
    # Rename dimension for clarity
    avg_preds = avg_preds.rename({"y_dim_0": "Year"})
    modelled_trend = az.summary(avg_preds, hdi_prob=0.95)
    modelled_trend["Year"] = avg_preds.coords["Year"].values

    ## Summarise the sampled data
    sample_data_summary = (
        newdata_10b.groupby("Year").agg(
            Mean=("HCC", "mean"),
            Median=("HCC", "median")
        ).reset_index()
    )

    years = np.unique(year_array)
    cover_mean = avg_preds.mean(dim=["draw", "chain"])
    cover = avg_preds
    ## ----end
    summary_dict = {
        "modelled_trend": modelled_trend,
        "sample_data_summary": sample_data_summary,
        "years": years,
        "cover_mean": cover_mean,
        "cover": cover
    }
    pickle.dump(summary_dict, open(product, "wb"))

def sampled_reefs_bart_plot_13_3_(product, upstream):
    # paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_13_3_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python sampled data pymc_barts plot 13 3
    def multiply_by_100(x, pos):
        return f"{x * 100:.0f}"
    plt.figure(figsize=(8, 6))
    plt.plot(years, cover_mean, "w", lw=3)
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Mean"],
        label="Sample mean",
        color="blue",
        linestyle="--"
    )
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Median"],
        label="Sample median",
        color="green",
        linestyle="--"
    )
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"95% HDI"}, hdi_prob=0.95)
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"50% HDI", "alpha": 1}, hdi_prob=0.50)
    # Add labels and title
    # Apply the custom formatter to the y-axis
    plt.gca().yaxis.set_major_formatter(FuncFormatter(multiply_by_100))
    plt.xlabel("Year")
    plt.ylabel("Mean")
    plt.title("Mean with HDI per Year")
    plt.legend()
    ## ----end
    plt.savefig(f"{paths['fig_path']}python_pdp_3_mod_pymc_bart_13b.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_mse_13b_3_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_13_"], "rb"))
    ppc4 = trace_dict["ppc4"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_13_"], "rb"))
    newdata_10b = data_dict["newdata_10b"]
    ## ---- python sampled data pymc_barts mse 12b 3 
    posterior_predictive = ppc4.posterior_predictive["y"]
    mse = (newdata_10b["HCC"] - posterior_predictive.mean(dim=["draw", "chain"]))**2
    acc = np.abs(np.exp(np.log(newdata_10b["HCC"]) - np.log(posterior_predictive.mean(dim=["draw", "chain"]))) - 1)
    df = pd.DataFrame({
        "mse_mean": [np.mean(mse)],
        "mse_median": [np.median(mse)],
        "acc_mean": [np.mean(acc)],
        "acc_median": [np.median(acc)],
        "model": ["pymc-bart"],
        "lower": [np.quantile(mse, 0.025)],
        "upper": [np.quantile(mse, 0.975)],
        "type": 3,
        "model_type": "covariates"
    })
    df.to_csv(f"{paths['data_path']}/modelled/pymc_bart_13b_mse_3.csv", index=False)
    ## ----end
    results=1
    pickle.dump(results, open(product, "wb"))

## Western sampled reefs --------------------------------------------------

def sampled_reefs_bart_data_16_(product, upstream):
    ## This function prepares the data for the BART model
    ## It prepares two datasets:
    ## 1. `sample_data`: The sampled data (25 sites (etc), 16 years) that
    ##    will be used to fit the model.
    ## 2. `full_data`: The full dataset (all sites, all years) that will
    ##    be used to make predictions on the scale of the full
    ##    spatio-temporal domain.
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_16 = pickle.load(open(upstream["incomplete_spatial_data_prep_16_"], "rb"))
    newdata_14b = pickle.load(open(upstream["incomplete_spatial_newdata_14b_"], "rb"))
    newdata_15b = pickle.load(open(upstream["incomplete_spatial_newdata_15b_"], "rb"))
    newdata_16b = pickle.load(open(upstream["incomplete_spatial_newdata_16b_"], "rb"))
    newdata_17b = pickle.load(open(upstream["incomplete_spatial_newdata_17b_"], "rb"))
    ## ---- python sampled data barts data 16
    upstream = {"incomplete_spatial_data_prep_16_": "../python/cache/incomplete_spatial_data_prep_16_.pkl",
                "incomplete_spatial_global_parameters_": "../python/cache/incomplete_spatial_global_parameters_.pkl",
                "incomplete_spatial_newdata_14b_": "../python/cache/incomplete_spatial_newdata_14b_.pkl",
                "incomplete_spatial_newdata_15b_": "../python/cache/incomplete_spatial_newdata_15b_.pkl",
                "incomplete_spatial_newdata_16b_": "../python/cache/incomplete_spatial_newdata_16b_.pkl",
                "incomplete_spatial_newdata_17b_": "../python/cache/incomplete_spatial_newdata_17b_.pkl"}


    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_16 = pickle.load(open(upstream["incomplete_spatial_data_prep_16_"], "rb"))
    newdata_16b = pickle.load(open(upstream["incomplete_spatial_newdata_16b_"], "rb"))
    # Prepare the training data 
    benthos_fixed_locs_obs_16_prep = benthos_fixed_locs_obs_16.copy()
    for col in ["Reef", "Site", "Transect"]:
        benthos_fixed_locs_obs_16_prep[col] = benthos_fixed_locs_obs_16_prep[col].astype("category").cat.remove_unused_categories()
    
    sample_data = benthos_fixed_locs_obs_16_prep
    # prepare the northern reefs prediction data (newdata_16b)
    newdata_16b = newdata_16b
    newdata_17b = newdata_17b
    ## ----end
    data_dict = {
        "sample_data": sample_data,
        "newdata_14b": newdata_14b,
        "newdata_15b": newdata_15b,
        "newdata_16b": newdata_16b,
        "newdata_17b": newdata_17b
    }
    pickle.dump(data_dict, open(product, "wb"))

def sampled_reefs_bart_data_16_sample_data_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_16 = pickle.load(open(upstream["incomplete_spatial_data_prep_16_"], "rb"))
    ## ---- python sampled data barts data 16 sample data
    upstream = {"incomplete_spatial_data_prep_16_": "../python/cache/incomplete_spatial_data_prep_16_.pkl",
                "incomplete_spatial_global_parameters_": "../python/cache/incomplete_spatial_global_parameters_.pkl"
                }

    benthos_fixed_locs_obs_16 = pickle.load(open(upstream["incomplete_spatial_data_prep_16_"], "rb"))
    benthos_fixed_locs_obs_16
    benthos_fixed_locs_obs_16.groupby("Year").agg(
        Mean=("HCC", "mean"),
        Median=("HCC", "median")
    ).reset_index()
    ## ----end
    result = 1
    pickle.dump(result, open(product, "wb"))

def sampled_reefs_bart_data_16_newdata_16b_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    newdata_16b = pickle.load(open(upstream["incomplete_spatial_newdata_16b_"], "rb"))
    ## ---- python sampled data barts data 16 newdata 16b
    upstream = {"incomplete_spatial_newdata_16b_": "../python/cache/incomplete_spatial_newdata_16b_.pkl",
                "incomplete_spatial_global_parameters_": "../python/cache/incomplete_spatial_global_parameters_.pkl"
                }
    newdata_16b = pickle.load(open(upstream["incomplete_spatial_newdata_16b_"], "rb"))
    newdata_16b
    newdata_16b.groupby("Year").agg(
        Mean=("HCC", "mean"),
        Median=("HCC", "median")
    ).reset_index()
    ## ----end
    result = 1
    pickle.dump(result, open(product, "wb"))

def sampled_reefs_bart_data_16_newdata_17b_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    newdata_17b = pickle.load(open(upstream["incomplete_spatial_newdata_17b_"], "rb"))
    ## ---- python sampled data barts data 16 newdata 17b
    upstream = {"incomplete_spatial_newdata_17b_": "../python/cache/incomplete_spatial_newdata_17b_.pkl",
                "incomplete_spatial_global_parameters_": "../python/cache/incomplete_spatial_global_parameters_.pkl"
                }
    newdata_17b = pickle.load(open(upstream["incomplete_spatial_newdata_17b_"], "rb"))
    newdata_17b
    newdata_17b.groupby("Year").agg(
        Mean=("HCC", "mean"),
        Median=("HCC", "median")
    ).reset_index()
    ## ----end
    result = 1
    pickle.dump(result, open(product, "wb"))

def sampled_reefs_bart_fit_16_(product, upstream):
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_16_"], "rb"))
    sample_data = data_dict["sample_data"]
    newdata_14b = data_dict["newdata_14b"]
    newdata_15b = data_dict["newdata_15b"]
    newdata_16b = data_dict["newdata_16b"]
    newdata_17b = data_dict["newdata_17b"]
    ## ---- python sampled data barts fit 16
    X = sample_data[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = sample_data["cover"]
    with pm.Model() as model_16:
        data_X = pm.Data("data_X", X)
        data_Y = pm.Data("data_Y", Y)
        μ = pmb.BART("μ", data_X, Y, m=100)
        alpha = pm.HalfNormal("alpha", sigma=2)
        beta = pm.Deterministic("beta", alpha * (1 - μ) / μ)
        y = pm.Beta("y", alpha=alpha, beta=beta, observed=data_Y)
    with model_16:
        idata = pm.sample(random_seed=163)
    with model_16:
        ppc1 = pm.sample_posterior_predictive(idata)

    another_X = newdata_16b[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]].copy()
    X = another_X[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = newdata_16b["HCC"]
    with model_16:
            pm.set_data({"data_X": another_X,
                        "data_Y": np.arange(another_X.shape[0])
                        })
            ppc2 = pm.sample_posterior_predictive(idata)

    another_X = newdata_14b[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]].copy()
    X = another_X[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = newdata_14b["HCC"]
    with model_16:
            pm.set_data({"data_X": another_X,
                        "data_Y": np.arange(another_X.shape[0])
                        })
            ppc3 = pm.sample_posterior_predictive(idata)

    another_X = newdata_15b[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]].copy()
    X = another_X[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = newdata_15b["HCC"]
    with model_16:
            pm.set_data({"data_X": another_X,
                        "data_Y": np.arange(another_X.shape[0])
                        })
            ppc4 = pm.sample_posterior_predictive(idata)
    ## ----end
    trace_dict = {
        "idata": idata,
        "ppc1": ppc1,
        "ppc2": ppc2,
        "ppc3": ppc3,
        "ppc4": ppc4
    }
    pickle.dump(trace_dict, open(product, "wb"))

def sampled_reefs_bart_prep_plot_16_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_16_"], "rb"))
    ppc1 = trace_dict["ppc1"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_16_"], "rb"))
    sample_data = data_dict["sample_data"]
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ---- python sampled data pymc_barts prep plot 16 
    posterior_predictive = ppc1.posterior_predictive["y"]  # replace "y_pred" with your variable name
    year_array = sample_data["Year"].values  # must match the shape of predictions

    # Add 'Year' as a coordinate to the xarray
    posterior_predictive = posterior_predictive.assign_coords({"y_dim_0": year_array})
    # Group by Year and compute mean per chain/draw
    avg_preds = posterior_predictive.groupby("y_dim_0").mean("y_dim_0")
    # Rename dimension for clarity
    avg_preds = avg_preds.rename({"y_dim_0": "Year"})
    modelled_trend = az.summary(avg_preds, hdi_prob=0.95)
    modelled_trend["Year"] = avg_preds.coords["Year"].values

    ## Summarise the sampled data
    sample_data_summary = (
        sample_data.groupby("Year").agg(
            Mean=("HCC", "mean"),
            Median=("HCC", "median")
        ).reset_index()
    )

    ## Summarise the full grid
    benthos_reefs_temporal_summary = (
        benthos_reefs_sf.groupby("Year").agg(
            Mean=("HCC", "mean"),
            Median=("HCC", "median"),
            SD=("HCC", "std"),
            Lower=("HCC", lambda x: x.quantile(0.025)),
            Upper=("HCC", lambda x: x.quantile(0.975))
        ).reset_index()
    )
    years = np.unique(year_array)
    cover_mean = avg_preds.mean(dim=["draw", "chain"])
    cover = avg_preds
    ## ----end
    summary_dict = {
        "modelled_trend": modelled_trend,
        "sample_data_summary": sample_data_summary,
        "benthos_reefs_temporal_summary": benthos_reefs_temporal_summary,
        "years": years,
        "cover_mean": cover_mean,
        "cover": cover
    }
    pickle.dump(summary_dict, open(product, "wb"))

def sampled_reefs_bart_plot_16_(product, upstream):
    # paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_16_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    benthos_reefs_temporal_summary = summary_dict["benthos_reefs_temporal_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python sampled data pymc_barts plot 16 
    def multiply_by_100(x, pos):
        return f"{x * 100:.0f}"
    plt.figure(figsize=(8, 6))
    plt.plot(years, cover_mean, "w", lw=3)
    ## Add true mean line
    #plt.plot(
    #    benthos_reefs_temporal_summary["Year"],
    #    benthos_reefs_temporal_summary["Mean"],
    #    label="True mean",
    #    color="blue"
    #)
    ## Add true median line
    #plt.plot(
    #    benthos_reefs_temporal_summary["Year"],
    #    benthos_reefs_temporal_summary["Median"],
    #    label="True median",
    #    color="green"
    #)
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Mean"],
        label="Sample mean",
        color="blue",
        linestyle="--"
    )
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Median"],
        label="Sample median",
        color="green",
        linestyle="--"
    )
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"95% HDI"}, hdi_prob=0.95)
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"50% HDI", "alpha": 1}, hdi_prob=0.50)
    # Add labels and title
    # Apply the custom formatter to the y-axis
    plt.gca().yaxis.set_major_formatter(FuncFormatter(multiply_by_100))
    plt.xlabel("Year")
    plt.ylabel("Mean")
    plt.title("Mean with HDI per Year")
    plt.legend()
    ## ----end
    plt.savefig(f"{paths['fig_path']}python_pdp_1_mod_pymc_bart_16b.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_mse_16b_1_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_16_"], "rb"))
    ppc2 = trace_dict["ppc2"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_16_"], "rb"))
    newdata_16b = data_dict["newdata_16b"]
    ## ---- python sampled data pymc_barts mse 16b 1 
    posterior_predictive = ppc2.posterior_predictive["y"]
    mse = (newdata_16b["HCC"] - posterior_predictive.mean(dim=["draw", "chain"]))**2
    acc = np.abs(np.exp(np.log(newdata_16b["HCC"]) - np.log(posterior_predictive.mean(dim=["draw", "chain"]))) - 1)
    df = pd.DataFrame({
        "mse_mean": [np.mean(mse)],
        "mse_median": [np.median(mse)],
        "acc_mean": [np.mean(acc)],
        "acc_median": [np.median(acc)],
        "model": ["pymc-bart"],
        "lower": [np.quantile(mse, 0.025)],
        "upper": [np.quantile(mse, 0.975)],
        "type": 1,
        "model_type": "covariates"
    })
    df.to_csv(f"{paths['data_path']}/modelled/pymc_bart_16b_mse_1.csv", index=False)
    ## ----end
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_prep_plot_16_2_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_16_"], "rb"))
    ppc3 = trace_dict["ppc3"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_16_"], "rb"))
    ## sample_data = data_dict["sample_data"]
    ## benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    newdata_14b = data_dict["newdata_14b"]
    ## ---- python sampled data pymc_barts prep plot 16 2
    posterior_predictive = ppc3.posterior_predictive["y"]  # replace "y_pred" with your variable name
    year_array = newdata_14b["Year"].values  # must match the shape of predictions

    # Add 'Year' as a coordinate to the xarray
    posterior_predictive = posterior_predictive.assign_coords({"y_dim_0": year_array})
    # Group by Year and compute mean per chain/draw
    avg_preds = posterior_predictive.groupby("y_dim_0").mean("y_dim_0")
    # Rename dimension for clarity
    avg_preds = avg_preds.rename({"y_dim_0": "Year"})
    modelled_trend = az.summary(avg_preds, hdi_prob=0.95)
    modelled_trend["Year"] = avg_preds.coords["Year"].values

    ## Summarise the sampled data
    sample_data_summary = (
        newdata_14b.groupby("Year").agg(
            Mean=("HCC", "mean"),
            Median=("HCC", "median")
        ).reset_index()
    )

    years = np.unique(year_array)
    cover_mean = avg_preds.mean(dim=["draw", "chain"])
    cover = avg_preds
    ## ----end
    summary_dict = {
        "modelled_trend": modelled_trend,
        "sample_data_summary": sample_data_summary,
        "years": years,
        "cover_mean": cover_mean,
        "cover": cover
    }
    pickle.dump(summary_dict, open(product, "wb"))

def sampled_reefs_bart_plot_16_2_(product, upstream):
    # paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_16_2_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python sampled data pymc_barts plot 16 2
    def multiply_by_100(x, pos):
        return f"{x * 100:.0f}"
    plt.figure(figsize=(8, 6))
    plt.plot(years, cover_mean, "w", lw=3)
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Mean"],
        label="Sample mean",
        color="blue",
        linestyle="--"
    )
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Median"],
        label="Sample median",
        color="green",
        linestyle="--"
    )
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"95% HDI"}, hdi_prob=0.95)
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"50% HDI", "alpha": 1}, hdi_prob=0.50)
    # Add labels and title
    # Apply the custom formatter to the y-axis
    plt.gca().yaxis.set_major_formatter(FuncFormatter(multiply_by_100))
    plt.xlabel("Year")
    plt.ylabel("Mean")
    plt.title("Mean with HDI per Year")
    plt.legend()
    ## ----end
    plt.savefig(f"{paths['fig_path']}python_pdp_2_mod_pymc_bart_16b.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_mse_16b_2_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_16_"], "rb"))
    ppc3 = trace_dict["ppc3"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_16_"], "rb"))
    newdata_14b = data_dict["newdata_14b"]
    ## ---- python sampled data pymc_barts mse 16b 2 
    posterior_predictive = ppc3.posterior_predictive["y"]
    mse = (newdata_14b["HCC"] - posterior_predictive.mean(dim=["draw", "chain"]))**2
    acc = np.abs(np.exp(np.log(newdata_14b["HCC"]) - np.log(posterior_predictive.mean(dim=["draw", "chain"]))) - 1)
    df = pd.DataFrame({
        "mse_mean": [np.mean(mse)],
        "mse_median": [np.median(mse)],
        "acc_mean": [np.mean(acc)],
        "acc_median": [np.median(acc)],
        "model": ["pymc-bart"],
        "lower": [np.quantile(mse, 0.025)],
        "upper": [np.quantile(mse, 0.975)],
        "type": 2,
        "model_type": "covariates"
    })
    df.to_csv(f"{paths['data_path']}/modelled/pymc_bart_16b_mse_2.csv", index=False)
    ## ----end
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_prep_plot_16_3_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_16_"], "rb"))
    ppc4 = trace_dict["ppc4"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_16_"], "rb"))
    ## sample_data = data_dict["sample_data"]
    ## benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    newdata_15b = data_dict["newdata_15b"]
    ## ---- python sampled data pymc_barts prep plot 16 3
    posterior_predictive = ppc4.posterior_predictive["y"]  # replace "y_pred" with your variable name
    year_array = newdata_15b["Year"].values  # must match the shape of predictions

    # Add 'Year' as a coordinate to the xarray
    posterior_predictive = posterior_predictive.assign_coords({"y_dim_0": year_array})
    # Group by Year and compute mean per chain/draw
    avg_preds = posterior_predictive.groupby("y_dim_0").mean("y_dim_0")
    # Rename dimension for clarity
    avg_preds = avg_preds.rename({"y_dim_0": "Year"})
    modelled_trend = az.summary(avg_preds, hdi_prob=0.95)
    modelled_trend["Year"] = avg_preds.coords["Year"].values

    ## Summarise the sampled data
    sample_data_summary = (
        newdata_15b.groupby("Year").agg(
            Mean=("HCC", "mean"),
            Median=("HCC", "median")
        ).reset_index()
    )

    years = np.unique(year_array)
    cover_mean = avg_preds.mean(dim=["draw", "chain"])
    cover = avg_preds
    ## ----end
    summary_dict = {
        "modelled_trend": modelled_trend,
        "sample_data_summary": sample_data_summary,
        "years": years,
        "cover_mean": cover_mean,
        "cover": cover
    }
    pickle.dump(summary_dict, open(product, "wb"))

def sampled_reefs_bart_plot_16_3_(product, upstream):
    # paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_16_3_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python sampled data pymc_barts plot 16 3
    def multiply_by_100(x, pos):
        return f"{x * 100:.0f}"
    plt.figure(figsize=(8, 6))
    plt.plot(years, cover_mean, "w", lw=3)
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Mean"],
        label="Sample mean",
        color="blue",
        linestyle="--"
    )
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Median"],
        label="Sample median",
        color="green",
        linestyle="--"
    )
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"95% HDI"}, hdi_prob=0.95)
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"50% HDI", "alpha": 1}, hdi_prob=0.50)
    # Add labels and title
    # Apply the custom formatter to the y-axis
    plt.gca().yaxis.set_major_formatter(FuncFormatter(multiply_by_100))
    plt.xlabel("Year")
    plt.ylabel("Mean")
    plt.title("Mean with HDI per Year")
    plt.legend()
    ## ----end
    plt.savefig(f"{paths['fig_path']}python_pdp_3_mod_pymc_bart_16b.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_mse_16b_3_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_16_"], "rb"))
    ppc4 = trace_dict["ppc4"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_16_"], "rb"))
    newdata_15b = data_dict["newdata_15b"]
    ## ---- python sampled data pymc_barts mse 16b 3 
    posterior_predictive = ppc4.posterior_predictive["y"]
    mse = (newdata_15b["HCC"] - posterior_predictive.mean(dim=["draw", "chain"]))**2
    acc = np.abs(np.exp(np.log(newdata_15b["HCC"]) - np.log(posterior_predictive.mean(dim=["draw", "chain"]))) - 1)
    df = pd.DataFrame({
        "mse_mean": [np.mean(mse)],
        "mse_median": [np.median(mse)],
        "acc_mean": [np.mean(acc)],
        "acc_median": [np.median(acc)],
        "model": ["pymc-bart"],
        "lower": [np.quantile(mse, 0.025)],
        "upper": [np.quantile(mse, 0.975)],
        "type": 3,
        "model_type": "covariates"
    })
    df.to_csv(f"{paths['data_path']}/modelled/pymc_bart_16b_mse_3.csv", index=False)
    ## ----end
    results=1
    pickle.dump(results, open(product, "wb"))

## Eastern sampled reefs --------------------------------------------------

def sampled_reefs_bart_data_17_(product, upstream):
    ## This function prepares the data for the BART model
    ## It prepares two datasets:
    ## 1. `sample_data`: The sampled data (25 sites (etc), 17 years) that
    ##    will be used to fit the model.
    ## 2. `full_data`: The full dataset (all sites, all years) that will
    ##    be used to make predictions on the scale of the full
    ##    spatio-temporal domain.
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_17 = pickle.load(open(upstream["incomplete_spatial_data_prep_17_"], "rb"))
    newdata_14b = pickle.load(open(upstream["incomplete_spatial_newdata_14b_"], "rb"))
    newdata_15b = pickle.load(open(upstream["incomplete_spatial_newdata_15b_"], "rb"))
    newdata_17b = pickle.load(open(upstream["incomplete_spatial_newdata_17b_"], "rb"))
    newdata_17b = pickle.load(open(upstream["incomplete_spatial_newdata_17b_"], "rb"))
    ## ---- python sampled data barts data 17
    upstream = {"incomplete_spatial_data_prep_17_": "../python/cache/incomplete_spatial_data_prep_17_.pkl",
                "incomplete_spatial_global_parameters_": "../python/cache/incomplete_spatial_global_parameters_.pkl",
                "incomplete_spatial_newdata_14b_": "../python/cache/incomplete_spatial_newdata_14b_.pkl",
                "incomplete_spatial_newdata_15b_": "../python/cache/incomplete_spatial_newdata_15b_.pkl",
                "incomplete_spatial_newdata_17b_": "../python/cache/incomplete_spatial_newdata_17b_.pkl",
                "incomplete_spatial_newdata_17b_": "../python/cache/incomplete_spatial_newdata_17b_.pkl"}


    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_17 = pickle.load(open(upstream["incomplete_spatial_data_prep_17_"], "rb"))
    newdata_17b = pickle.load(open(upstream["incomplete_spatial_newdata_17b_"], "rb"))
    # Prepare the training data 
    benthos_fixed_locs_obs_17_prep = benthos_fixed_locs_obs_17.copy()
    for col in ["Reef", "Site", "Transect"]:
        benthos_fixed_locs_obs_17_prep[col] = benthos_fixed_locs_obs_17_prep[col].astype("category").cat.remove_unused_categories()
    
    sample_data = benthos_fixed_locs_obs_17_prep
    # prepare the northern reefs prediction data (newdata_17b)
    newdata_17b = newdata_17b
    newdata_17b = newdata_17b
    ## ----end
    data_dict = {
        "sample_data": sample_data,
        "newdata_14b": newdata_14b,
        "newdata_15b": newdata_15b,
        "newdata_17b": newdata_17b,
        "newdata_17b": newdata_17b
    }
    pickle.dump(data_dict, open(product, "wb"))

def sampled_reefs_bart_data_17_sample_data_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_17 = pickle.load(open(upstream["incomplete_spatial_data_prep_17_"], "rb"))
    ## ---- python sampled data barts data 17 sample data
    upstream = {"incomplete_spatial_data_prep_17_": "../python/cache/incomplete_spatial_data_prep_17_.pkl",
                "incomplete_spatial_global_parameters_": "../python/cache/incomplete_spatial_global_parameters_.pkl"
                }

    benthos_fixed_locs_obs_17 = pickle.load(open(upstream["incomplete_spatial_data_prep_17_"], "rb"))
    benthos_fixed_locs_obs_17
    benthos_fixed_locs_obs_17.groupby("Year").agg(
        Mean=("HCC", "mean"),
        Median=("HCC", "median")
    ).reset_index()
    ## ----end
    result = 1
    pickle.dump(result, open(product, "wb"))

def sampled_reefs_bart_data_17_newdata_17b_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    newdata_17b = pickle.load(open(upstream["incomplete_spatial_newdata_17b_"], "rb"))
    ## ---- python sampled data barts data 17 newdata 17b
    upstream = {"incomplete_spatial_newdata_17b_": "../python/cache/incomplete_spatial_newdata_17b_.pkl",
                "incomplete_spatial_global_parameters_": "../python/cache/incomplete_spatial_global_parameters_.pkl"
                }
    newdata_17b = pickle.load(open(upstream["incomplete_spatial_newdata_17b_"], "rb"))
    newdata_17b
    newdata_17b.groupby("Year").agg(
        Mean=("HCC", "mean"),
        Median=("HCC", "median")
    ).reset_index()
    ## ----end
    result = 1
    pickle.dump(result, open(product, "wb"))

def sampled_reefs_bart_data_17_newdata_17b_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    newdata_17b = pickle.load(open(upstream["incomplete_spatial_newdata_17b_"], "rb"))
    ## ---- python sampled data barts data 17 newdata 17b
    upstream = {"incomplete_spatial_newdata_17b_": "../python/cache/incomplete_spatial_newdata_17b_.pkl",
                "incomplete_spatial_global_parameters_": "../python/cache/incomplete_spatial_global_parameters_.pkl"
                }
    newdata_17b = pickle.load(open(upstream["incomplete_spatial_newdata_17b_"], "rb"))
    newdata_17b
    newdata_17b.groupby("Year").agg(
        Mean=("HCC", "mean"),
        Median=("HCC", "median")
    ).reset_index()
    ## ----end
    result = 1
    pickle.dump(result, open(product, "wb"))

def sampled_reefs_bart_fit_17_(product, upstream):
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_17_"], "rb"))
    sample_data = data_dict["sample_data"]
    newdata_14b = data_dict["newdata_14b"]
    newdata_15b = data_dict["newdata_15b"]
    newdata_17b = data_dict["newdata_17b"]
    newdata_17b = data_dict["newdata_17b"]
    ## ---- python sampled data barts fit 17
    X = sample_data[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = sample_data["cover"]
    with pm.Model() as model_17:
        data_X = pm.Data("data_X", X)
        data_Y = pm.Data("data_Y", Y)
        μ = pmb.BART("μ", data_X, Y, m=100)
        alpha = pm.HalfNormal("alpha", sigma=2)
        beta = pm.Deterministic("beta", alpha * (1 - μ) / μ)
        y = pm.Beta("y", alpha=alpha, beta=beta, observed=data_Y)
    with model_17:
        idata = pm.sample(random_seed=173)
    with model_17:
        ppc1 = pm.sample_posterior_predictive(idata)

    another_X = newdata_17b[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]].copy()
    X = another_X[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = newdata_17b["HCC"]
    with model_17:
            pm.set_data({"data_X": another_X,
                        "data_Y": np.arange(another_X.shape[0])
                        })
            ppc2 = pm.sample_posterior_predictive(idata)

    another_X = newdata_15b[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]].copy()
    X = another_X[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = newdata_15b["HCC"]
    with model_17:
            pm.set_data({"data_X": another_X,
                        "data_Y": np.arange(another_X.shape[0])
                        })
            ppc3 = pm.sample_posterior_predictive(idata)

    another_X = newdata_14b[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]].copy()
    X = another_X[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = newdata_14b["HCC"]
    with model_17:
            pm.set_data({"data_X": another_X,
                        "data_Y": np.arange(another_X.shape[0])
                        })
            ppc4 = pm.sample_posterior_predictive(idata)
    ## ----end
    trace_dict = {
        "idata": idata,
        "ppc1": ppc1,
        "ppc2": ppc2,
        "ppc3": ppc3,
        "ppc4": ppc4
    }
    pickle.dump(trace_dict, open(product, "wb"))

def sampled_reefs_bart_prep_plot_17_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_17_"], "rb"))
    ppc1 = trace_dict["ppc1"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_17_"], "rb"))
    sample_data = data_dict["sample_data"]
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ---- python sampled data pymc_barts prep plot 17 
    posterior_predictive = ppc1.posterior_predictive["y"]  # replace "y_pred" with your variable name
    year_array = sample_data["Year"].values  # must match the shape of predictions

    # Add 'Year' as a coordinate to the xarray
    posterior_predictive = posterior_predictive.assign_coords({"y_dim_0": year_array})
    # Group by Year and compute mean per chain/draw
    avg_preds = posterior_predictive.groupby("y_dim_0").mean("y_dim_0")
    # Rename dimension for clarity
    avg_preds = avg_preds.rename({"y_dim_0": "Year"})
    modelled_trend = az.summary(avg_preds, hdi_prob=0.95)
    modelled_trend["Year"] = avg_preds.coords["Year"].values

    ## Summarise the sampled data
    sample_data_summary = (
        sample_data.groupby("Year").agg(
            Mean=("HCC", "mean"),
            Median=("HCC", "median")
        ).reset_index()
    )

    ## Summarise the full grid
    benthos_reefs_temporal_summary = (
        benthos_reefs_sf.groupby("Year").agg(
            Mean=("HCC", "mean"),
            Median=("HCC", "median"),
            SD=("HCC", "std"),
            Lower=("HCC", lambda x: x.quantile(0.025)),
            Upper=("HCC", lambda x: x.quantile(0.975))
        ).reset_index()
    )
    years = np.unique(year_array)
    cover_mean = avg_preds.mean(dim=["draw", "chain"])
    cover = avg_preds
    ## ----end
    summary_dict = {
        "modelled_trend": modelled_trend,
        "sample_data_summary": sample_data_summary,
        "benthos_reefs_temporal_summary": benthos_reefs_temporal_summary,
        "years": years,
        "cover_mean": cover_mean,
        "cover": cover
    }
    pickle.dump(summary_dict, open(product, "wb"))

def sampled_reefs_bart_plot_17_(product, upstream):
    # paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_17_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    benthos_reefs_temporal_summary = summary_dict["benthos_reefs_temporal_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python sampled data pymc_barts plot 17 
    def multiply_by_100(x, pos):
        return f"{x * 100:.0f}"
    plt.figure(figsize=(8, 6))
    plt.plot(years, cover_mean, "w", lw=3)
    ## Add true mean line
    #plt.plot(
    #    benthos_reefs_temporal_summary["Year"],
    #    benthos_reefs_temporal_summary["Mean"],
    #    label="True mean",
    #    color="blue"
    #)
    ## Add true median line
    #plt.plot(
    #    benthos_reefs_temporal_summary["Year"],
    #    benthos_reefs_temporal_summary["Median"],
    #    label="True median",
    #    color="green"
    #)
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Mean"],
        label="Sample mean",
        color="blue",
        linestyle="--"
    )
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Median"],
        label="Sample median",
        color="green",
        linestyle="--"
    )
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"95% HDI"}, hdi_prob=0.95)
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"50% HDI", "alpha": 1}, hdi_prob=0.50)
    # Add labels and title
    # Apply the custom formatter to the y-axis
    plt.gca().yaxis.set_major_formatter(FuncFormatter(multiply_by_100))
    plt.xlabel("Year")
    plt.ylabel("Mean")
    plt.title("Mean with HDI per Year")
    plt.legend()
    ## ----end
    plt.savefig(f"{paths['fig_path']}python_pdp_1_mod_pymc_bart_17b.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_mse_17b_1_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_17_"], "rb"))
    ppc2 = trace_dict["ppc2"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_17_"], "rb"))
    newdata_17b = data_dict["newdata_17b"]
    ## ---- python sampled data pymc_barts mse 17b 1 
    posterior_predictive = ppc2.posterior_predictive["y"]
    mse = (newdata_17b["HCC"] - posterior_predictive.mean(dim=["draw", "chain"]))**2
    acc = np.abs(np.exp(np.log(newdata_17b["HCC"]) - np.log(posterior_predictive.mean(dim=["draw", "chain"]))) - 1)
    df = pd.DataFrame({
        "mse_mean": [np.mean(mse)],
        "mse_median": [np.median(mse)],
        "acc_mean": [np.mean(acc)],
        "acc_median": [np.median(acc)],
        "model": ["pymc-bart"],
        "lower": [np.quantile(mse, 0.025)],
        "upper": [np.quantile(mse, 0.975)],
        "type": 1,
        "model_type": "covariates"
    })
    df.to_csv(f"{paths['data_path']}/modelled/pymc_bart_17b_mse_1.csv", index=False)
    ## ----end
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_prep_plot_17_2_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_17_"], "rb"))
    ppc3 = trace_dict["ppc3"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_17_"], "rb"))
    ## sample_data = data_dict["sample_data"]
    ## benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    newdata_15b = data_dict["newdata_15b"]
    ## ---- python sampled data pymc_barts prep plot 17 2
    posterior_predictive = ppc3.posterior_predictive["y"]  # replace "y_pred" with your variable name
    year_array = newdata_15b["Year"].values  # must match the shape of predictions

    # Add 'Year' as a coordinate to the xarray
    posterior_predictive = posterior_predictive.assign_coords({"y_dim_0": year_array})
    # Group by Year and compute mean per chain/draw
    avg_preds = posterior_predictive.groupby("y_dim_0").mean("y_dim_0")
    # Rename dimension for clarity
    avg_preds = avg_preds.rename({"y_dim_0": "Year"})
    modelled_trend = az.summary(avg_preds, hdi_prob=0.95)
    modelled_trend["Year"] = avg_preds.coords["Year"].values

    ## Summarise the sampled data
    sample_data_summary = (
        newdata_15b.groupby("Year").agg(
            Mean=("HCC", "mean"),
            Median=("HCC", "median")
        ).reset_index()
    )

    years = np.unique(year_array)
    cover_mean = avg_preds.mean(dim=["draw", "chain"])
    cover = avg_preds
    ## ----end
    summary_dict = {
        "modelled_trend": modelled_trend,
        "sample_data_summary": sample_data_summary,
        "years": years,
        "cover_mean": cover_mean,
        "cover": cover
    }
    pickle.dump(summary_dict, open(product, "wb"))

def sampled_reefs_bart_plot_17_2_(product, upstream):
    # paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_17_2_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python sampled data pymc_barts plot 17 2
    def multiply_by_100(x, pos):
        return f"{x * 100:.0f}"
    plt.figure(figsize=(8, 6))
    plt.plot(years, cover_mean, "w", lw=3)
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Mean"],
        label="Sample mean",
        color="blue",
        linestyle="--"
    )
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Median"],
        label="Sample median",
        color="green",
        linestyle="--"
    )
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"95% HDI"}, hdi_prob=0.95)
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"50% HDI", "alpha": 1}, hdi_prob=0.50)
    # Add labels and title
    # Apply the custom formatter to the y-axis
    plt.gca().yaxis.set_major_formatter(FuncFormatter(multiply_by_100))
    plt.xlabel("Year")
    plt.ylabel("Mean")
    plt.title("Mean with HDI per Year")
    plt.legend()
    ## ----end
    plt.savefig(f"{paths['fig_path']}python_pdp_2_mod_pymc_bart_17b.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_mse_17b_2_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_17_"], "rb"))
    ppc3 = trace_dict["ppc3"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_17_"], "rb"))
    newdata_15b = data_dict["newdata_15b"]
    ## ---- python sampled data pymc_barts mse 17b 2 
    posterior_predictive = ppc3.posterior_predictive["y"]
    mse = (newdata_15b["HCC"] - posterior_predictive.mean(dim=["draw", "chain"]))**2
    acc = np.abs(np.exp(np.log(newdata_15b["HCC"]) - np.log(posterior_predictive.mean(dim=["draw", "chain"]))) - 1)
    df = pd.DataFrame({
        "mse_mean": [np.mean(mse)],
        "mse_median": [np.median(mse)],
        "acc_mean": [np.mean(acc)],
        "acc_median": [np.median(acc)],
        "model": ["pymc-bart"],
        "lower": [np.quantile(mse, 0.025)],
        "upper": [np.quantile(mse, 0.975)],
        "type": 2,
        "model_type": "covariates"
    })
    df.to_csv(f"{paths['data_path']}/modelled/pymc_bart_17b_mse_2.csv", index=False)
    ## ----end
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_prep_plot_17_3_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_17_"], "rb"))
    ppc4 = trace_dict["ppc4"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_17_"], "rb"))
    ## sample_data = data_dict["sample_data"]
    ## benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    newdata_14b = data_dict["newdata_14b"]
    ## ---- python sampled data pymc_barts prep plot 17 3
    posterior_predictive = ppc4.posterior_predictive["y"]  # replace "y_pred" with your variable name
    year_array = newdata_14b["Year"].values  # must match the shape of predictions

    # Add 'Year' as a coordinate to the xarray
    posterior_predictive = posterior_predictive.assign_coords({"y_dim_0": year_array})
    # Group by Year and compute mean per chain/draw
    avg_preds = posterior_predictive.groupby("y_dim_0").mean("y_dim_0")
    # Rename dimension for clarity
    avg_preds = avg_preds.rename({"y_dim_0": "Year"})
    modelled_trend = az.summary(avg_preds, hdi_prob=0.95)
    modelled_trend["Year"] = avg_preds.coords["Year"].values

    ## Summarise the sampled data
    sample_data_summary = (
        newdata_14b.groupby("Year").agg(
            Mean=("HCC", "mean"),
            Median=("HCC", "median")
        ).reset_index()
    )

    years = np.unique(year_array)
    cover_mean = avg_preds.mean(dim=["draw", "chain"])
    cover = avg_preds
    ## ----end
    summary_dict = {
        "modelled_trend": modelled_trend,
        "sample_data_summary": sample_data_summary,
        "years": years,
        "cover_mean": cover_mean,
        "cover": cover
    }
    pickle.dump(summary_dict, open(product, "wb"))

def sampled_reefs_bart_plot_17_3_(product, upstream):
    # paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_17_3_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python sampled data pymc_barts plot 17 3
    def multiply_by_100(x, pos):
        return f"{x * 100:.0f}"
    plt.figure(figsize=(8, 6))
    plt.plot(years, cover_mean, "w", lw=3)
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Mean"],
        label="Sample mean",
        color="blue",
        linestyle="--"
    )
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Median"],
        label="Sample median",
        color="green",
        linestyle="--"
    )
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"95% HDI"}, hdi_prob=0.95)
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"50% HDI", "alpha": 1}, hdi_prob=0.50)
    # Add labels and title
    # Apply the custom formatter to the y-axis
    plt.gca().yaxis.set_major_formatter(FuncFormatter(multiply_by_100))
    plt.xlabel("Year")
    plt.ylabel("Mean")
    plt.title("Mean with HDI per Year")
    plt.legend()
    ## ----end
    plt.savefig(f"{paths['fig_path']}python_pdp_3_mod_pymc_bart_17b.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_mse_17b_3_(product, upstream):
    paths = pickle.load(open(upstream["incomplete_spatial_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_17_"], "rb"))
    ppc4 = trace_dict["ppc4"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_17_"], "rb"))
    newdata_14b = data_dict["newdata_14b"]
    ## ---- python sampled data pymc_barts mse 17b 3 
    posterior_predictive = ppc4.posterior_predictive["y"]
    mse = (newdata_14b["HCC"] - posterior_predictive.mean(dim=["draw", "chain"]))**2
    acc = np.abs(np.exp(np.log(newdata_14b["HCC"]) - np.log(posterior_predictive.mean(dim=["draw", "chain"]))) - 1)
    df = pd.DataFrame({
        "mse_mean": [np.mean(mse)],
        "mse_median": [np.median(mse)],
        "acc_mean": [np.mean(acc)],
        "acc_median": [np.median(acc)],
        "model": ["pymc-bart"],
        "lower": [np.quantile(mse, 0.025)],
        "upper": [np.quantile(mse, 0.975)],
        "type": 3,
        "model_type": "covariates"
    })
    df.to_csv(f"{paths['data_path']}/modelled/pymc_bart_17b_mse_3.csv", index=False)
    ## ----end
    results=1
    pickle.dump(results, open(product, "wb"))



