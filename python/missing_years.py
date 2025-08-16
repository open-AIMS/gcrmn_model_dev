
# + tags=["parameters"]
# declare a list tasks whose products you want to use as inputs
upstream = None
# -

import multiprocessing

## ---- python missing years libraries
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

def missing_years_global_parameters_(product):
    ## ---- python missing years global parameters
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

def missing_years_read_all_reefs_data_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    ## ---- python missing years read all reefs data
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ----end
    pickle.dump(benthos_reefs_sf, open(product, "wb"))

def missing_years_read_all_reefs_data_show_(product, upstream):
    benthos_reefs_sf = pickle.load(open(upstream["missing_years_read_all_reefs_data_"], "rb"))
    ## ---- python missing years read all reefs data show
    pd.set_option("display.max_columns", None)  # Show all columns
    benthos_reefs_sf 
    pd.reset_option("display.max_columns")
    ## ----end
    pickle.dump("done", open(product, "wb"))

def missing_years_read_all_temporal_summary_(product, upstream):
    benthos_reefs_sf = pickle.load(open(upstream["missing_years_read_all_reefs_data_"], "rb"))
    ## ---- python all reefs temporal summary
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

def missing_years_read_all_temporal_summary_plot_(product, upstream):
    benthos_reefs_temporal_summary = pickle.load(open(upstream["missing_years_read_all_temporal_summary_"], "rb"))
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    ## ---- python all reefs temporal summary plot
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

## All sampled reefs --------------------------------------------------

def missing_years_read_sampled_reefs_data_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    ## ---- python missing years read sampled reefs data
    benthos_fixed_locs_obs = pd.read_csv(f"{paths['data_path']}synthetic/benthos_fixed_locs_obs.csv")
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    # Group and summarize the benthos_reefs_sf DataFrame
    benthos_reefs_summary = (
        benthos_reefs_sf
        .drop(columns=["geometry"])  # Equivalent to st_drop_geometry()
        .groupby(["Year", "Reef"], as_index=False)
        .agg({
            "CYC": "mean",
            "DHW": "mean",
            "OTHER": "mean"
        })
    )
    # Perform the left join
    benthos_fixed_locs_obs = benthos_fixed_locs_obs.merge(
        benthos_reefs_summary,
        on=["Year", "Reef"],
        how="left"
    )
    ## ----end
    pickle.dump(benthos_fixed_locs_obs, open(product, "wb"))

def missing_years_read_sampled_reefs_data_show_(product, upstream):
    benthos_fixed_locs_obs = pickle.load(open(upstream["missing_years_read_sampled_reefs_data_"], "rb"))
    ## ---- python missing years read sampled reefs data show
    pd.set_option("display.max_columns", None)  # Show all columns
    benthos_fixed_locs_obs
    pd.reset_option("display.max_columns")
    ## ----end
    pickle.dump("done", open(product, "wb"))

def missing_years_sampled_simple_raw_means_(product, upstream):
    benthos_fixed_locs_obs = pickle.load(open(upstream["missing_years_read_sampled_reefs_data_"], "rb"))
    ## ---- python missing years sampled simple raw means
    benthos_fixed_locs_obs = benthos_fixed_locs_obs.assign(
        fYear = benthos_fixed_locs_obs['Year'].astype('category'),
        Reef = benthos_fixed_locs_obs['Reef'].astype('category'),
        Site = benthos_fixed_locs_obs['Reef'].astype(str) + "_" +
        benthos_fixed_locs_obs['Site'].astype(str),
        Transect = benthos_fixed_locs_obs['Site'] + "_" +
        benthos_fixed_locs_obs['Transect'].astype(str),
        cover = benthos_fixed_locs_obs['HCC'] / 100
    )
    all_sampled_sum = (
        benthos_fixed_locs_obs.groupby("fYear").agg(
            mean_mean_response = ("cover", "mean"),
            median_median_response = ("cover", "median"),
            mean_sd = ("cover", "std"),
            median_lower_lower = ("cover", lambda x: x.quantile(0.025)),
            median_upper_upper = ("cover", lambda x: x.quantile(0.975))
        ).reset_index()
    )
    # Calculate additional columns
    all_sampled_sum['mean_lower_lower'] = all_sampled_sum['mean_mean_response'] - 1.96 * all_sampled_sum['mean_sd']
    all_sampled_sum['mean_upper_upper'] = all_sampled_sum['mean_mean_response'] + 1.96 * all_sampled_sum['mean_sd']
    # Drop the 'mean_sd' column
    all_sampled_sum = all_sampled_sum.drop(columns=['mean_sd'])
    # Pivot longer
    all_sampled_sum = all_sampled_sum.melt(
        id_vars=['fYear'],
        value_vars=[
            'mean_mean_response', 'median_median_response', 'mean_lower_lower',
            'mean_upper_upper', 'median_lower_lower', 'median_upper_upper'
        ],
        var_name='names',
        value_name='values'
    )
    # Split the 'names' column into 'type', 'variable', and 'stat'
    all_sampled_sum[['type', 'variable', 'stat']] = all_sampled_sum['names'].str.split('_', expand=True)
    # Modify the 'type' column
    all_sampled_sum['type'] = 'all_sampled_' + all_sampled_sum['type']
    # Drop the 'variable' column
    all_sampled_sum = all_sampled_sum.drop(columns=['variable', 'names'])
    # Pivot wider
    all_sampled_sum = all_sampled_sum.pivot(
        index=['fYear', 'type'],
        columns='stat',
        values='values'
    ).reset_index()
    # Convert 'fYear' to numeric
    all_sampled_sum['Year'] = pd.to_numeric(all_sampled_sum['fYear'])
    all_sampled_sum
    ## ----end
    pickle.dump(all_sampled_sum, open(product, "wb"))

def missing_years_sampled_simple_raw_means_plot_(product, upstream):
    all_sampled_sum = pickle.load(open(upstream["missing_years_sampled_simple_raw_means_"], "rb"))
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    ## ---- python missing years sampled simple raw means plot
    # Create the plot
    plt.figure(figsize=(10, 6))
    # Add ribbon (shaded area)
    for _, group in all_sampled_sum.groupby('type'):
        plt.fill_between(
            group['Year'],
            group['lower'],
            group['upper'],
            alpha=0.2,
            label=_
        )
    # Add mean line
    for _, group in all_sampled_sum.groupby('type'):
        plt.plot(
            group["Year"],
            group["response"],
            label=_
        )
    # Customize the y-axis
    plt.ylim(0, 1)
    plt.ylabel("Hard Coral Cover (%)")
    plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x*100:.0f}%"))
    # Add theme and legend
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.legend(title='Type')
    plt.title("Hard Coral Cover Over Time")
    plt.tight_layout()
    # Show the plot
    #plt.show()
    # Save the plot as a PNG file
    plt.savefig(f"{paths['fig_path']}Python_full_simple_raw_means_plot2.png",
                dpi=72, bbox_inches="tight")
    ## ----end
    plt.savefig(product, dpi=72, bbox_inches="tight")

## Temporal gap of some reefs -----------------------------------------

def read_sampled_reefs_data_3_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    ## ---- python read sampled reefs data 3
    benthos_fixed_locs_obs_3 = pd.read_csv(f"{paths['data_path']}synthetic/benthos_fixed_locs_obs_3.csv")
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    # Group and summarize the benthos_reefs_sf DataFrame
    benthos_reefs_summary = (
        benthos_reefs_sf
        .drop(columns=["geometry"])  # Equivalent to st_drop_geometry()
        .groupby(["Year", "Reef"], as_index=False)
        .agg({
            "CYC": "mean",
            "DHW": "mean",
            "OTHER": "mean"
        })
    )
    # Perform the left join
    benthos_fixed_locs_obs_3 = benthos_fixed_locs_obs_3.merge(
        benthos_reefs_summary,
        on=["Year", "Reef"],
        how="left"
    )
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_3, open(product, "wb"))

def read_sampled_reefs_data_3_show_(product, upstream):
    benthos_fixed_locs_obs_3 = pickle.load(open(upstream["read_sampled_reefs_data_3_"], "rb"))
    ## ---- python read sampled reefs data 3 show
    pd.set_option("display.max_columns", None)  # Show all columns
    benthos_fixed_locs_obs_3
    pd.reset_option("display.max_columns")
    ## ----end
    pickle.dump("done", open(product, "wb"))

def sampled_reefs_data_3_plot_(product, upstream):
    benthos_fixed_locs_obs_3 = pickle.load(open(upstream["read_sampled_reefs_data_3_"], "rb"))
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    ## ---- python sampled reefs data 3 plot
    # Create the plot
    all_years = pd.DataFrame({"Year": range(benthos_fixed_locs_obs_3["Year"].min(), benthos_fixed_locs_obs_3["Year"].max() + 1)})
    data_with_gaps = (
        benthos_fixed_locs_obs_3.groupby(["Reef", "Site", "Transect"])
        .apply(lambda group: group.set_index("Year").reindex(all_years["Year"]).reset_index())
        .reset_index(drop=True)
    )
    # Fill missing columns with appropriate values
    data_with_gaps["Reef"] = data_with_gaps["Reef"].fillna(method="ffill")
    data_with_gaps["Site"] = data_with_gaps["Site"].fillna(method="ffill")
    data_with_gaps["Transect"] = data_with_gaps["Transect"].fillna(method="ffill")
    data_with_gaps["Trans"] = data_with_gaps["Site"] + "_" + data_with_gaps["Transect"]

    g = sns.FacetGrid(
        data=data_with_gaps,
        col="Reef",
        col_wrap=5,  # Adjust the number of columns in the facet grid
        height=2,
        aspect=1.6,
        sharey=True,
        sharex=True
    )
    
    # Add line plots to each facet.
    ## Note, have to fudge a pointplot because a lineplot silently drops nan values
    g.map_dataframe(
        sns.pointplot,
        "Year",
        "HCC",
        markers=None,
        hue="Trans",
        errorbar=None,
        palette=[color for color in ["#1f77b4", "#ff7f0e"] for _ in range(5)],
        units="Transect"
    )

    # Customize the plot
    g.set_axis_labels("Year", "Hard Coral Cover (HCC)")
    g.set_titles(col_template="{col_name}")
    # g.add_legend(title="Site")
    g.fig.tight_layout()

    # Show the plot
    #plt.show()
    
    # Save the plot as a PNG file
    plt.savefig(f"{paths['fig_path']}Python_sampled_reefs_3_plot.png", dpi=72, bbox_inches="tight")
    ## ----end
    plt.savefig(product, dpi=72, bbox_inches="tight")

## Temporal gap of all reefs -----------------------------------------

def read_sampled_reefs_data_4_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    ## ---- python read sampled reefs data 4
    benthos_fixed_locs_obs_4 = pd.read_csv(f"{paths['data_path']}synthetic/benthos_fixed_locs_obs_4.csv")
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    # Group and summarize the benthos_reefs_sf DataFrame
    benthos_reefs_summary = (
        benthos_reefs_sf
        .drop(columns=["geometry"])  # Equivalent to st_drop_geometry()
        .groupby(["Year", "Reef"], as_index=False)
        .agg({
            "CYC": "mean",
            "DHW": "mean",
            "OTHER": "mean"
        })
    )
    # Perform the left join
    benthos_fixed_locs_obs_4 = benthos_fixed_locs_obs_4.merge(
        benthos_reefs_summary,
        on=["Year", "Reef"],
        how="left"
    )
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_4, open(product, "wb"))

def read_sampled_reefs_data_4_show_(product, upstream):
    benthos_fixed_locs_obs_4 = pickle.load(open(upstream["read_sampled_reefs_data_4_"], "rb"))
    ## ---- python read sampled reefs data 4 show
    pd.set_option("display.max_columns", None)  # Show all columns
    benthos_fixed_locs_obs_4
    pd.reset_option("display.max_columns")
    ## ----end
    pickle.dump("done", open(product, "wb"))

def sampled_reefs_data_4_plot_(product, upstream):
    benthos_fixed_locs_obs_4 = pickle.load(open(upstream["read_sampled_reefs_data_4_"], "rb"))
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    ## ---- python sampled reefs data 4 plot
    # Create the plot
    all_years = pd.DataFrame({"Year": range(benthos_fixed_locs_obs_4["Year"].min(), benthos_fixed_locs_obs_4["Year"].max() + 1)})
    data_with_gaps = (
        benthos_fixed_locs_obs_4.groupby(["Reef", "Site", "Transect"])
        .apply(lambda group: group.set_index("Year").reindex(all_years["Year"]).reset_index())
        .reset_index(drop=True)
    )
    # Fill missing columns with appropriate values
    data_with_gaps["Reef"] = data_with_gaps["Reef"].fillna(method="ffill")
    data_with_gaps["Site"] = data_with_gaps["Site"].fillna(method="ffill")
    data_with_gaps["Transect"] = data_with_gaps["Transect"].fillna(method="ffill")
    data_with_gaps["Trans"] = data_with_gaps["Site"] + "_" + data_with_gaps["Transect"]

    g = sns.FacetGrid(
        data=data_with_gaps,
        col="Reef",
        col_wrap=5,  # Adjust the number of columns in the facet grid
        height=2,
        aspect=1.6,
        sharey=True,
        sharex=True
    )
    
    # Add line plots to each facet.
    ## Note, have to fudge a pointplot because a lineplot silently drops nan values
    g.map_dataframe(
        sns.pointplot,
        "Year",
        "HCC",
        markers=None,
        hue="Trans",
        errorbar=None,
        palette=[color for color in ["#1f77b4", "#ff7f0e"] for _ in range(5)],
        units="Transect"
    )

    # Customize the plot
    g.set_axis_labels("Year", "Hard Coral Cover (HCC)")
    g.set_titles(col_template="{col_name}")
    # g.add_legend(title="Site")
    g.fig.tight_layout()

    # Show the plot
    #plt.show()
    
    # Save the plot as a PNG file
    plt.savefig(f"{paths['fig_path']}Python_sampled_reefs_4_plot.png", dpi=72, bbox_inches="tight")
    ## ----end
    plt.savefig(product, dpi=72, bbox_inches="tight")

## Data transformations ============================================

## All sampled reefs --------------------------------------------------

## All reefs -----------------------------------------

def missing_years_benthos_reefs_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    benthos_reefs_sf = pickle.load(open(upstream["missing_years_read_all_reefs_data_"], "rb"))
    ## ---- python benthos reefs
    pd.set_option("display.max_columns", None)  # Show all columns
    benthos_reefs = benthos_reefs_sf.copy()
    benthos_reefs[['Longitude', 'Latitude']] = (
      benthos_reefs['geometry']
      .str.strip('c()')
      .str.split(',', expand = True)
      .astype(float)
    )
    # newdata_b = (
    #   benthos_fixed_locs_obs_0.groupby(["Year", "Reef"], as_index=False, observed=True).agg(
    #     HCC = ("HCC", "mean"),
    #     CYC = ("CYC", "mean"),
    #     DHW = ("DHW", "mean"),
    #     OTHER = ("OTHER", "mean"),
    #     Longitude = ("Longitude", "mean"),
    #     Latitude = ("Latitude", "mean")
    #   ).reset_index()
    # )
    benthos_reefs
    ## ----end
    pickle.dump(benthos_reefs, open(product, "wb"))

def missing_years_newdata_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    benthos_reefs = pickle.load(open(upstream["missing_years_benthos_reefs_"], "rb"))
    ## ---- python newdata
    newdata = pd.MultiIndex.from_product(
        [benthos_reefs["Year"].unique()],
        names=["Year"]
    ).to_frame(index=False)
    
    newdata
    ## ----end
    pickle.dump(newdata, open(product, "wb"))

def missing_years_newdata_b_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    benthos_reefs = pickle.load(open(upstream["missing_years_benthos_reefs_"], "rb"))
    ## ---- python newdata b
    pd.set_option("display.max_columns", None)  # Show all columns
    newdata_b = (
      benthos_reefs.groupby(["Year", "Reef"], as_index=False, observed=True).agg(
        HCC = ("HCC", "mean"),
        CYC = ("CYC", "mean"),
        DHW = ("DHW", "mean"),
        OTHER = ("OTHER", "mean"),
        Longitude = ("Longitude", "mean"),
        Latitude = ("Latitude", "mean")
      ).reset_index()
    )
    newdata_b
    ## ----end
    pickle.dump(newdata_b, open(product, "wb"))

## All sampled reefs -----------------------------------------

def missing_years_sampled_reefs_data_prep_0_(product, upstream):
    benthos_fixed_locs_obs = pickle.load(open(upstream["missing_years_read_sampled_reefs_data_"], "rb"))
    ## ---- python missing years sampled data prep 0
    benthos_fixed_locs_obs_0 = (
        benthos_fixed_locs_obs.assign(
            fYear=lambda df: df["Year"].astype("category"),
            Reef=lambda df: df["Reef"].astype("category"),
            Site=lambda df: (df["Reef"].astype("str") + "_" +
                             df["Site"]).astype("category"),  # Interaction of Reef and Site
            Transect=lambda df: (df["Site"].astype("str") + "_" +
                                 df["Transect"]).astype("category"),  # Interaction of Site and Transect
            cover=lambda df: df["HCC"] / 100  # Calculate cover
        )
    )
    benthos_fixed_locs_obs_0 
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_0, open(product, "wb"))

def missing_years_newdata_0_(product, upstream):
    benthos_fixed_locs_obs_0 = pickle.load(open(upstream["missing_years_sampled_reefs_data_prep_0_"], "rb"))
    ## ---- python missing years newdata 0
    newdata_0 = pd.MultiIndex.from_product(
        [benthos_fixed_locs_obs_0["fYear"].unique()],
        names=["fYear"]
    ).to_frame(index=False)
    newdata_0
    ## ----end
    pickle.dump(newdata_0, open(product, "wb"))

## Temporal gap of some reefs -----------------------------------------

def sampled_reefs_data_prep_3_(product, upstream):
    benthos_fixed_locs_obs_3 = pickle.load(open(upstream["read_sampled_reefs_data_3_"], "rb"))
    ## ---- python sampled data prep 3
    benthos_fixed_locs_obs_3 = (
        benthos_fixed_locs_obs_3.assign(
            fYear=lambda df: df["Year"].astype("category"),
            Reef=lambda df: df["Reef"].astype("category"),
            Site=lambda df: (df["Reef"].astype("str") + "_" +
                             df["Site"]).astype("category"),  # Interaction of Reef and Site
            Transect=lambda df: (df["Site"].astype("str") + "_" +
                                 df["Transect"]).astype("category"),  # Interaction of Site and Transect
            cover=lambda df: df["HCC"] / 100  # Calculate cover
        )
    )
    benthos_fixed_locs_obs_3
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_3, open(product, "wb"))

def newdata_3_(product, upstream):
    benthos_fixed_locs_obs_3 = pickle.load(open(upstream["sampled_reefs_data_prep_3_"], "rb"))
    ## ---- python newdata 3
    newdata_3 = pd.MultiIndex.from_product(
        [benthos_fixed_locs_obs_3["fYear"].unique()],
        names=["fYear"]
    ).to_frame(index=False)
    newdata_3
    ## ----end
    pickle.dump(newdata_3, open(product, "wb"))

def newdata_3b_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_3 = pickle.load(open(upstream["sampled_reefs_data_prep_3_"], "rb"))
    ## ---- python newdata 3b
    pd.set_option("display.max_columns", None)  # Show all columns
    newdata_3b = (
      benthos_fixed_locs_obs_3.groupby(["Year", "Reef"], as_index=False, observed=True).agg(
        HCC = ("HCC", "mean"),
        CYC = ("CYC", "mean"),
        DHW = ("DHW", "mean"),
        OTHER = ("OTHER", "mean"),
        Longitude = ("Longitude", "mean"),
        Latitude = ("Latitude", "mean")
      ).reset_index()
    )
    newdata_3b
    ## ----end
    pickle.dump(newdata_3b, open(product, "wb"))

    
## Temporal gap of all reefs -----------------------------------------

def sampled_reefs_data_prep_4_(product, upstream):
    benthos_fixed_locs_obs_4 = pickle.load(open(upstream["read_sampled_reefs_data_4_"], "rb"))
    ## ---- python sampled data prep 4
    benthos_fixed_locs_obs_4 = (
        benthos_fixed_locs_obs_4.assign(
            fYear=lambda df: df["Year"].astype("category"),
            Reef=lambda df: df["Reef"].astype("category"),
            Site=lambda df: (df["Reef"].astype("str") + "_" +
                             df["Site"]).astype("category"),  # Interaction of Reef and Site
            Transect=lambda df: (df["Site"].astype("str") + "_" +
                                 df["Transect"]).astype("category"),  # Interaction of Site and Transect
            cover=lambda df: df["HCC"] / 100  # Calculate cover
        )
    )
    benthos_fixed_locs_obs_4
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_4, open(product, "wb"))

def newdata_4_(product, upstream):
    benthos_fixed_locs_obs_4 = pickle.load(open(upstream["sampled_reefs_data_prep_4_"], "rb"))
    ## ---- python newdata 4
    newdata_4 = pd.MultiIndex.from_product(
        [benthos_fixed_locs_obs_4["fYear"].unique()],
        names=["fYear"]
    ).to_frame(index=False)
    newdata_4
    ## ----end
    pickle.dump(newdata_4, open(product, "wb"))

def newdata_4b_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_4 = pickle.load(open(upstream["sampled_reefs_data_prep_4_"], "rb"))
    ## ---- python newdata 4b
    pd.set_option("display.max_columns", None)  # Show all columns
    newdata_4b = (
      benthos_fixed_locs_obs_4.groupby(["Year", "Reef"], as_index=False, observed=True).agg(
        HCC = ("HCC", "mean"),
        CYC = ("CYC", "mean"),
        DHW = ("DHW", "mean"),
        OTHER = ("OTHER", "mean"),
        Longitude = ("Longitude", "mean"),
        Latitude = ("Latitude", "mean")
      ).reset_index()
    )
    newdata_4b
    ## ----end
    pickle.dump(newdata_4b, open(product, "wb"))

## Fit models ============================================

## All sampled reefs --------------------------------------------------

def missing_years_sampled_reefs_bart_data_0_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_0 = pickle.load(open(upstream["missing_years_sampled_reefs_data_prep_0_"], "rb"))
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ---- python missing years sampled data barts data 0
    benthos_fixed_locs_obs_0_prep = benthos_fixed_locs_obs_0.copy()
    for col in ["Reef", "Site", "Transect"]:
        benthos_fixed_locs_obs_0_prep[col] = benthos_fixed_locs_obs_0_prep[col].astype("category").cat.remove_unused_categories()
        
    benthos_fixed_locs_obs_0_prep.groupby(["Year"]).agg(
        {"cover": "mean"})
    sample_data = benthos_fixed_locs_obs_0_prep.copy()

    full_data = benthos_reefs_sf
    full_data[["Longitude", "Latitude"]] = full_data["geometry"].str.extract(r"c\(([^,]+),\s*([^)]+)\)")
    # Convert to numeric
    full_data["Longitude"] = pd.to_numeric(full_data["Longitude"])
    full_data["Latitude"] = pd.to_numeric(full_data["Latitude"])
    ## ----end
    data_dict = {
        "sample_data": sample_data,
        "full_data": full_data
    }
    pickle.dump(data_dict, open(product, "wb"))

def missing_years_sampled_reefs_bart_fit_0_(product, upstream):
    data_dict = pickle.load(open(upstream["missing_years_sampled_reefs_bart_data_0_"], "rb"))
    sample_data = data_dict["sample_data"]
    full_data = data_dict["full_data"]
    ## ---- python missing years sampled data barts fit 0
    X = sample_data[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = sample_data["cover"]
    with pm.Model() as model_1:
        data_X = pm.Data("data_X", X)
        data_Y = pm.Data("data_Y", Y)
        # σ = pm.HalfNormal("σ", 1)
        μ = pmb.BART("μ", data_X, Y, m=100)
        # y = pm.Normal("y", mu=μ, sigma=σ, observed=data_Y)
        # Alpha and Beta shape parameters can be set or estimated
        alpha = pm.HalfNormal("alpha", sigma=2)
        beta = pm.Deterministic("beta", alpha * (1 - μ) / μ)
        y = pm.Beta("y", alpha=alpha, beta=beta, observed=data_Y)
#     model_1 = bart_model_0(X, Y) 
    with model_1:
        idata = pm.sample(random_seed=123)
    with model_1:
        ppc1 = pm.sample_posterior_predictive(idata)

    another_X = full_data[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]].copy()
    X = another_X[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = full_data["HCC"]
    with model_1:
            pm.set_data({"data_X": another_X,
                        "data_Y": np.arange(another_X.shape[0])
                        })
            ppc2 = pm.sample_posterior_predictive(idata)
    ## ----end
    trace_dict = {
        "idata": idata,
        "ppc1": ppc1,
        "ppc2": ppc2
    }
    pickle.dump(trace_dict, open(product, "wb"))

def missing_years_sampled_reefs_bart_prep_plot_0_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["missing_years_sampled_reefs_bart_fit_0_"], "rb"))
    ppc1 = trace_dict["ppc1"]
    data_dict = pickle.load(open(upstream["missing_years_sampled_reefs_bart_data_0_"], "rb"))
    sample_data = data_dict["sample_data"]
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ---- python missing years sampled data barts prep plot 0 
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

def missing_years_sampled_reefs_bart_plot_0_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["missing_years_sampled_reefs_bart_prep_plot_0_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    benthos_reefs_temporal_summary = summary_dict["benthos_reefs_temporal_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python missing years sampled data barts plot 0 
    def multiply_by_100(x, pos):
        return f"{x * 100:.0f}"
    plt.figure(figsize=(8, 6))
    plt.plot(years, cover_mean, "w", lw=3)
    # Add true mean line
    plt.plot(
        benthos_reefs_temporal_summary["Year"],
        benthos_reefs_temporal_summary["Mean"],
        label="True mean",
        color="blue"
    )
    # Add true median line
    plt.plot(
        benthos_reefs_temporal_summary["Year"],
        benthos_reefs_temporal_summary["Median"],
        label="True median",
        color="green"
    )
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Mean"]/100,
        label="Sample mean",
        color="blue",
        linestyle="--"
    )
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Median"]/100,
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
    plt.savefig(f"{paths['fig_path']}python_pdp_pymc_bart_0b2.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

def missing_years_sampled_reefs_bart_prep_plot_0c_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["missing_years_sampled_reefs_bart_fit_0_"], "rb"))
    ppc2 = trace_dict["ppc2"]
    data_dict = pickle.load(open(upstream["missing_years_sampled_reefs_bart_data_0_"], "rb"))
    sample_data = data_dict["sample_data"]
    full_data = data_dict["full_data"]
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ---- python missing years sampled data barts prep plot 0c 
    posterior_predictive = ppc2.posterior_predictive["y"]  # replace "y_pred" with your variable name
    year_array = full_data["Year"].values  # must match the shape of predictions

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

def missing_years_sampled_reefs_bart_plot_0c_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["missing_years_sampled_reefs_bart_prep_plot_0c_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    benthos_reefs_temporal_summary = summary_dict["benthos_reefs_temporal_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python missing years sampled data barts plot 0c 
    def multiply_by_100(x, pos):
        return f"{x * 100:.0f}"
    plt.figure(figsize=(8, 6))
    plt.plot(years, cover_mean, "w", lw=3)
    # Add true mean line
    plt.plot(
        benthos_reefs_temporal_summary["Year"],
        benthos_reefs_temporal_summary["Mean"],
        label="True mean",
        color="blue"
    )
    # Add true median line
    plt.plot(
        benthos_reefs_temporal_summary["Year"],
        benthos_reefs_temporal_summary["Median"],
        label="True median",
        color="green"
    )
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Mean"]/100,
        label="Sample mean",
        color="blue",
        linestyle="--"
    )
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Median"]/100,
        label="Sample median",
        color="green",
        linestyle="--"
    )
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"95% HDI"}, hdi_prob=0.95)
    az.plot_hdi(years, cover, smooth=False,
                fill_kwargs={"label":"50% HDI", "alpha": 1}, hdi_prob=0.50)
    # Add labels and title
    plt.xlabel("Year")
    plt.ylabel("Mean")
    plt.title("Mean with HDI per Year")
    plt.legend()
    ## ----end
    plt.savefig(f"{paths['fig_path']}python_pdp_pymc_bart_0c2.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

## Temporal gap of some reefs -----------------------------------------

def sampled_reefs_bart_data_3_(product, upstream):
    ## This function prepares the data for the BART model
    ## It prepares two datasets:
    ## 1. `sample_data`: The sampled data (25 sites (etc), 12 years) that
    ##    will be used to fit the model.
    ## 2. `full_data`: The full dataset (all sites, all years) that will
    ##    be used to make predictions on the scale of the full
    ##    spatio-temporal domain.
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_3 = pickle.load(open(upstream["sampled_reefs_data_prep_3_"], "rb"))
    newdata_3b = pickle.load(open(upstream["newdata_3b_"], "rb"))
    newdata_b = pickle.load(open(upstream["missing_years_newdata_b_"], "rb"))
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ---- python sampled data barts data 3
    benthos_fixed_locs_obs_3_prep = benthos_fixed_locs_obs_3.copy()
    for col in ["Reef", "Site", "Transect"]:
        benthos_fixed_locs_obs_3_prep[col] = benthos_fixed_locs_obs_3_prep[col].astype("category").cat.remove_unused_categories()
        
    benthos_fixed_locs_obs_3_prep.groupby(["Year"]).agg(
        {"cover": "mean"})
    sample_data = benthos_fixed_locs_obs_3_prep.copy()

    full_data = benthos_reefs_sf
    full_data[["Longitude", "Latitude"]] = full_data["geometry"].str.extract(r"c\(([^,]+),\s*([^)]+)\)")
    # Convert to numeric
    full_data["Longitude"] = pd.to_numeric(full_data["Longitude"])
    full_data["Latitude"] = pd.to_numeric(full_data["Latitude"])

    newdata_3b = newdata_3b
    newdata_b = newdata_b
    ## ----end
    data_dict = {
        "sample_data": sample_data,
        "full_data": full_data,
        "newdata_3b": newdata_3b,
        "newdata_b": newdata_b
    }
    pickle.dump(data_dict, open(product, "wb"))

def sampled_reefs_bart_fit_3_(product, upstream):
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_3_"], "rb"))
    sample_data = data_dict["sample_data"]
    full_data = data_dict["full_data"]
    newdata_3b = data_dict["newdata_3b"]
    newdata_b = data_dict["newdata_b"]
    ## ---- python sampled data barts fit 3
    X = sample_data[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = sample_data["cover"]
    with pm.Model() as model_3:
        data_X = pm.Data("data_X", X)
        data_Y = pm.Data("data_Y", Y)
        μ = pmb.BART("μ", data_X, Y, m=100)
        alpha = pm.HalfNormal("alpha", sigma=2)
        beta = pm.Deterministic("beta", alpha * (1 - μ) / μ)
        y = pm.Beta("y", alpha=alpha, beta=beta, observed=data_Y)
    with model_3:
        idata = pm.sample(random_seed=123)
    with model_3:
        ppc1 = pm.sample_posterior_predictive(idata)

    # another_X = full_data[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]].copy()
    # X = another_X[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    # Y = full_data["HCC"]
    another_X = newdata_b[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]].copy()
    X = another_X[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = newdata_b["HCC"]
    with model_3:
        pm.set_data({"data_X": another_X,
                     "data_Y": np.arange(another_X.shape[0])
                     })
        ppc2 = pm.sample_posterior_predictive(idata)

    another_X = newdata_3b[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]].copy()
    X = another_X[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = newdata_3b["HCC"]
    with model_3:
        pm.set_data({"data_X": another_X,
                    "data_Y": np.arange(another_X.shape[0])
                    })
        ppc3 = pm.sample_posterior_predictive(idata)
    ## ----end
    trace_dict = {
        "idata": idata,
        "ppc1": ppc1,
        "ppc2": ppc2,
        "ppc3": ppc3
    }
    pickle.dump(trace_dict, open(product, "wb"))

def sampled_reefs_bart_prep_plot_3_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_3_"], "rb"))
    ppc3 = trace_dict["ppc3"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_3_"], "rb"))
    newdata_3b = data_dict["newdata_3b"]
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    sample_data = data_dict["sample_data"]
    ## ---- python sampled data barts prep plot 3 
    posterior_predictive = ppc3.posterior_predictive["y"]  # replace "y_pred" with your variable name
    # year_array = sample_data["Year"].values  # must match the shape of predictions
    year_array = newdata_3b["Year"].values  # must match the shape of predictions

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

def sampled_reefs_bart_plot_3_(product, upstream):
    # paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_3_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    benthos_reefs_temporal_summary = summary_dict["benthos_reefs_temporal_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python sampled data barts plot 3 
    def multiply_by_100(x, pos):
        return f"{x * 100:.0f}"
    plt.figure(figsize=(8, 6))
    plt.plot(years, cover_mean, "w", lw=3)
    # # Add true mean line
    # plt.plot(
    #     benthos_reefs_temporal_summary["Year"],
    #     benthos_reefs_temporal_summary["Mean"],
    #     label="True mean",
    #     color="blue"
    # )
    # # Add true median line
    # plt.plot(
    #     benthos_reefs_temporal_summary["Year"],
    #     benthos_reefs_temporal_summary["Median"],
    #     label="True median",
    #     color="green"
    # )
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Mean"]/100,
        label="Sample mean",
        color="blue",
        linestyle="--"
    )
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Median"]/100,
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
    plt.savefig(f"{paths['fig_path']}python_pdp_pymc_bart_3b.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_mse_3b_1_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_3_"], "rb"))
    ppc3 = trace_dict["ppc3"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_3_"], "rb"))
    newdata_3b = data_dict["newdata_3b"]
    ## ---- python sampled data pymc_barts mse 3b 1 
    posterior_predictive = ppc3.posterior_predictive["y"]
    mse = (newdata_3b["HCC"]/100 - posterior_predictive.mean(dim=["draw", "chain"]))**2
    acc = np.abs(np.exp(np.log(newdata_3b["HCC"]/100) - np.log(posterior_predictive.mean(dim=["draw", "chain"]))) - 1)
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
    df.to_csv(f"{paths['data_path']}/modelled/pymc_bart_3b_mse_1.csv", index=False)
    ## ----end
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_prep_plot_3c_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_3_"], "rb"))
    ppc2 = trace_dict["ppc2"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_3_"], "rb"))
    sample_data = data_dict["sample_data"]
    full_data = data_dict["full_data"]
    newdata_b = data_dict["newdata_b"]
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ---- python sampled data barts prep plot 3c 
    posterior_predictive = ppc2.posterior_predictive["y"]  # replace "y_pred" with your variable name
    year_array = newdata_b["Year"].values  # must match the shape of predictions

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

def sampled_reefs_bart_plot_3c_(product, upstream):
    # paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_3c_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    benthos_reefs_temporal_summary = summary_dict["benthos_reefs_temporal_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python sampled data barts plot 3c 
    def multiply_by_100(x, pos):
        return f"{x * 100:.0f}"
    plt.figure(figsize=(8, 6))
    plt.plot(years, cover_mean, "w", lw=3)
    # Add true mean line
    plt.plot(
        benthos_reefs_temporal_summary["Year"],
        benthos_reefs_temporal_summary["Mean"],
        label="True mean",
        color="blue"
    )
    # Add true median line
    plt.plot(
        benthos_reefs_temporal_summary["Year"],
        benthos_reefs_temporal_summary["Median"],
        label="True median",
        color="green"
    )
    # # Add sample mean line
    # plt.plot(
    #     sample_data_summary["Year"],
    #     sample_data_summary["Mean"]/100,
    #     label="Sample mean",
    #     color="blue",
    #     linestyle="--"
    # )
    # # Add sample mean line
    # plt.plot(
    #     sample_data_summary["Year"],
    #     sample_data_summary["Median"]/100,
    #     label="Sample median",
    #     color="green",
    #     linestyle="--"
    # )
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
    plt.savefig(f"{paths['fig_path']}python_pdp_pymc_bart_3c.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_mse_3b_2_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_3_"], "rb"))
    ppc2 = trace_dict["ppc2"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_3_"], "rb"))
    newdata_b = data_dict["newdata_b"]
    ## ---- python sampled data pymc_barts mse 3b 2 
    posterior_predictive = ppc2.posterior_predictive["y"]
    mse = (newdata_b["HCC"] - posterior_predictive.mean(dim=["draw", "chain"]))**2
    acc = np.abs(np.exp(np.log(newdata_b["HCC"]) - np.log(posterior_predictive.mean(dim=["draw", "chain"]))) - 1)
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
    df.to_csv(f"{paths['data_path']}/modelled/pymc_bart_3b_mse_2.csv", index=False)
    ## ----end
    results=1
    pickle.dump(results, open(product, "wb"))

## Temporal gap of all reefs -----------------------------------------

def sampled_reefs_bart_data_4_(product, upstream):
    ## This function prepares the data for the BART model
    ## It prepares two datasets:
    ## 1. `sample_data`: The sampled data (25 sites (etc), 12 years) that
    ##    will be used to fit the model.
    ## 2. `full_data`: The full dataset (all sites, all years) that will
    ##    be used to make predictions on the scale of the full
    ##    spatio-temporal domain.
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_4 = pickle.load(open(upstream["sampled_reefs_data_prep_4_"], "rb"))
    newdata_4b = pickle.load(open(upstream["newdata_4b_"], "rb"))
    newdata_b = pickle.load(open(upstream["missing_years_newdata_b_"], "rb"))
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ---- python sampled data barts data 4
    benthos_fixed_locs_obs_4_prep = benthos_fixed_locs_obs_4.copy()
    for col in ["Reef", "Site", "Transect"]:
        benthos_fixed_locs_obs_4_prep[col] = benthos_fixed_locs_obs_4_prep[col].astype("category").cat.remove_unused_categories()
        
    benthos_fixed_locs_obs_4_prep.groupby(["Year"]).agg(
        {"cover": "mean"})
    sample_data = benthos_fixed_locs_obs_4_prep.copy()

    full_data = benthos_reefs_sf
    full_data[["Longitude", "Latitude"]] = full_data["geometry"].str.extract(r"c\(([^,]+),\s*([^)]+)\)")
    # Convert to numeric
    full_data["Longitude"] = pd.to_numeric(full_data["Longitude"])
    full_data["Latitude"] = pd.to_numeric(full_data["Latitude"])

    newdata_4b = newdata_4b
    newdata_b = newdata_b
    ## ----end
    data_dict = {
        "sample_data": sample_data,
        "full_data": full_data,
        "newdata_4b": newdata_4b,
        "newdata_b": newdata_b
    }
    pickle.dump(data_dict, open(product, "wb"))

def sampled_reefs_bart_fit_4_(product, upstream):
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_4_"], "rb"))
    sample_data = data_dict["sample_data"]
    full_data = data_dict["full_data"]
    newdata_4b = data_dict["newdata_4b"]
    newdata_b = data_dict["newdata_b"]
    ## ---- python sampled data barts fit 4
    X = sample_data[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = sample_data["cover"]
    with pm.Model() as model_4:
        data_X = pm.Data("data_X", X)
        data_Y = pm.Data("data_Y", Y)
        μ = pmb.BART("μ", data_X, Y, m=100)
        alpha = pm.HalfNormal("alpha", sigma=2)
        beta = pm.Deterministic("beta", alpha * (1 - μ) / μ)
        y = pm.Beta("y", alpha=alpha, beta=beta, observed=data_Y)
    with model_4:
        idata = pm.sample(random_seed=123)
    with model_4:
        ppc1 = pm.sample_posterior_predictive(idata)

    another_X = newdata_b[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]].copy()
    X = another_X[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = newdata_b["HCC"]
    with model_4:
        pm.set_data({"data_X": another_X,
                     "data_Y": np.arange(another_X.shape[0])
                     })
        ppc2 = pm.sample_posterior_predictive(idata)
    
    another_X = newdata_4b[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]].copy()
    X = another_X[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = newdata_4b["HCC"]
    with model_4:
        pm.set_data({"data_X": another_X,
                    "data_Y": np.arange(another_X.shape[0])
                    })
        ppc3 = pm.sample_posterior_predictive(idata)
    ## ----end
    trace_dict = {
        "idata": idata,
        "ppc1": ppc1,
        "ppc2": ppc2,
        "ppc3": ppc3
    }
    pickle.dump(trace_dict, open(product, "wb"))

def sampled_reefs_bart_prep_plot_4_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_4_"], "rb"))
    ppc3 = trace_dict["ppc3"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_4_"], "rb"))
    newdata_4b = data_dict["newdata_4b"]
    sample_data = data_dict["sample_data"]
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ---- python sampled data barts prep plot 4 
    posterior_predictive = ppc3.posterior_predictive["y"]  # replace "y_pred" with your variable name
    year_array = newdata_4b["Year"].values  # must match the shape of predictions

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

def sampled_reefs_bart_plot_4_(product, upstream):
    # paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_4_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    benthos_reefs_temporal_summary = summary_dict["benthos_reefs_temporal_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python sampled data barts plot 4 
    def multiply_by_100(x, pos):
        return f"{x * 100:.0f}"
    plt.figure(figsize=(8, 6))
    plt.plot(years, cover_mean, "w", lw=3)
    # # Add true mean line
    # plt.plot(
    #     benthos_reefs_temporal_summary["Year"],
    #     benthos_reefs_temporal_summary["Mean"],
    #     label="True mean",
    #     color="blue"
    # )
    # # Add true median line
    # plt.plot(
    #     benthos_reefs_temporal_summary["Year"],
    #     benthos_reefs_temporal_summary["Median"],
    #     label="True median",
    #     color="green"
    # )
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Mean"]/100,
        label="Sample mean",
        color="blue",
        linestyle="--"
    )
    # Add sample mean line
    plt.plot(
        sample_data_summary["Year"],
        sample_data_summary["Median"]/100,
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
    plt.savefig(f"{paths['fig_path']}python_pdp_pymc_bart_4b.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_mse_4b_1_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_4_"], "rb"))
    ppc3 = trace_dict["ppc3"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_4_"], "rb"))
    newdata_4b = data_dict["newdata_4b"]
    ## ---- python sampled data pymc_barts mse 4b 1 
    posterior_predictive = ppc3.posterior_predictive["y"]
    mse = (newdata_4b["HCC"]/100 - posterior_predictive.mean(dim=["draw", "chain"]))**2
    acc = np.abs(np.exp(np.log(newdata_4b["HCC"]/100) - np.log(posterior_predictive.mean(dim=["draw", "chain"]))) - 1)
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
    df.to_csv(f"{paths['data_path']}/modelled/pymc_bart_4b_mse_1.csv", index=False)
    ## ----end
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_prep_plot_4c_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_4_"], "rb"))
    ppc2 = trace_dict["ppc2"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_4_"], "rb"))
    sample_data = data_dict["sample_data"]
    full_data = data_dict["full_data"]
    newdata_b = data_dict["newdata_b"]
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ---- python sampled data barts prep plot 4c 
    posterior_predictive = ppc2.posterior_predictive["y"]  # replace "y_pred" with your variable name
    year_array = newdata_b["Year"].values  # must match the shape of predictions

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

def sampled_reefs_bart_plot_4c_(product, upstream):
    # paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_4c_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    benthos_reefs_temporal_summary = summary_dict["benthos_reefs_temporal_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python sampled data barts plot 4c 
    def multiply_by_100(x, pos):
        return f"{x * 100:.0f}"
    plt.figure(figsize=(8, 6))
    plt.plot(years, cover_mean, "w", lw=3)
    # Add true mean line
    plt.plot(
        benthos_reefs_temporal_summary["Year"],
        benthos_reefs_temporal_summary["Mean"],
        label="True mean",
        color="blue"
    )
    # Add true median line
    plt.plot(
        benthos_reefs_temporal_summary["Year"],
        benthos_reefs_temporal_summary["Median"],
        label="True median",
        color="green"
    )
    # # Add sample mean line
    # plt.plot(
    #     sample_data_summary["Year"],
    #     sample_data_summary["Mean"]/100,
    #     label="Sample mean",
    #     color="blue",
    #     linestyle="--"
    # )
    # # Add sample mean line
    # plt.plot(
    #     sample_data_summary["Year"],
    #     sample_data_summary["Median"]/100,
    #     label="Sample median",
    #     color="green",
    #     linestyle="--"
    # )
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
    plt.savefig(f"{paths['fig_path']}python_pdp_pymc_bart_4c.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_mse_4b_2_(product, upstream):
    paths = pickle.load(open(upstream["missing_years_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_4_"], "rb"))
    ppc2 = trace_dict["ppc2"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_4_"], "rb"))
    newdata_b = data_dict["newdata_b"]
    ## ---- python sampled data pymc_barts mse 4b 2 
    posterior_predictive = ppc2.posterior_predictive["y"]
    mse = (newdata_b["HCC"] - posterior_predictive.mean(dim=["draw", "chain"]))**2
    acc = np.abs(np.exp(np.log(newdata_b["HCC"]) - np.log(posterior_predictive.mean(dim=["draw", "chain"]))) - 1)
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
    df.to_csv(f"{paths['data_path']}/modelled/pymc_bart_4b_mse_2.csv", index=False)
    ## ----end
    results=1
    pickle.dump(results, open(product, "wb"))
