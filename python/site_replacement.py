# + tags=["parameters"]
# declare a list tasks whose products you want to use as inputs
upstream = None
# -

import multiprocessing

## ---- python site replacement libraries
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

# script_dir = "/home/murray/Work/AIMS/Projects/GCRMN/2025/dev/python/"
# os.chdir(script_dir)

# def change_cwd_(product):
#     os.chdir(os.path.join(os.getcwd(), "python"))
#     pickle.dump("done", open(product, "wb"))

# os.chdir(os.path.join(os.getcwd(), "..", "python"))
# print(os.getcwd())

def site_replacement_global_parameters_(product):
    ## ---- python site replacement global parameters
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

def read_all_reefs_data_(product, upstream):
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    ## ---- python read all reefs data
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ----end
    pickle.dump(benthos_reefs_sf, open(product, "wb"))

def read_all_reefs_data_show_(product, upstream):
    benthos_reefs_sf = pickle.load(open(upstream["read_all_reefs_data_"], "rb"))
    ## ---- python read all reefs data show
    pd.set_option("display.max_columns", None)  # Show all columns
    benthos_reefs_sf 
    pd.reset_option("display.max_columns")
    ## ----end
    pickle.dump("done", open(product, "wb"))

def read_all_temporal_summary_(product, upstream):
    benthos_reefs_sf = pickle.load(open(upstream["read_all_reefs_data_"], "rb"))
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

def read_all_temporal_summary_plot_(product, upstream):
    benthos_reefs_temporal_summary = pickle.load(open(upstream["read_all_temporal_summary_"], "rb"))
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
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
    plt.savefig(f"{paths['fig_path']}Python_all_temporal_summary_plot.png", dpi=300, bbox_inches="tight")
    ## ----end
    plt.savefig(product, dpi=72, bbox_inches="tight")

## All sampled reefs --------------------------------------------------

def read_sampled_reefs_data_(product, upstream):
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    ## ---- python read sampled reefs data
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

def read_sampled_reefs_data_show_(product, upstream):
    benthos_fixed_locs_obs = pickle.load(open(upstream["read_sampled_reefs_data_"], "rb"))
    ## ---- python read sampled reefs data show
    pd.set_option("display.max_columns", None)  # Show all columns
    benthos_fixed_locs_obs
    pd.reset_option("display.max_columns")
    ## ----end
    pickle.dump("done", open(product, "wb"))

def sampled_simple_raw_means_(product, upstream):
    benthos_fixed_locs_obs = pickle.load(open(upstream["read_sampled_reefs_data_"], "rb"))
    ## ---- python sampled simple raw means
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

def sampled_simple_raw_means_plot_(product, upstream):
    all_sampled_sum = pickle.load(open(upstream["sampled_simple_raw_means_"], "rb"))
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    ## ---- python sampled simple raw means plot
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
    plt.savefig(f"{paths['fig_path']}Python_full_simple_raw_means_plot.png",
                dpi=72, bbox_inches="tight")
    ## ----end
    plt.savefig(product, dpi=72, bbox_inches="tight")

## Replace a reef --------------------------------------------------

def read_sampled_reefs_data_1_(product, upstream):
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    ## ---- python read sampled reefs data 1
    benthos_fixed_locs_obs_1 = pd.read_csv(f"{paths['data_path']}synthetic/benthos_fixed_locs_obs_1.csv")
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
    benthos_fixed_locs_obs_1 = benthos_fixed_locs_obs_1.merge(
        benthos_reefs_summary,
        on=["Year", "Reef"],
        how="left"
    )
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_1, open(product, "wb"))

def read_sampled_reefs_data_1_show_(product, upstream):
    benthos_fixed_locs_obs_1 = pickle.load(open(upstream["read_sampled_reefs_data_1_"], "rb"))
    ## ---- python read sampled reefs data 1 show
    pd.set_option("display.max_columns", None)  # Show all columns
    benthos_fixed_locs_obs_1
    pd.reset_option("display.max_columns")
    ## ----end
    pickle.dump("done", open(product, "wb"))

def sampled_reefs_data_1_plot_(product, upstream):
    benthos_fixed_locs_obs_1 = pickle.load(open(upstream["read_sampled_reefs_data_1_"], "rb"))
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    ## ---- python sampled reefs data 1 plot
    # Create the plot
    g = sns.FacetGrid(
        data=benthos_fixed_locs_obs_1,
        col="Reef",
        col_wrap=5,  # Adjust the number of columns in the facet grid
        height=2,
        aspect=1.6,
        sharey=True,
        sharex=True
    )
    
    # Add line plots to each facet
    g.map_dataframe(
        sns.lineplot,
        x="Year",
        y="HCC",
        hue="Site",
        style="Transect",
        estimator=None,
        units="Transect",
        legend=False
    )

    # Customize the plot
    g.set_axis_labels("Year", "Hard Coral Cover (HCC)")
    g.set_titles(col_template="{col_name}")
    g.add_legend(title="Site")
    g.fig.tight_layout()

    # Show the plot
    #plt.show()
    
    # Save the plot as a PNG file
    plt.savefig(f"{paths['fig_path']}Python_sampled_reefs_1_plot.png", dpi=72, bbox_inches="tight")
    ## ----end
    plt.savefig(product, dpi=72, bbox_inches="tight")

## Replace a reef (V2) ---------------------------------------------

def read_sampled_reefs_data_2_(product, upstream):
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    ## ---- python read sampled reefs data 2
    benthos_fixed_locs_obs_2 = pd.read_csv(f"{paths['data_path']}synthetic/benthos_fixed_locs_obs_2.csv")
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
    benthos_fixed_locs_obs_2 = benthos_fixed_locs_obs_2.merge(
        benthos_reefs_summary,
        on=["Year", "Reef"],
        how="left"
    )
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_2, open(product, "wb"))

def read_sampled_reefs_data_2_show_(product, upstream):
    benthos_fixed_locs_obs_2 = pickle.load(open(upstream["read_sampled_reefs_data_2_"], "rb"))
    ## ---- python read sampled reefs data 2 show
    pd.set_option("display.max_columns", None)  # Show all columns
    benthos_fixed_locs_obs_2
    pd.reset_option("display.max_columns")
    ## ----end
    pickle.dump("done", open(product, "wb"))

def sampled_reefs_data_2_plot_(product, upstream):
    benthos_fixed_locs_obs_2 = pickle.load(open(upstream["read_sampled_reefs_data_2_"], "rb"))
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    ## ---- python sampled reefs data 2 plot
    # Create the plot
    g = sns.FacetGrid(
        data=benthos_fixed_locs_obs_2,
        col="Reef",
        col_wrap=3,  # Adjust the number of columns in the facet grid
        height=3,
        aspect=1.6,
        sharey=True,
        sharex=True
    )
    
    # Add line plots to each facet
    g.map_dataframe(
        sns.lineplot,
        x="Year",
        y="HCC",
        hue="Site",
        style="Transect",
        estimator=None,
        units="Transect",
        legend=False
    )

    # Customize the plot
    g.set_axis_labels("Year", "Hard Coral Cover (HCC)")
    g.set_titles(col_template="{col_name}")
    g.add_legend(title="Site")
    g.fig.tight_layout()

    # Show the plot
    #plt.show()
    
    # Save the plot as a PNG file
    plt.savefig(f"{paths['fig_path']}Python_sampled_reefs_2_plot.png", dpi=72, bbox_inches="tight")
    ## ----end
    plt.savefig(product, dpi=72, bbox_inches="tight")

## Data transformations ============================================

## All sampled reefs --------------------------------------------------

def sampled_reefs_data_prep_0_(product, upstream):
    benthos_fixed_locs_obs = pickle.load(open(upstream["read_sampled_reefs_data_"], "rb"))
    ## ---- python sampled data prep 0
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
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_0, open(product, "wb"))

def newdata_0_(product, upstream):
    benthos_fixed_locs_obs_0 = pickle.load(open(upstream["sampled_reefs_data_prep_0_"], "rb"))
    ## ---- python newdata 0
    newdata_0 = pd.MultiIndex.from_product(
        [benthos_fixed_locs_obs_0["fYear"].unique()],
        names=["fYear"]
    ).to_frame(index=False)
    ## ----end
    pickle.dump(newdata_0, open(product, "wb"))

## Replace a reef --------------------------------------------------

def sampled_reefs_data_prep_1_(product, upstream):
    benthos_fixed_locs_obs_1 = pickle.load(open(upstream["read_sampled_reefs_data_1_"], "rb"))
    ## ---- python sampled data prep 1
    benthos_fixed_locs_obs_1 = (
        benthos_fixed_locs_obs_1.assign(
            fYear=lambda df: df["Year"].astype("category"),
            Reef=lambda df: df["Reef"].astype("category"),
            Site=lambda df: (df["Reef"].astype("str") + "_" +
                             df["Site"]).astype("category"),  # Interaction of Reef and Site
            Transect=lambda df: (df["Site"].astype("str") + "_" +
                                 df["Transect"]).astype("category"),  # Interaction of Site and Transect
            cover=lambda df: df["HCC"] / 100  # Calculate cover
        )
    )
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_1, open(product, "wb"))

def newdata_1_(product, upstream):
    benthos_fixed_locs_obs_1 = pickle.load(open(upstream["sampled_reefs_data_prep_1_"], "rb"))
    ## ---- python newdata 1
    newdata_1 = pd.MultiIndex.from_product(
        [benthos_fixed_locs_obs_1["fYear"].unique()],
        names=["fYear"]
    ).to_frame(index=False)
    ## ----end
    pickle.dump(newdata_1, open(product, "wb"))

## Replace a reef (V2) ---------------------------------------------
def sampled_reefs_data_prep_2_(product, upstream):
    benthos_fixed_locs_obs_2 = pickle.load(open(upstream["read_sampled_reefs_data_2_"], "rb"))
    ## ---- python sampled data prep 2
    benthos_fixed_locs_obs_2 = (
        benthos_fixed_locs_obs_2.assign(
            fYear=lambda df: df["Year"].astype("category"),
            Reef=lambda df: df["Reef"].astype("category"),
            Site=lambda df: (df["Reef"].astype("str") + "_" +
                             df["Site"]).astype("category"),  # Interaction of Reef and Site
            Transect=lambda df: (df["Site"].astype("str") + "_" +
                                 df["Transect"]).astype("category"),  # Interaction of Site and Transect
            cover=lambda df: df["HCC"] / 100  # Calculate cover
        )
    )
    ## ----end
    pickle.dump(benthos_fixed_locs_obs_2, open(product, "wb"))

def newdata_2_(product, upstream):
    benthos_fixed_locs_obs_2 = pickle.load(open(upstream["sampled_reefs_data_prep_2_"], "rb"))
    ## ---- python newdata 2
    newdata_2 = pd.MultiIndex.from_product(
        [benthos_fixed_locs_obs_2["fYear"].unique()],
        names=["fYear"]
    ).to_frame(index=False)
    ## ----end
    pickle.dump(newdata_2, open(product, "wb"))

## Data transformations ============================================

## All sampled reefs --------------------------------------------------


def sampled_reefs_bart_data_0_(product, upstream):
    ## This function prepares the data for the BART model
    ## It prepares two datasets:
    ## 1. `sample_data`: The sampled data (25 sites (etc), 12 years) that
    ##    will be used to fit the model.
    ## 2. `full_data`: The full dataset (all sites, all years) that will
    ##    be used to make predictions on the scale of the full
    ##    spatio-temporal domain.
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_0 = pickle.load(open(upstream["sampled_reefs_data_prep_0_"], "rb"))
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ---- python sampled data barts data 0
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

def sampled_reefs_bart_fit_0_(product, upstream):
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_0_"], "rb"))
    sample_data = data_dict["sample_data"]
    full_data = data_dict["full_data"]
#     benthos_reefs_sf = pickle.load(open(upstream["sampled_reefs_bart_prep_0_"], "rb"))
#     benthos_fixed_locs_obs_0_prep = pickle.load(open(upstream["sampled_reefs_bart_prep_0_"], "rb"))
    ## ---- python sampled data barts fit 0
    X = sample_data[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = sample_data["cover"]
    with pm.Model() as model_1:
        data_X = pm.Data("data_X", X)
        data_Y = pm.Data("data_Y", Y)
        σ = pm.HalfNormal("σ", 1)
        μ = pmb.BART("μ", data_X, Y, m=100)
        y = pm.Normal("y", mu=μ, sigma=σ, observed=data_Y)
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

def sampled_reefs_bart_prep_plot_0_(product, upstream):
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_0_"], "rb"))
    ppc1 = trace_dict["ppc1"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_0_"], "rb"))
    sample_data = data_dict["sample_data"]
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ---- python sampled data barts prep plot 0 
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

def sampled_reefs_bart_plot_0_(product, upstream):
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_0_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    benthos_reefs_temporal_summary = summary_dict["benthos_reefs_temporal_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python sampled data barts plot 0 
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
    plt.savefig(f"{paths['fig_path']}python_pdp_pymc_bart_0b.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_plot_0_old(product, upstream):
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_0_"], "rb"))
    ppc1 = trace_dict["ppc1"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_0_"], "rb"))
    sample_data = data_dict["sample_data"]
    
    X = sample_data[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = sample_data["cover"]
    another_X = X[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]].copy()
    ## ---- pyth
    n_chains, n_draws, n_observations = ppc1.posterior_predictive["y"].shape
    flattened_predictions = ppc1.posterior_predictive["y"].values#.reshape(-1, n_observations)

    ## Repeat the prediction grid for each draw
    repeated_grid = pd.concat([another_X] * (n_chains * n_draws), ignore_index=False)
    ## Add predictions to the grid
    repeated_grid["Posterior_Predictions"] = flattened_predictions.flatten()
    repeated_grid["Draw"] = np.tile(np.repeat(np.arange(n_draws), n_observations), n_chains)
    repeated_grid["Chain"] = np.repeat(np.arange(n_chains), n_draws * n_observations)
    aa = repeated_grid.groupby(["Year","Chain","Draw"])["Posterior_Predictions"].mean().reset_index()

    grouped = aa.groupby("Year")["Posterior_Predictions"]
    years = []
    means = []
    hdis = []

    # Loop through each group
    for year, predictions in grouped:
        years.append(year)
        means.append(predictions.mean())  # Calculate mean
        hdi = az.hdi(predictions.values, hdi_prob=0.95)  # Calculate 95% HDI
        hdis.append(hdi)
    results = pd.DataFrame({
        "Year": years,
        "Mean": means,
        "HDI_Lower": [hdi[0] for hdi in hdis],
        "HDI_Upper": [hdi[1] for hdi in hdis]
    })

    # print(results)

    df = results
    plt.figure(figsize=(8, 6))
    # Plot the mean values with error bars
    plt.errorbar(
        df["Year"], df["Mean"], 
        yerr=[df["Mean"] - df["HDI_Lower"], df["HDI_Upper"] - df["Mean"]], 
        fmt='o', capsize=4, label="Mean with HDI"
    )
    plt.fill_between(
        df["Year"], df["HDI_Lower"], df["HDI_Upper"], 
        color="blue", alpha=0.2, label="HDI (95%)"
    )

    # Add labels and title
    plt.xlabel("Year")
    plt.ylabel("Mean")
    plt.title("Mean with HDI per Year")
    plt.legend()

    plt.savefig(f"{paths['fig_path']}python_pdp_pymc_bart_0b.png", dpi=300, bbox_inches='tight')
    ## ----end
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_post_0c_(product, upstream):
    # benthos_reefs_sf = pickle.load(open(upstream["sampled_reefs_bart_prep_0_"], "rb"))
    idata = pickle.load(open(upstream["sampled_reefs_bart_fit_0_"], "rb"))
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    ## ---- python sampled data barts post 0c
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    df = benthos_reefs_sf
    df[["Longitude", "Latitude"]] = df["geometry"].str.extract(r"c\(([^,]+),\s*([^)]+)\)")
    # Convert to numeric
    df["Longitude"] = pd.to_numeric(df["Longitude"])
    df["Latitude"] = pd.to_numeric(df["Latitude"])

    benthos_reefs_summary = (
       df
        .drop(columns=["geometry"])  # Equivalent to st_drop_geometry()
        .groupby(["Year", "Reef"], as_index=False)
        .agg({
            "Latitude": "mean",
            "Longitude": "mean",
            "CYC": "mean",
            "DHW": "mean",
            "OTHER": "mean"
        })
    )
    # benthos_reefs_summary

    another_X = df[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]].copy()
    X = another_X[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = df["HCC"]
    model_1 = bart_model_0(X, Y) 
    with model_1:
            pm.set_data({"data_X": another_X,
                        "data_Y": np.arange(another_X.shape[0])
                        })
            ppc2 = pm.sample_posterior_predictive(idata)
    ## ----end
    pickle.dump(ppc2, open(product, "wb"))

def sampled_reefs_bart_plot_0c_old(product, upstream):
    # paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_0_"], "rb"))
    ppc2 = trace_dict["ppc2"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_0_"], "rb"))
    full_data = data_dict["full_data"]
    ## ---- pyth
    # benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    # df = benthos_reefs_sf
    # df[["Longitude", "Latitude"]] = df["geometry"].str.extract(r"c\(([^,]+),\s*([^)]+)\)")
    # # Convert to numeric
    # df["Longitude"] = pd.to_numeric(df["Longitude"])
    # df["Latitude"] = pd.to_numeric(df["Latitude"])
    another_X = full_data[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]].copy()

    n_chains, n_draws, n_observations = ppc2.posterior_predictive["y"].shape
    flattened_predictions = ppc2.posterior_predictive["y"].values#.reshape(-1, n_observations)

    ## Repeat the prediction grid for each draw
    repeated_grid = pd.concat([another_X] * (n_chains * n_draws), ignore_index=False)
    ## Add predictions to the grid
    repeated_grid["Posterior_Predictions"] = flattened_predictions.flatten()
    repeated_grid["Draw"] = np.tile(np.repeat(np.arange(n_draws), n_observations), n_chains)
    repeated_grid["Chain"] = np.repeat(np.arange(n_chains), n_draws * n_observations)
    aa = repeated_grid.groupby(["Year","Chain","Draw"])["Posterior_Predictions"].mean().reset_index()

    grouped = aa.groupby("Year")["Posterior_Predictions"]
    years = []
    means = []
    hdis = []

    # Loop through each group
    for year, predictions in grouped:
        years.append(year)
        means.append(predictions.mean())  # Calculate mean
        hdi = az.hdi(predictions.values, hdi_prob=0.95)  # Calculate 95% HDI
        hdis.append(hdi)
    results = pd.DataFrame({
        "Year": years,
        "Mean": means,
        "HDI_Lower": [hdi[0] for hdi in hdis],
        "HDI_Upper": [hdi[1] for hdi in hdis]
    })

    # print(results)

    df = results
    plt.figure(figsize=(8, 6))
    # Plot the mean values with error bars
    plt.errorbar(
        df["Year"], df["Mean"], 
        yerr=[df["Mean"] - df["HDI_Lower"], df["HDI_Upper"] - df["Mean"]], 
        fmt='o', capsize=4, label="Mean with HDI"
    )
    plt.fill_between(
        df["Year"], df["HDI_Lower"], df["HDI_Upper"], 
        color="blue", alpha=0.2, label="HDI (95%)"
    )

    # Add labels and title
    plt.xlabel("Year")
    plt.ylabel("Mean")
    plt.title("Mean with HDI per Year")
    plt.legend()

    plt.savefig(f"{paths['fig_path']}python_pdp_pymc_bart_0c.png", dpi=300, bbox_inches='tight')
    ## ----end
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_prep_plot_0c_(product, upstream):
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_0_"], "rb"))
    ppc2 = trace_dict["ppc2"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_0_"], "rb"))
    sample_data = data_dict["sample_data"]
    full_data = data_dict["full_data"]
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ---- python sampled data barts prep plot 0c 
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
    
def sampled_reefs_bart_plot_0c_(product, upstream):
    # paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_0c_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    benthos_reefs_temporal_summary = summary_dict["benthos_reefs_temporal_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python sampled data barts plot 0c 
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
    plt.savefig(f"{paths['fig_path']}python_pdp_pymc_bart_0c.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

## Replace a reef --------------------------------------------------

def sampled_reefs_bart_data_1_(product, upstream):
    ## This function prepares the data for the BART model
    ## It prepares two datasets:
    ## 1. `sample_data`: The sampled data (25 sites (etc), 12 years) that
    ##    will be used to fit the model.
    ## 2. `full_data`: The full dataset (all sites, all years) that will
    ##    be used to make predictions on the scale of the full
    ##    spatio-temporal domain.
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_1 = pickle.load(open(upstream["sampled_reefs_data_prep_1_"], "rb"))
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ---- python sampled data barts data 1
    benthos_fixed_locs_obs_1_prep = benthos_fixed_locs_obs_1.copy()
    for col in ["Reef", "Site", "Transect"]:
        benthos_fixed_locs_obs_1_prep[col] = benthos_fixed_locs_obs_1_prep[col].astype("category").cat.remove_unused_categories()
        
    benthos_fixed_locs_obs_1_prep.groupby(["Year"]).agg(
        {"cover": "mean"})
    sample_data = benthos_fixed_locs_obs_1_prep.copy()

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

def sampled_reefs_bart_fit_1_(product, upstream):
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_1_"], "rb"))
    sample_data = data_dict["sample_data"]
    full_data = data_dict["full_data"]
    ## ---- python sampled data barts fit 1
    X = sample_data[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = sample_data["cover"]
    with pm.Model() as model_1:
        data_X = pm.Data("data_X", X)
        data_Y = pm.Data("data_Y", Y)
        σ = pm.HalfNormal("σ", 1)
        μ = pmb.BART("μ", data_X, Y, m=100)
        y = pm.Normal("y", mu=μ, sigma=σ, observed=data_Y)
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

def sampled_reefs_bart_prep_plot_1_(product, upstream):
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_1_"], "rb"))
    ppc1 = trace_dict["ppc1"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_1_"], "rb"))
    sample_data = data_dict["sample_data"]
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ---- python sampled data barts prep plot 1 
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

def sampled_reefs_bart_plot_1_(product, upstream):
    # paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_1_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    benthos_reefs_temporal_summary = summary_dict["benthos_reefs_temporal_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python sampled data barts plot 1 
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
    plt.savefig(f"{paths['fig_path']}python_pdp_pymc_bart_1b.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_prep_plot_1c_(product, upstream):
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_1_"], "rb"))
    ppc2 = trace_dict["ppc2"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_1_"], "rb"))
    sample_data = data_dict["sample_data"]
    full_data = data_dict["full_data"]
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ---- python sampled data barts prep plot 1c 
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

def sampled_reefs_bart_plot_1c_(product, upstream):
    # paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_1c_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    benthos_reefs_temporal_summary = summary_dict["benthos_reefs_temporal_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python sampled data barts plot 1c 
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
    plt.savefig(f"{paths['fig_path']}python_pdp_pymc_bart_1c.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

## Replace a reef (V2) ---------------------------------------------
def sampled_reefs_bart_data_2_(product, upstream):
    ## This function prepares the data for the BART model
    ## It prepares two datasets:
    ## 1. `sample_data`: The sampled data (25 sites (etc), 12 years) that
    ##    will be used to fit the model.
    ## 2. `full_data`: The full dataset (all sites, all years) that will
    ##    be used to make predictions on the scale of the full
    ##    spatio-temporal domain.
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    benthos_fixed_locs_obs_2 = pickle.load(open(upstream["sampled_reefs_data_prep_2_"], "rb"))
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ---- python sampled data barts data 2
    benthos_fixed_locs_obs_2_prep = benthos_fixed_locs_obs_2.copy()
    for col in ["Reef", "Site", "Transect"]:
        benthos_fixed_locs_obs_2_prep[col] = benthos_fixed_locs_obs_2_prep[col].astype("category").cat.remove_unused_categories()
        
    benthos_fixed_locs_obs_2_prep.groupby(["Year"]).agg(
        {"cover": "mean"})
    sample_data = benthos_fixed_locs_obs_2_prep.copy()

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

def sampled_reefs_bart_fit_2_(product, upstream):
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_2_"], "rb"))
    sample_data = data_dict["sample_data"]
    full_data = data_dict["full_data"]
    ## ---- python sampled data barts fit 2
    X = sample_data[["Year", "Latitude", "Longitude", "CYC", "DHW", "OTHER"]]
    Y = sample_data["cover"]
    with pm.Model() as model_1:
        data_X = pm.Data("data_X", X)
        data_Y = pm.Data("data_Y", Y)
        σ = pm.HalfNormal("σ", 1)
        μ = pmb.BART("μ", data_X, Y, m=100)
        y = pm.Normal("y", mu=μ, sigma=σ, observed=data_Y)
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

def sampled_reefs_bart_prep_plot_2_(product, upstream):
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_2_"], "rb"))
    ppc1 = trace_dict["ppc1"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_2_"], "rb"))
    sample_data = data_dict["sample_data"]
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ---- python sampled data barts prep plot 2 
    ## Summarise the modelled trend
    posterior_predictive = ppc1.posterior_predictive["y"] 
    year_array = sample_data["Year"].values  
    # Add 'Year' as a coordinate to the xarray
    posterior_predictive = posterior_predictive.assign_coords({"y_dim_0": year_array})
    # Group by Year and compute mean per chain/draw
    avg_preds = posterior_predictive.groupby("y_dim_0").mean("y_dim_0")
    # Rename dimension for clarity
    avg_preds = avg_preds.rename({"y_dim_0": "Year"})
    summary = az.summary(avg_preds, hdi_prob=0.95)
    summary["Year"] = avg_preds.coords["Year"].values
    print(summary)

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
        "modelled_trend": summary,
        "sample_data_summary": sample_data_summary,
        "benthos_reefs_temporal_summary": benthos_reefs_temporal_summary,
        "years": years,
        "cover_mean": cover_mean,
        "cover": cover
    }
    pickle.dump(summary_dict, open(product, "wb"))
    
def sampled_reefs_bart_plot_2_(product, upstream):
    # paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    # trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_2_"], "rb"))
    # ppc1 = trace_dict["ppc1"]
    # data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_2_"], "rb"))
    # sample_data = data_dict["sample_data"]
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_2_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    benthos_reefs_temporal_summary = summary_dict["benthos_reefs_temporal_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    ## ---- python sampled data barts plot 2 
    # Define a function to multiply the y-axis values by 100
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
    plt.savefig(f"{paths['fig_path']}python_pdp_pymc_bart_2b.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))

def sampled_reefs_bart_prep_plot_2c_(product, upstream):
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_2_"], "rb"))
    ppc2 = trace_dict["ppc2"]
    data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_2_"], "rb"))
    sample_data = data_dict["sample_data"]
    full_data = data_dict["full_data"]
    benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
    ## ---- python sampled data barts prep plot 2c 
    ## Summarise the modelled trend
    posterior_predictive = ppc2.posterior_predictive["y"] 
    year_array = full_data["Year"].values  
    # Add 'Year' as a coordinate to the xarray
    posterior_predictive = posterior_predictive.assign_coords({"y_dim_0": year_array})
    # Group by Year and compute mean per chain/draw
    avg_preds = posterior_predictive.groupby("y_dim_0").mean("y_dim_0")
    # Rename dimension for clarity
    avg_preds = avg_preds.rename({"y_dim_0": "Year"})
    summary = az.summary(avg_preds, hdi_prob=0.95)
    summary["Year"] = avg_preds.coords["Year"].values
    print(summary)

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
        "modelled_trend": summary,
        "sample_data_summary": sample_data_summary,
        "benthos_reefs_temporal_summary": benthos_reefs_temporal_summary,
        "years": years,
        "cover_mean": cover_mean,
        "cover": cover
    }
    pickle.dump(summary_dict, open(product, "wb"))

def sampled_reefs_bart_plot_2c_(product, upstream):
    # paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    # ppc2 = pickle.load(open(upstream["sampled_reefs_bart_post_0c_"], "rb"))
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    summary_dict = pickle.load(open(upstream["sampled_reefs_bart_prep_plot_2c_"], "rb"))
    modelled_trend = summary_dict["modelled_trend"]
    sample_data_summary = summary_dict["sample_data_summary"]
    benthos_reefs_temporal_summary = summary_dict["benthos_reefs_temporal_summary"]
    years = summary_dict["years"]
    cover_mean = summary_dict["cover_mean"]
    cover = summary_dict["cover"]
    # trace_dict = pickle.load(open(upstream["sampled_reefs_bart_fit_2_"], "rb"))
    # ppc2 = trace_dict["ppc2"]
    # data_dict = pickle.load(open(upstream["sampled_reefs_bart_data_2_"], "rb"))
    # full_data = data_dict["full_data"]
    ## ---- python sampled data barts plot 2c 
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
    # posterior_predictive = ppc2.posterior_predictive["y"]  # replace "y_pred" with your variable name
    # year_array = full_data["Year"].values  # must match the shape of predictions

    # # Add 'Year' as a coordinate to the xarray
    # posterior_predictive = posterior_predictive.assign_coords({"y_dim_0": year_array})
    # # Group by Year and compute mean per chain/draw
    # avg_preds = posterior_predictive.groupby("y_dim_0").mean("y_dim_0")
    # # Rename dimension for clarity
    # avg_preds = avg_preds.rename({"y_dim_0": "Year"})
    # summary = az.summary(avg_preds, hdi_prob=0.95)
    # summary["Year"] = avg_preds.coords["Year"].values
    # print(summary)

    # plt.figure(figsize=(8, 6))
    # # Plot the mean values with error bars
    # plt.errorbar(
    #     summary["Year"], summary["mean"], 
    #     yerr=[summary["mean"] - summary["hdi_2.5%"], summary["hdi_97.5%"] - summary["mean"]], 
    #     fmt='o', capsize=4, label="Mean with HDI"
    # )
    # plt.fill_between(
    #     summary["Year"], summary["hdi_2.5%"], summary["hdi_97.5%"], 
    #     color="blue", alpha=0.2, label="HDI (95%)"
    # )

    # # Add labels and title
    # plt.xlabel("Year")
    # plt.ylabel("Mean")
    # plt.title("Mean with HDI per Year")
    # plt.legend()
    ## ----end
    plt.savefig(f"{paths['fig_path']}python_pdp_pymc_bart_2c.png", dpi=300, bbox_inches='tight')
    results=1
    pickle.dump(results, open(product, "wb"))
