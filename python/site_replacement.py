# + tags=["parameters"]
# declare a list tasks whose products you want to use as inputs
upstream = None
# -

## ---- python site replacement libraries
import pandas as pd
import numpy as np
import xarray as xr
import pymc as pm
import arviz as az
import preliz as pz
import pymc_bart as pmb
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pickle
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
    ## ---- python read sampled reefs data show
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
    ## ---- python read sampled reefs data 1
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
    ## ---- python read sampled reefs data show
    pd.set_option("display.max_columns", None)  # Show all columns
    benthos_fixed_locs_obs_2
    pd.reset_option("display.max_columns")
    ## ----end
    pickle.dump("done", open(product, "wb"))

def sampled_reefs_data_2_plot_(product, upstream):
    benthos_fixed_locs_obs_2 = pickle.load(open(upstream["read_sampled_reefs_data_2_"], "rb"))
    paths = pickle.load(open(upstream["site_replacement_global_parameters_"], "rb"))
    ## ---- python sampled reefs data 1 plot
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
