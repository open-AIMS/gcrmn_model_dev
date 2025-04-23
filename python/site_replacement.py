## ---- python site replacement libraries
import pandas as pd
import numpy as np
import xarray as xr
import pymc as pm
import arviz as az
import preliz as pz
import pymc_bart as pmb
import matplotlib.pyplot as plt
import os
## ----end

# script_dir = "/home/murray/Work/AIMS/Projects/GCRMN/2025/dev/python/"
# os.chdir(script_dir)

os.chdir(os.path.join(os.getcwd(), "python"))
# os.chdir(os.path.join(os.getcwd(), "..", "python"))
# print(os.getcwd())

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

## print("Working directory set to:", os.getcwd())
# benthos_reefs_sf_np = pd.read_csv("../data/synthetic/benthos_reefs_sf.csv")

## ---- python read all reefs data
benthos_reefs_sf = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
## ----end

## ---- python read all reefs data show
pd.set_option("display.max_columns", None)  # Show all columns
benthos_reefs_sf 
pd.reset_option("display.max_columns")
## ----end

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

## All sampled reefs --------------------------------------------------
## ---- python read sampled reefs data
benthos_fixed_locs_obs = pd.read_csv(f"{paths['data_path']}synthetic/benthos_fixed_locs_obs.csv")
## ----end
## ---- python read sampled reefs data show
pd.set_option("display.max_columns", None)  # Show all columns
benthos_fixed_locs_obs
pd.reset_option("display.max_columns")
## ----end

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

## ---- python sampled simple raw means plot
# Create the plot
plt.figure(figsize=(10, 6))

# Add the ribbon (confidence interval)
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
plt.savefig(f"{paths['fig_path']}Python_full_simple_raw_means_plot.png", dpi=72, bbox_inches="tight")
## ----end


# ## Data transformations

# ## All sampled reefs

# ## ---- python sampled data prep 0
# benthos_fixed_locs_obs_0 = (
#     benthos_fixed_locs_obs_0.assign(
#         fYear=lambda df: df["Year"].astype("category"),
#         Reef=lambda df: df["Reef"].astype("category"),
#         Site=lambda df: (df["Reef"].astype("str") + "_" +
#                          df["Site"]).astype("category"),  # Interaction of Reef and Site
#         Transect=lambda df: (df["Site"].astype("str") + "_" +
#                              df["Transect"]).astype("category"),  # Interaction of Site and Transect
#         cover=lambda df: df["HCC"] / 100  # Calculate cover
#     )
# )
# ## ----end

# ## ---- python newdata 0
# newdata_0 = pd.MultiIndex.from_product(
#     [benthos_fixed_locs_obs_0["fYear"].unique()],
#     names=["fYear"]
# ).to_frame(index=False)
# ## ----end

# ## ---- python newdata2 0
# # Step 1: Expand unique combinations of fYear and Transect
# newdata2_0 = pd.MultiIndex.from_product(
#     [benthos_fixed_locs_obs_0["fYear"].unique(),
#      benthos_fixed_locs_obs_0["Transect"].unique()],
#     names=["fYear", "Transect"]
# ).to_frame(index=False)

# # Step 2: Select distinct rows for Latitude, Longitude, Reef, Site, and Transect
# distinct = (
#     benthos_fixed_locs_obs_0[["Latitude", "Longitude", "Reef", "Site", "Transect"]]
#     .drop_duplicates()
# )

# # Step 3: Perform a left join
# newdata2_0 = newdata2_0.merge(distinct, on="Transect",
#                               how="left", validate="many_to_many")

# ## ----end

# # with pm.Model() as test_model:
# #     # Priors
# #     a = pm.Normal("a", mu=0, sigma=10)
# #     b = pm.Normal("b", mu=0, sigma=10, shape=benthos_fixed_locs_obs_0["fYear"].nunique())
# #     σ = pm.HalfNormal("σ", sigma=10)

# #     # Expected value of outcome
# #     μ = a + b[benthos_fixed_locs_obs_0["fYear"].cat.codes]

# #     # Likelihood (sampling distribution) of observations
# #     y_obs = pm.Normal("y_obs", mu=μ, sigma=σ, observed=benthos_fixed_locs_obs_0["cover"])

# #     idata = pm.sample()
    
# #     # # Posterior predictive checks
# #     # y_pred = pm.sample_posterior_predictive(test_model)

# # Define the logit function
# def qlogis(x):
#     return np.log(x / (1 - x))

# summary = (
#     benthos_fixed_locs_obs_0.groupby("fYear").agg(
#         mean_abundance=("cover", lambda x: qlogis(x.mean())),
#         median_abundance=("cover", lambda x: qlogis(x.median())),
#         sd_abundance=("cover", lambda x: x.pipe(qlogis).std())
#         # mad_abundance=("cover", lambda x: x.pipe(qlogis).mad()),
#     )
#     .reset_index()
# )

    
# with pm.Model() as test_model:
#     # Priors
#     a = pm.Normal("a", mu=-0.34, sigma=0.6)
#     b = pm.Normal("b", mu=0, sigma=0.4, shape=benthos_fixed_locs_obs_0["fYear"].nunique())
#     kappa = pm.Gamma("kappa", alpha=0.01, beta=0.01)  # Dispersion parameter

#     # Expected value of outcome
#     μ = pm.math.sigmoid(a + b[benthos_fixed_locs_obs_0["fYear"].cat.codes])

#     # Convert μ and kappa to alpha and beta
#     alpha = μ * kappa
#     beta = (1 - μ) * kappa
    
#     # Likelihood (sampling distribution) of observations
#     y_obs = pm.Beta("y_obs", alpha = alpha, beta = beta, observed=benthos_fixed_locs_obs_0["cover"])

#     idata = pm.sample(draws=2000,
#                        tune=1000,
#                        chains = 3,
#                        cores = 3,
#                        thin = 5,
#                        target_accept=0.95,
#                        return_inferencedata=True,
#                        random_seed=123)

# with pm.Model() as model_0:
#     μ_ = pmb.BART("μ_", X=x_data, Y=np.log(y_data), m=20)
#     μ = pm.Deterministic("μ", pm.math.exp(μ_))
#     y_pred = pm.Poisson("y_pred", mu=μ, observed=y_data)
#     idata_coal = pm.sample(random_seed=RANDOM_SEED)
