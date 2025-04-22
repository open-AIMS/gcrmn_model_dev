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

# ## ---- python read all reefs data
# print("Working directory set to:", os.getcwd())
# # benthos_reefs_sf_np = pd.read_csv("../data/synthetic/benthos_reefs_sf.csv")
# benthos_reefs_sf_np = pd.read_csv(f"{paths['data_path']}synthetic/benthos_reefs_sf.csv")
# pd.set_option("display.max_columns", None)  # Show all columns
# benthos_reefs_sf_np 
# pd.reset_option("display.max_columns")
# ## ----end

# ## ---- python all reefs temporal summary
# benthos_reefs_temporal_summary = (
#     benthos_reefs_sf_np.groupby("Year").agg(
#         Mean=("HCC", "mean"),
#         Median=("HCC", "median"),
#         SD=("HCC", "std"),
#         Lower=("HCC", lambda x: x.quantile(0.025)),
#         Upper=("HCC", lambda x: x.quantile(0.975))
#     ).reset_index()
# )
# ## ----end

# ## ---- python all reefs temporal summary plot
# # Create the plot
# plt.figure(figsize=(10, 6))

# # Add ribbon (shaded area)
# plt.fill_between(
#     benthos_reefs_temporal_summary["Year"],
#     benthos_reefs_temporal_summary["Lower"],
#     benthos_reefs_temporal_summary["Upper"],
#     alpha=0.2,
#     label="Confidence Interval"
# )

# # Add mean line
# plt.plot(
#     benthos_reefs_temporal_summary["Year"],
#     benthos_reefs_temporal_summary["Mean"],
#     label="Mean",
#     color="blue"
# )

# # Add median line
# plt.plot(
#     benthos_reefs_temporal_summary["Year"],
#     benthos_reefs_temporal_summary["Median"],
#     label="Median",
#     color="orange"
# )

# # Add labels and theme
# plt.xlabel("Year")
# plt.ylabel("Value")
# plt.title("Temporal Summary")
# plt.legend()
# plt.grid(True)
# # plt.show()

# # Save the plot as a PNG file
# plt.savefig("../data/synthetic/benthos_reefs_sf_temporal_summary.png", dpi=300, bbox_inches="tight")
# ## ----end

# ## ---- python sampled read all reefs data
# benthos_fixed_locs_obs_0 = pd.read_csv(f"{paths['data_path']}synthetic/benthos_fixed_locs_obs_0.csv")
# pd.set_option("display.max_columns", None)  # Show all columns
# benthos_fixed_locs_obs_0 
# pd.reset_option("display.max_columns")
# ## ----end

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
