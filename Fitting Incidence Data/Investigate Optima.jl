
using MAT
using CSV
using DataFrames
using Plots

# -------------------------- Import and Extract Data ------------------------- #
data_raw = matread("Fitting Incidence Data/Data/XX3.mat")
data = data_raw["XX"]

println(keys(data_raw))

all_par_inits = data[:,1]
all_par_hats = data[:,2]

all_incidence_hats = data[:,5]




# ------------------------- Import observed incidence ------------------------ #
incidence_obs_raw = CSV.read("Fitting Incidence Data/Data/Stats Can Incidence.csv", DataFrame)

#* Observed incidence has 3 years. Estimated incidence only has 12, and I can't find a "year" variable in XX3.mat to ensure they line-up. I'm going to just remove the earliest year from the observed incidence. Alternatively, I guess I could try both and see which is better. That sounds like a task for future William.

incidence_obs = incidence_obs_raw[2:end,2]



# ------------------- Compute optimized error for all runs ------------------- #

function obj_fun(inc_hat, inc_obs)
    return sum(((inc_hat .- inc_obs)./inc_obs).^2)
end

all_errs = [obj_fun(incidence_hat, incidence_obs) for incidence_hat in all_incidence_hats]



# ----------------------- Histogram of optimized errors ---------------------- #



