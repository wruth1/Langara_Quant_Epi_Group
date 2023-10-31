
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
    # return sum(((inc_hat .- inc_obs)).^2)               # Unnormalized
    return sum(((inc_hat .- inc_obs)./inc_obs).^2)    # Normalized
end

all_errs = [obj_fun(incidence_hat, incidence_obs) for incidence_hat in all_incidence_hats]






# ----------------------- Histogram of optimized errors ---------------------- #

histogram(all_errs, legend=nothing, xlabel="Error", ylabel="Frequency", title="Histogram of Optimized Errors", size = (1000, 800))


# Extract configurations with large error
large_err_check = all_errs .> 0.1
large_err_ind = findall(large_err_check)
large_err_pars = all_par_inits[large_err_ind,:]

# Count frequency of each par value among large errors
all_freq_dists = []
for i in eachindex(large_err_pars[1])
    this_par_list = getindex.(large_err_pars, i)
    freq_dist = countmap(this_par_list)
    push!(all_freq_dists, freq_dist)
end



# -------------- Extract a few configurations with small errors -------------- #
n_best = 5

# Find n_best best errors
best_err_inds = sortperm(all_errs)[1:n_best]
best_errs = all_errs[best_err_inds]
best_pars = all_par_inits[best_err_inds,:]
best_par_hats = all_par_hats[best_err_inds,:]





# ------------------ Histograms of optimal parameter values ------------------ #
all_q1_hats = getindex.(all_par_hats, 1)
all_q2_hats = getindex.(all_par_hats, 2)
all_E0_hats = getindex.(all_par_hats, 3)
all_L0_hats = getindex.(all_par_hats, 4)

using LaTeXStrings
plot_q1 = histogram(all_q1_hats, legend=nothing, title=L"$q_1$");
vline!(plot_q1, getindex.(best_par_hats, 1));
plot_q2 = histogram(all_q2_hats, legend=nothing, title=L"$q_2$");
vline!(plot_q2, getindex.(best_par_hats, 2));
plot_E0 = histogram(all_E0_hats, legend=nothing, title=L"$\mathbf{E}_0$");
vline!(plot_E0, getindex.(best_par_hats, 3));
plot_L0 = histogram(all_L0_hats, legend=nothing, title=L"$\mathbf{L}_0$");
vline!(plot_L0, getindex.(best_par_hats, 4));
par_est_hists = plot(plot_q1, plot_q2, plot_E0, plot_L0, layout=(2,2), size = (1000, 800))
savefig(par_est_hists, "Fitting Incidence Data/Plots/Parameter_Estimate_Histograms.pdf")