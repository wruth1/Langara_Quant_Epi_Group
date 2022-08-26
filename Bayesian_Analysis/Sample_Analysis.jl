
using Distributions, Random
using Plots

# Consider a very simple model: We start with an initial number of cases, and each day the number changes by a factor of theta.
# This is an applied math type model, where case counts are best thought of as the mean of some stochastic process.
# We will then generate observed data from this model (integer-valued).
# Finally, we will perform a data analysis as though the simulated observations were gathered in the wild.

Random.seed!(1)

theta_0 = 2   # Parameter of epidemic process
            # Think R0
X_1 = 1      # Initial number of cases

t_max = 5   # Duration of trajectory


# --------------------------- Generate mean process -------------------------- #
X = Vector(undef, t_max)
X[1] = X_1
for t in 2:t_max
    X[t] = theta_0 * X[t-1]
end


# ----------------------- Generate observation process ----------------------- #
data = [rand(Poisson(X[i])) for i in 1:t_max]



# ---------------------------------------------------------------------------- #
#                                 Data Analysis                                #
# ---------------------------------------------------------------------------- #

"""
Evaluate the prior for theta at the specified value.
"""
function prior(theta)
    prior_dist = Chisq(1)
    prior_val = pdf(prior_dist, theta)

    return prior_val
end

"""
Compute the posterior at theta for the given observed data.
"""
function posterior(theta, X_obs)
    log_lik = log_lik_from_theta(theta, X_obs)
    lik = exp(log_lik)

    prior_val = prior(theta)

    return lik*prior_val
end


# -------------- Compute the posterior at a grid of theta values ------------- #
theta_seq = 1.5:0.01:2.5
posterior_seq = posterior.(theta_seq, Ref(X_obs))

# ---------------------------- Plot the posterior ---------------------------- #
plot(theta_seq, posterior_seq)











# ---------------------------------------------------------------------------- #
#                               Helper Functions                               #
# ---------------------------------------------------------------------------- #
"""
Construct the mean number of cases for the given parameter values.
"""
function get_mean_trajectory(theta, X_1, t_max)
    X = [X_1 * theta^(t-1) for t in 1:t_max]
    return(X)
end

"""
Compute the log likelihood given observed data and a vector of means
"""
function log_lik_from_means(X_means, X_obs)
    all_dists = Poisson.(X_means)

    all_log_liks = [logpdf(all_dists[i], X_obs[i]) for i in 1:t_max]

    log_lik = sum(all_log_liks)
    return log_lik
end

"""
Compute the log likelihood given observed data and a value for theta.
"""
function log_lik_from_theta(theta, X_obs)
    X_means = get_mean_trajectory(theta, X_1, t_max)
    
    log_lik = log_lik_from_means(X_means, X_obs)
    return(log_lik)
end