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