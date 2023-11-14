%% Big Data Structure Experiment_Output
%{
Experiment_Output stores all the data output.  Each row corresponds to a combination of
parameters / initial conditions.

- 1st column is vector |BP|+|IC|, starting parameters in order [beta p w v a d n q1 q2, TP0, E0, L0, T0, R0]
- 2nd column is vector in R5, optimal [q1 q2 E0 L0 R0]
- 3rd column is hessian (output of fmincon)
- 4th column is matrix, population XELTR vs time (using optimized params)
- 5th column is TB incidence (using optimized params)
- 6th column is error of TB incidence using optimized params
%}

bioParameters = [beta p w v a d n q1 q2]
IC = [TP0, E0, L0, T0, R0] % TP0 is not X0, X0 = TP0 - the rest

