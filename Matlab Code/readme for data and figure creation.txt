%{
XX stores all the data.  Each row corresponds to a combination of
parameters / initial conditions.

- 1st cell is vector |BP|+|IC|, starting parameters
    Parameter order:
     1 beta = bioParameters(1); %TB infectivity
     2 p  = bioParameters(2); % ~probability someone in E goes straight into E; pi in Guo-Wu
     3 w = bioParameters(3); % period of time new infectee considered E rather than L
     4 v = bioParameters(4); % rate people in L develop TB
     5 a = bioParameters(5); % TB death rate for people in T
     6 d = bioParameters(6); % recovery rate; delta in Guo-Wu
     7 n = bioParameters(7); % natural removal rate, i.e. death rate
     8 u = bioParameters(8); % percentage people in R develop active TB 
     9 qE = q1 = bioParameters(9); % percent of immigrants into E
     10 qL = q2 = bioParameters(10); % percent of immigrants into L
     11 qR = q3 = bioParameters(11); % percent of immigrants into R
     12 X0, 13 E0, 14 L0, 15 T0, 16 R0 % initial conditions
     17 power_relapse, the power to raise error_incidence proportion by
- 2nd cell is solution vector in R9, optimal 
	x= [qE qL qR sigma, X0 E0 L0 T0 R0]
- 3rd cell is output of fmincon, {exitflag, output, lambda, grad, hessian};
- 4th cell is matrix, population XELTR vs time (using optimized params)
- 5th cell is TB incidence (using optimized params)
- 6th cell is TB prevalence
- 7th cell is error of TB incidence using optimized params
%}


