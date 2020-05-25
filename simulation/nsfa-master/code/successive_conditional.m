function [g param_set]=successive_conditional(param_set, settings, N, D)
% function g=successive_conditional(param_set, stat_set, N, D)
% Sample data ~ model|params, params ~ MCMC q|params,data,
% g = test_functions

Y = sample_data(param_set, N, D);
param_set = nsfa(Y,ones(D,N),param_set,settings);
g = test_functions(param_set, Y);