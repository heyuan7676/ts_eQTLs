function g = marginal_conditional(settings, N, D)
% function g = marginal_conditional(stat_set, N, D)
% Sample param ~ prior, data ~ model|param, 
% g = test_functions(params, data) 

    param_set = init_nsfa(settings);
    Y=sample_data(param_set, N, D);
    g = test_functions(param_set, Y); 