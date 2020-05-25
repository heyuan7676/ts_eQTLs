function g = test_functions(param_set, X)
% function g = test_functions(param_set, X)
% Calculate test function values from param_set and X

% g=X(:);
% cross_terms=X(:)*X(:)';
% for i=1:(prod(size(X))-1)
%     g=[g;cross_terms((i+1):end,i)];
% end
% g=[g;size(param_set.Z,2)];


g=size(param_set.Z,2);
g=[g;sum(sum(param_set.Z))]; 
g=[g;sum(param_set.Z,2)]; 

% g=[param_set.X(:);sum(var(param_set.X))]; 