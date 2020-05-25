function y=logmvnpdf(x,mu,prec)
% Columns of x are vectors
% Should do this with cholesky really...
y=(-.5*size(x,1)*log(2*pi)+.5*log(det(prec)))*size(x,2);
for i=1:size(x,2)
    y=y-.5*(x(:,i)-mu)'*prec*(x(:,i)-mu);
end