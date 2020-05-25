addpath('utils');

M1=1e4;
M2=1e3;
thinout=10;
burnin=1e2;

settings=defaultsettings();
settings.N=2;
settings.D=2;
settings.verbose=0;
settings.iterations=1; 

% marginal conditional simulator
g1=[];
for m=1:M1
    if (mod(m,1000)==0)
        fprintf(1,'marginal conditional iteration %d\n',m);
    end
    g1(m,:) = marginal_conditional(settings, settings.N, settings.D);
end 

% successive conditional simulator
param_set=init_nsfa(settings);
g2=[];
for m=1:(burnin+M2*thinout)
    if (mod(m,1000)==0)
        fprintf(1,'successive conditional iteration %d, K=%d\n'...
            ,m,size(param_set.Z,2));
    end
    
    [temp param_set]=successive_conditional(param_set, settings, settings.N, settings.D);
    if (m>burnin && mod(m,thinout)==0)
        g2=[g2;temp'];
    end
end

 
g_bar1=mean(g1)
g_bar2=mean(g2)
sigma_g1=std(g1)
sigma_g2=std(g2)
size(g1)
size(g2)
test_stat=(g_bar1-g_bar2)./  ...
    (sigma_g1.^2/M1+sigma_g2.^2/M2).^.5; 
p_values=2*normcdf(-abs(test_stat));
p_values