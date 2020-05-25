function fit(tau, seed)

fn=strcat('X_tau',num2str(tau),'_seed', num2str(seed));
addpath('utils');
load(strcat('../../input/',fn,'.mat'))
Y = struct2array(X);
settings=defaultsettings();
[settings.D,settings.N]=size(Y); 
settings.iterations=100;
settings.K = 5;
mvmask=binornd(1,1-0.1,settings.D,settings.N);
initialsample=init_nsfa(settings);
[finalsample,resultstable]=nsfa(Y,mvmask,initialsample,settings);

L = finalsample.G;
F = finalsample.X;
save(strcat('results/',fn,'.mat'), 'L', 'F');


end
