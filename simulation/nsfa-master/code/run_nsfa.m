addpath('utils');
load example.mat
Ycentered=Y-repmat(mean(Y,2),1,size(Y,2));
settings=defaultsettings();
[settings.D,settings.N]=size(Ycentered); 
settings.iterations=100;
mvmask=binornd(1,1-0.1,settings.D,settings.N);
initialsample=init_nsfa(settings);
[finalsample,resultstable]=nsfa(Ycentered,mvmask,initialsample,settings);
plot(resultstable(:,1),resultstable(:,2),'k+-'); 
xlabel('cpu time');
ylabel('log joint probability');