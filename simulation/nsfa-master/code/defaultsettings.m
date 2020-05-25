function settings=defaultsettings()

settings.a=1; % noise variance hyperparameters
settings.b=1;
settings.a0=1; % hyperparameters on b
settings.b0=1;
settings.d=1; % factor variance hyperparameters
settings.c0=1;
settings.d0=1;
settings.e=1; % IBP alpha hyperparamters
settings.f=1;
settings.alpha=1; 
settings.Beta=1;
settings.iterations=1000;
settings.verbose=1;
settings.addNoiseToMV=0; 
settings.store_samples=0;
settings.sampleNoiseVariance=0;
settings.lambdae=10; 
settings.sampleMixingMatrixHyperparameters=0; 
settings.lambdag=1; 

settings.sampleBeta=0;
settings.sampleNewG=1;
settings.sampleZG=1; 
settings.sampleX=1;
settings.ibpUniform=0;

settings.output='nsfa_results.mat';

% Interesting settings
settings.sampleAlpha=1; % whether to sample the IBP concentration parameter
settings.isotropic=0; % enforce isotropic noise?
settings.sparse=1; % sparse prior on mixing coefficients?
settings.np=1; % nonparametric-use IBP? Only valid if sparse
settings.perfactor=1; % Per factor covariances? With c=1 corresponds to ARD prior
settings.learnscale=1; % Hierarchical prior on the factor covariances? Only valid if perfactor is true
settings.c=1; % 1 corresponds to an ARD prior. 2 corresponds to not. 
settings.K=0; % Whether we should initialise with a specific number of features
settings.betaFlag=false; % use two-parameter ibp?
settings.shareNoise=false; % share power across noise dimensions? 
settings.fokoue=0; % use Fokoue's method