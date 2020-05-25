function sample=init_nsfa(settings)
% Make an initial sample for the MCMC sampler. 

    sample.alpha=1;
    if ~settings.sampleAlpha
        sample.alpha=settings.alpha;
    end
    sample.Beta=1;
    D=settings.D;
    if settings.sampleZG
        if settings.np
            if settings.K
                sample.Z=rand(settings.D,settings.K) < 0.5;
            else
                % Initialise sample.Z to a random sample from the IBP
                sample.Z = ibpgen(settings.D,sample.alpha);
            end
            % K is the number of non-zero rows in sample.Z
            [D K] = size(sample.Z);
        else
            K=settings.K;
            D=settings.D;
            sample.Z=rand(settings.D,settings.K) < 0.5;
        end
        
    if settings.fokoue
        sample.lambdag=ones(D,K);
    else
        if settings.sampleMixingMatrixHyperparameters
            sample.lambdag=ones(1,K)*gamrnd(settings.c,settings.d);
        else
            sample.lambdag=ones(1,K)*settings.lambdag; 
        end
    end

        sample.G = normrnd(0,1,[settings.D K])*diag(sample.lambdag.^-.5);
        sample.Gold = sample.G;

        sample.G = sample.G.*sample.Z;
    else
        K=settings.K;
        sample.Z=settings.Z; 
        sample.G = settings.G; 
    end
    
    
    
    sample.b=settings.b;
    % Initialise alpha,sigmae,sigmag, and sample.Beta
    if settings.sampleNoiseVariance
        sample.lambdae=ones(1,D)*gamrnd(settings.a,sample.b);
    else
        sample.lambdae=ones(1,D)*settings.lambdae;
    end

    % Randomly initialise X and G from their priors
    if settings.sampleX
        sample.X = normrnd(0,1,[K settings.N]);
    else
        sample.X = settings.X;
    end
    
    sample.d=settings.d;
end