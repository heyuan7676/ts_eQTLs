function [samp,results]=nsfa(Y,mvmask,samp,settings,results)
% Infinite Sparse Factor Analysis
% Inputs:
%   Y - Dimensions x observations data matrix
%   mvmask - 1 for elements of Y considered to be observed
%            0 for elements that should be considered to be missing
%   samp   - Initial parameters. These can be set randomly using initisFA.m
%   results - If restarting the simulation this should be the old results.
%            Otherwise it should be []
%   settings - Model/algorithm settings. Default settings can be found in
%            defaultsettings.m
%
% Outputs:
%   samp - last sample from the sampler
%   results - table of results, each row is an iteration
%
%  The columns of results are:
% 1. Cumulative CPU time
% 2. log joint
% 3. Alpha
% 4. Beta
% 5. K (number of latent factors)
% 6. Mean square predictive error on missing values
% 7. Predictive log likelihood of missing values
% 8. Sparsity of Z
% 9. Sparsity of T

if nargin==4
    results=[];
end

[D N]=size(Y);
[D K] = size(samp.G);

any_mv=any(any(mvmask==0));
% store the "original" Y for calculating predictive performance
if any_mv
    Yoriginal=Y;
    Ymv=Y(mvmask==0);
end
predError=0;
logpredlikelihood=0;

% Used to monitor acceptance ratio for MH steps
proposals=0;
acceptances=0;

totalObservations=sum(sum(mvmask));

if ~isempty(results)
    basetime=results(end,1);
else
    basetime=0;
end

if settings.store_samples
    samples=cell(settings.store_samples,1);
end

tic

M=[];
XX=[];

% samp.Gibbs sampling loop
for r=1:settings.iterations
    % Prediction
    Ypred=(samp.Z.*samp.G)*samp.X;

    % Calculate log predictive likelihood over all missing values
    if any_mv
        calculatePredictivePerformance();
    end

    if any_mv
        if settings.addNoiseToMV
            Ypred=Ypred+mvnrnd(zeros(N,D),diag(samp.lambdae.^-1))';
        end
        Y(mvmask==0)=Ypred(mvmask==0);
    end

    % Calculate error, likelihood, and posterior prob terms
    [E,loglikelihood,logposterior]=calculatePosterior();

    Z_sparsity=mean(mean(samp.Z));
    % Record some statistics about this iteration...
    results=[results;basetime+toc,logposterior,samp.alpha,samp.Beta,K,predError, ...
        logpredlikelihood,mean(mean(samp.Z))];

    % ... and output them to the screen
    if settings.verbose
        if proposals>0
            acceptance_ratio=acceptances/proposals;
        else
            acceptance_ratio=0;
        end
        fprintf(1,'%d/%d: ll:%-3.3f lp:%-3.3f b=%-3.3f alpha=%-3.3f K=%d a/r=%-3.3f lpl=%-3.3f pe=%-3.3f Zs=%-3.3f\n' ...
            ,r,settings.iterations,loglikelihood,logposterior,samp.b,samp.alpha,K,acceptance_ratio, ...
            logpredlikelihood,predError,Z_sparsity);
    end

    probZ=[];

    if settings.sampleZG
        % Sample Z and G
        if settings.sparse && ~settings.fokoue
            sampleZG();
        else
            sampleG();
        end
    end
    
    if K>0 && settings.sampleX
        sampleX();
    end

    if settings.sampleNoiseVariance
        sampleNoiseHyperparameters();
    end

    if settings.sampleMixingMatrixHyperparameters
        sampleMixingMatrixHyperparameters();
    end
    
    if settings.np && settings.sampleAlpha
        sampleIBPHyperparameterAlpha()
    end
    if settings.sampleBeta
       sampleIBPHyperparameterBeta()
    end

    if settings.store_samples>0
        for s=1:(settings.store_samples-1)
            samples{s}=samples{s+1};
        end
        samples{settings.store_samples}=samp;
    else
        samples=samp;
    end

    save(settings.output,'results','samples','settings');
end

    function calculatePredictivePerformance()
        predError=mean((Ypred(mvmask==0)-Ymv).^2);
        logpredlikelihood=0;
        for d=1:D
            Yo=Yoriginal(d,:);
            Yp=Ypred(d,:);
            dmv=length(Yp(mvmask(d,:)==0));
             logpredlikelihood=logpredlikelihood+ ...
             logmvnpdf(Yo(mvmask(d,:)==0)',Yp(mvmask(d,:)==0)',samp.lambdae(d)*eye(dmv));
        end
    end

    function sampleZG()
        XX=sum(samp.X.^2,2);
        E = Y - samp.G*samp.X; 

        if settings.np
            priorExp=samp.alpha*samp.Beta/(samp.Beta+D-1);
            poissdraws=poissrnd(priorExp,D,1);
        end

        for d = 1:D
            for k = 1:K
                m = sum(samp.Z(:,k))-samp.Z(d,k);
                if m>0 || ~settings.np
                    sampleOneElement(d,k);
                end
            end
            E(d,:)=Y(d,:)-samp.G(d,:)*samp.X; % not strictly necessary, maybe good to refresh?
            if settings.np
                sampleK(d,poissdraws(d));
            end
        end

        % Sample z(k,t)
        function sampleOneElement(d,k)
            if samp.Z(d,k)
                ek=E(d,:)+samp.G(d,k)*samp.X(k,:);
            else
                ek=E(d,:);
            end
            lambda=samp.lambdae(d)*XX(k)+samp.lambdag(k);
            mu=samp.lambdae(d)*samp.X(k,:)*ek'/lambda;
            logrp=log(m/(samp.Beta+D-1-m));
            if settings.ibpUniform
                logrp=0;
            end
            logrl=0.5*(lambda*mu^2-log(lambda)+log(samp.lambdag(k)));
            logrprop=logrp+logrl;
            probzis1 = 1/(1+exp(-logrprop));
            assert(~isnan(probzis1)); 
            samp.Z(d,k) = (rand()<probzis1);
            if samp.Z(d,k)
                samp.G(d,k)=mu+randn()*lambda^-0.5;
            else
                samp.G(d,k)=0;
                samp.Gold(d,k)=0;
            end
            E(d,:)=ek-samp.G(d,k)*samp.X(k,:);
        end
    end

    % sample new features
    function sampleK(d,Knew)

        proposals=proposals+1;
        if settings.sampleMixingMatrixHyperparameters
            if settings.perfactor|(K==0) 
                lambdag=gamrnd(settings.c,samp.d,[1 Knew]);
            else
                lambdag=ones(1,Knew)*samp.lambdag(1);
            end
        else
            lambdag=ones(1,Knew)*settings.lambdag;
        end
        
        % Calculate Lambda and mu
        m = sum(samp.Z,1) - samp.Z(d,:);
        singletons=m==0; 
        current_kappa=sum(singletons);
        Gtemp = samp.G(d,:); 
        Gtemp(singletons)=0; 
        Ed = Y(d,:) - Gtemp * samp.X; 
        current_g=samp.G(d,singletons)'; 
        prec=samp.lambdae(d)*(current_g*current_g')+eye(current_kappa);
        M=samp.lambdae(d)*(prec \ current_g)*Ed;
        
        logrprop =  N/2*log(det(prec)) - .5*trace(M'*prec*M);
        
        g = randn(Knew,1).*(lambdag.^-.5)';
        prec=samp.lambdae(d)*(g*g')+eye(Knew);

        M=samp.lambdae(d)*(prec\g)*Ed;
        
        logrprop=logrprop - N/2*log(det(prec))+.5*trace(M'*prec*M);
        if rand()<exp(logrprop)
            
            E(d,:) = E(d,:)+samp.G(d,singletons)*samp.X(singletons,:); 
            
            samp.Z(d,singletons)=0;          
            deleteRedundantFeatures(); 
            
            acceptances=acceptances+1;
            if (Knew>0)
                % accept proposal
                samp.Z = [samp.Z zeros(D,Knew)];
                % Put ones in the column t of the new rows
                samp.Z(d,K+1:K+Knew)=ones(1,Knew);
                Xprime=M+chol(prec) \ randn(Knew,N); 
                samp.X = [samp.X;Xprime];
                XX=[XX;sum(Xprime.^2,2)];
                samp.G = [samp.G zeros(D,Knew)];
                samp.Gold = [samp.Gold zeros(D,Knew)];
                samp.G(d,K+1:K+Knew)=g';
                samp.Gold(d,K+1:K+Knew)=g';
                samp.lambdag=[samp.lambdag lambdag];
                K = K+Knew;
            
                E(d,:)=E(d,:)-g'*Xprime;
                if settings.sampleNewG
                    for k=(K-Knew+1):K
                        ek=E(d,:)+samp.G(d,k)*samp.X(k,:);
                        XdotH=samp.X(k,:);
                        lambda=samp.lambdae(d)*XX(k)+samp.lambdag(k);
                        mu=samp.lambdae(d)*XdotH*ek'/lambda;
                        samp.G(d,k)=normrnd(mu,lambda^(-0.5));
                    end
                end
            end
        end
        E=Y-samp.G*samp.X; 
    end

    function deleteRedundantFeatures()
        % Delete features which are not active at any data point
        singleton_set = sum( samp.Z , 1 ) == 0;
        % Delete these features
        samp.Z(:,singleton_set)=[];
        samp.X(singleton_set,:)=[];
        XX(singleton_set)=[];
        samp.G(:,singleton_set)=[];
        samp.Gold(:,singleton_set)=[];
        samp.lambdag(singleton_set)=[];
        K=size(samp.Z,2);
    end

    % sample the latent factors X
    function sampleX()
        lambdaG=bsxfun(@times,samp.G,samp.lambdae');
        prec=samp.G'*lambdaG+eye(K);
        cholPrec = chol(prec);
        mus=cholPrec \ (cholPrec' \ (lambdaG'*Y));
        samp.X=mus+cholPrec \ randn(K,N);
    end

    % sample the loading matrix in the non-sparse case
    function sampleG()
        XX=samp.X*samp.X';
        for d=1:D
            if settings.fokoue
                prec=samp.lambdae(d)*XX+diag(samp.lambdag(d,:));
            else
                prec=samp.lambdae(d)*XX+diag(samp.lambdag);
            end
            cholPrec=chol(prec); 
            samp.G(d,:)=cholPrec \ ((cholPrec' \ (samp.lambdae(d)*samp.X*Y(d,:)')) + randn(K,1)); 
        end
    end

    function sampleNoiseHyperparameters()
        E = Y - samp.G*samp.X;
        E = E .* mvmask;
        % Sample noise level
        if settings.isotropic
            samp.lambdae=ones(1,D)*gamrndi(settings.a+totalObservations/2,samp.b+.5*trace(E'*E));
        else
            for d=1:D
                samp.lambdae(d)=gamrndi(settings.a+sum(mvmask(d,:))/2,samp.b+.5*E(d,:)*E(d,:)');
            end
            if settings.shareNoise
                samp.b=gamrndi(settings.a0+D*settings.a,settings.b0+sum(samp.lambdae));
            end
        end
    end

    function sampleMixingMatrixHyperparameters()
        % Sample mixing matrix scale
        if settings.fokoue
            for d=1:D
                for k=1:K
                    samp.lambdag(d,k)=gamrndi(settings.c+1/2,settings.d+1/2*samp.G(d,k)^2);
                end
            end
        else
            if settings.perfactor
                for k=1:K
                    samp.lambdag(k)=gamrndi(settings.c+sum(samp.Z(:,k))/2,samp.d+.5*samp.G(:,k)'*samp.G(:,k));
                end
                if settings.learnscale
                    samp.d=gamrndi(settings.c0+settings.c*K,settings.d0+sum(samp.lambdag));
                end
            else
                samp.lambdag=ones(1,K)*gamrndi(settings.c+sum(sum(samp.Z,1))/2,samp.d+.5*trace(samp.G'*samp.G));
            end
        end
    end

    function sampleIBPHyperparameterAlpha()
        % Sample alpha
        samp.alpha=gamrndi(settings.e+K,settings.f+HD(samp.Beta));
    end

    function sampleIBPHyperparameterBeta()
        % Sample beta if using two parameter model using a simple MH update
        mk=sum(samp.Z,2);
        samp.BetaPrime=gamrnd(2,1);
        propa=prod(beta(mk(1:K),N-mk(1:K)+samp.BetaPrime)./beta(mk(1:K),N-mk(1:K)+samp.Beta));
        propb=(samp.BetaPrime/samp.Beta)^K*exp(-samp.alpha*(HD(samp.BetaPrime)-HD(samp.Beta)));
        prop=propa*propb;
        if (rand()<prop)
            samp.Beta=samp.BetaPrime;
        end
    end

    function answer=HD(x)
        answer=sum(x./(x+(1:D)-1));
    end

    function [E,loglikelihood,logposterior]=calculatePosterior
        E = Y - samp.G*samp.X;
        E = E.*mvmask;

        term(1)=1/2*sum(sum(mvmask,2)'.*log(.5*samp.lambdae/pi))-sum(samp.lambdae'.*sum(E.^2,2))/2;
        %P(lambdae|a,b)
        if settings.isotropic
            term(2)=log(gampdf(samp.lambdae(1),settings.a,samp.b^-1));
        else
            term(2)=sum(log(gampdf(samp.lambdae,settings.a,samp.b^-1)));
        end
        %P(X)
        term(3)=-K*N/2*log(2*pi)-sum(sum(samp.X.^2))/2;
        term(4)=0;
        if settings.fokoue
            for d=1:D
                term(4)=term(4)+sum(log(gampdf(samp.lambdag(d,:),settings.c,settings.d^-1)));
            end
        else
            if settings.perfactor
                % This should still work if K==0
                term(4)=sum(log(gampdf(samp.lambdag,settings.c,samp.d^-1)));%P(sigma_g|c,d)
                if settings.learnscale
                    term(4)=term(4)+log(gampdf(samp.d,settings.c0,settings.d0^-1));
                end
            else
                if K==0
                    term(4)=0;
                else
                    term(4)=log(gampdf(samp.lambdag(1),settings.c,samp.d^-1));%P(sigma_g|c,d)
                end
            end
        end
        term(5)=0;
        if K>0
            for d=1:D
                if settings.fokoue
                    term(5)=term(5)+log(mvnpdf(samp.G(d,:),0,diag(samp.lambdag(d,:).^-1)));
                else
                    term(5)=term(5)+log(mvnpdf(samp.G(d,:),0,diag(samp.lambdag.^-1)));
                end
            end
        end
        mk=sum(samp.Z);
        if settings.np
            term(6)=K*log(samp.alpha*samp.Beta)-samp.alpha*HD(samp.Beta)+sum(betaln(mk(1:K),D-mk(1:K)+samp.Beta));%TODO missing one term
        elseif settings.sparse
            term(6)=sum(betaln(mk+samp.alpha/K,D-mk+1)-betaln(samp.alpha/K,1));
        else
            term(6)=0;
        end
        term(7)=log(gampdf(samp.alpha,settings.e,settings.f^-1));
        term(8)=log(gampdf(samp.Beta,1,2));
        loglikelihood=term(1);
        logposterior=sum(term);

    end

end
