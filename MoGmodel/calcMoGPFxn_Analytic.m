function [pFxn] = calcMoGPFxn_Analytic(priorMus,priorSigs,priorWs,mu1,sig1,mu2,sig2)

% Given a 2D likelihood function, what's the probability of responding s2
% instad of s1? Analytical approximation via Eqn 51.

%% Correct dimensions

if ~isrow(mu1)
    mu1 = mu1';
end

if ~isrow(mu2)
    mu2 = mu2';
end

if ~isrow(sig1)
    sig1 = sig1';
end

if ~isrow(sig2)
    sig2 = sig2';
end


%% Get pars
% Collect modified weights and means and the shrinkage factors that define
% posteriors using approximation of measurements with expected value
[alphas1,muTildes1,wTildes1] = getMogPostPars(priorWs,priorSigs,priorMus,sig1,mu1);
[alphas2,muTildes2,wTildes2] = getMogPostPars(priorWs,priorSigs,priorMus,sig2,mu2);

% Output mat dim 1: component n of MoG 
% Output mat dim 2: stimulus m


%% Calculate numerator components
postMean1 = sum(wTildes1.*(alphas1.*mu1 + muTildes1),1);
postMean2 = sum(wTildes2.*(alphas2.*mu2 + muTildes2),1);

% Output dim 1: 1
% Output dim 2: stimulus m;


%% Calculate denominator
denom = sqrt( (sig1.^2).*sum(wTildes1.*alphas1,1).^2 + (sig2.^2).*sum(wTildes2.*alphas2,1).^2 );

% Output dim 1: 1
% Output dim 2: stimulus pair m


%% Calculate psychometric function value(s)
pFxn = normcdf( (postMean2 - postMean1)./denom );

% Output dim 1: 1
% Output dim 2: stimulus pair m


end