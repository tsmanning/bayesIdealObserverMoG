function [p] = calcMoGEstimDist_Analytic(priorMus,priorSigs,priorWs,mu,sig,r)

% Given a set of Bayesian observer parameters and a set of stimuli, what is
% the probability of a given estimate for each stimulus? Computes
% approximation of the true estimate distribution defined in Eqn. 56 in
% main text.

%% Get pars
% Collect modified weights and means and the shrinkage factors that define
% posteriors using approximation of measurements with expected value
[alphas,muTildes,wTildes] = getMogPostPars(priorWs,priorSigs,priorMus,sig,mu);

% Output mat dim 1: component n of MoG 
% Output mat dim 2: stimulus m


%% Calculate probabilities

x       = mu;
estMus  = sum(wTildes.*(alphas.*x + muTildes),1);
estSigs = sig.^2.*sum(wTildes.*alphas,1)^.2;

p = normpdf(r,estMus,estSigs);

% Output mat dim 1: 1
% Output mat dim 2: stimulus m

end