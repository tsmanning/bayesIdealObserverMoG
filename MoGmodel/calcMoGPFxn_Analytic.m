function [pFxn] = calcMoGPFxn_Analytic(priorMus,priorSigs,priorWs,mu1,sig1,mu2,sig2)

% Given a 2D likelihood function, what's the probability of responding s2
% instad of s1? Analytical approximation via Eqn 37.

% Collect modified weights and means and the shrinkage factors that define
% posteriors using approximation of measurements with expected value
[alphas1,muTildes1,wTildes1] = getMogPostPars(priorWs,priorSigs,priorMus,sig1,mu1);
[alphas2,muTildes2,wTildes2] = getMogPostPars(priorWs,priorSigs,priorMus,sig2,mu2);

% Calculate numerator components
postMean1 = wTildes1'*(alphas1*mu1 + muTildes1);
postMean2 = wTildes2'*(alphas2*mu2 + muTildes2);

% Calculate denominator
denom = sqrt( (sig1^2)*(wTildes1'*alphas1)^2 + (sig2^2)*(wTildes2'*alphas2)^2 );

% Calculate psychometric function value
pFxn = normcdf( (postMean2 - postMean1)/denom );

end