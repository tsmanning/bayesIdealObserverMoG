function [pNuHat,pGamHat,pWHat,pSigNseHat,priHat,nll] = fitEstimData_analytical()

% Fit Bayesian ideal observer model with fixed Gaussian noise from observer
% estimate data, consisting of pairs of stimuli 'x' and observer estimates
% 'x_hat'. Uses approximate analytical solution to p(xHat|x) defined in
% Eqn. 56 of the main text.
%
% Inputs:
% -------
%           x [N x 1]  - column vector of presented stimuli
%       x_hat [N x 1]  - column vector of observer estimates
%       prs0 [Nb+1 x 1]- initial parameters: [signse0; bwts0] (OPTIONAL)
%
% Output:
% -------
%    signse - estimate of stdev of encoding noise
%    prihat - estimated prior on x grid
%      bwts - estimated weights of prior in basis
%    Mlihat - matrix of inferred likelihood 
%  Mposthat - matrix of inferred posterior


% ----------------------------
% Extract sizes & initialize 
% ----------------------------

% Initialize parameters, if necessary
if nargin < 5
    pNu0      = zeros(nB,1);              % initial values of component centers
    pGam0     = 2*ones(nB,1);                  % initial values of component widths
    pW0       = normpdf(1:nB,(nB+1)/2,nB/4)';
    pW0       = pW0./sum(pW0);               % initial value of MoG weights
    pSigNse0  = log([1]);                  % initial estimate of noise stdev
    prs0      = [pNu0;pGam0;pW0;pSigNse0];
end

% ----------------------------
% Numerically optimize negative log-likelihood
% ----------------------------

% bounds for parameters
LB  = [-eps*ones(nB,1); -10*zeros(nB,1); zeros(nB,1);   -10*ones(1,1)]; % lower bound
UB  = [eps*ones(nB,1);  10*ones(nB,1);  ones(nB,1)+eps; 10*ones(1,1)]; % upper bound
Aeq = [zeros(1,2*nB) ones(1,nB) 0]; % equality constraint (prior weights sum to 1)
beq = 1;                              % equality constraint

% loss function
lossfun = @(prs)(neglogli_BaysObsModel(prs,stim,r,nB));

% optimization options
opts = optimset('display', 'iter','algorithm','interior-point'); 
% opts = optimset('display', 'iter','algorithm','sqp'); 

% perform optimization 
prshat = fmincon(lossfun,prs0,[],[],Aeq,beq,LB,UB,[],opts);

% ----------------------------
% Extract fitted parameters 
% ----------------------------

pNuHat  = prshat(1:nB);
pGamHat = prshat(nB + 1:2*nB);
pWHat   = prshat(2*nB + 1:3*nB);

pSigNseHat = prshat(3*nB + 1);

if nargout > 5
    
    % log-likelihood at optimum if requested
    nll = lossfun(prshat);
    
end

priHat = buildMoGPrior(pGamHat,pNuHat,pWHat,xgrid);


end

% ===================================================================
% LOSS FUNCTION: negative log-likelihood 
% ===================================================================
function [nll] = neglogli_BaysObsModel(prs,stim,r,nB)

% Computes negative log-likelihood of data (for optimization)

% Extract parameters
pNu  = prs(1:nB);
pGam = prs(nB + 1:2*nB);
pW   = prs(2*nB + 1:3*nB);

sig  = prs(3*nB + 1);

mu   = stim;

% Get probability of responses
p = calcMoGEstimDist_Analytic(pNu,pGam,pW,mu,sig,r);

% Calculate negative log-likelihood
nll = -log(p);

end