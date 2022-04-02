function [pNuHat,pGamHat,pWHat,pSigNse1Hat,pSigNse2Hat,priHat,nll] = fit2AFCData_numerical(stim,r,nB,xgrid,prs0)

% Fit Bayesian ideal observer model with fixed Gaussian noise from observer
% 2AFC data, consisting of triplets of stimuli [x1,x2] and observer choices
% 'x2>x1: yes/no'. Numerically calculates p("yes"|x1,x2) by summing
% probability of joint likelihood above decision boundary.
%
% Inputs:
% -------
%          stim [N x 2] - array of presented stimuli
%             r [N x 1] - column vector of observer choices
%                    nB - scalar defining number of Gaussians in mixture
%        xgrid [Nx x 1] - grid of x evenly-spaced values on which prior is defined
%       prs0 [Nb+1 x 1] - initial parameters: [signse0; bwts0] (OPTIONAL)
%
% Output:
% -------
%    pNuHat          - estimated prior on x grid
%    pGamHat         - estimated prior on x grid
%    pWHat           - estimated prior on x grid
%    pSigNse(1/2)Hat - estimates of stdev of encoding noise
%    priHat          - estimated prior on x grid
%    nll             - negative log-likelihood of data, given parameters (OPTIONAL)


% ----------------------------
% Extract sizes & initialize 
% ----------------------------

% Initialize parameters, if necessary
if nargin < 5
    pNu0      = zeros(nB,1);              % initial values of component centers
    pGam0     = 2*ones(nB,1);                  % initial values of component widths
    pW0       = normpdf(1:nB,(nB+1)/2,nB/4)';
    pW0       = pW0./sum(pW0);               % initial value of MoG weights
    pSigNse0  = log([1; 1]);                  % initial estimate of noise stdev
    prs0      = [pNu0;pGam0;pW0;pSigNse0];
end


% ----------------------------
% Numerically optimize negative log-likelihood
% ----------------------------

% bounds for parameters
LB  = [-eps*ones(nB,1); -10*zeros(nB,1); zeros(nB,1);   -10*ones(2,1)]; % lower bound
UB  = [eps*ones(nB,1);  10*ones(nB,1);  ones(nB,1)+eps; 10*ones(2,1)]; % upper bound
Aeq = [zeros(1,2*nB) ones(1,nB) 0 0]; % equality constraint (prior weights sum to 1)
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
pGamHat = exp(prshat(nB + 1:2*nB));
pWHat   = prshat(2*nB + 1:3*nB);

pSigNse1Hat = exp(prshat(3*nB + 1));
pSigNse2Hat = exp(prshat(3*nB + 2));

if nargout > 6
    
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
pGam = exp(prs(nB + 1:2*nB));
pW   = prs(2*nB + 1:3*nB);

sig1 = exp(prs(3*nB + 1));
sig2 = exp(prs(3*nB + 2));

mu1  = stim(:,1);
mu2  = stim(:,2);

% Generate supports
dx      = 100;

nTrials  = numel(mu1);

% Calculate psychometric function values
for ii = 1:nTrials
    suppLB1 = mu1(ii) - 4*sig1;
    suppUB1 = mu1(ii) + 4*sig1;
    suppLB2 = mu2(ii) - 4*sig2;
    suppUB2 = mu2(ii) + 4*sig2;
    
    sup1 = linspace(suppLB1,suppUB1,dx);
    sup2 = linspace(suppLB2,suppUB2,dx);
    
    p(ii) = calcMoGPFxn_Numeric(sup1,sup2,pNu,pGam,pW,mu1(ii),sig1,mu2(ii),sig2,0);
end

% Compute negative log-likelihood of data
nll = -r'*log(p') - (1 - r')*log(1-p');

end