function [prior] = buildMoGPrior(gam,nu,w,theseVals)

% Reconstruct estimated prior from a set of variables

if iscolumn(w)
    w   = w';
    nu  = nu';
    gam = gam';
end

% Define support in log units
if isempty(theseVals)
    xLin = [0.5 1 2 4 8 12];
    x    = getLogXform(xLin,0.3);
else
    x = theseVals;
end

% Number of components in prior
numComp = size(w,2);

% Assume a second dimension in input parameters defines number of different
% trials
numRuns = size(w,1);

% Build prior as a weighted sum of Gaussians
prior = zeros(numRuns,numel(x));
for jj = 1:numRuns
    for ii = 1:numComp
        prior(jj,:) = prior(jj,:) + w(jj,ii)*normpdf(x,nu(jj,ii),gam(jj,ii));
    end
end

% Normalize to obtain a probability density
diffVel = diff(x);
auc     = sum(repmat(diffVel,[numRuns 1]).*(prior(:,2:end)+prior(:,1:end-1))/2,2);
prior   = prior./repmat(auc,[1 numel(x)]);

end