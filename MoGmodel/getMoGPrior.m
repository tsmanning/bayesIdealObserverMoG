function [prior] = getMoGPrior(gam,nu,w,thisSupport)

% Ensure weights sum to 1
w = w/sum(w);

numComp = numel(w);

% Make a prior from parameters [w,nu,gam]
prior = zeros(1,numel(thisSupport));

for ii = 1:numComp
    prior = prior + w(ii)*normpdf(thisSupport,nu(ii),gam(ii));
end

end