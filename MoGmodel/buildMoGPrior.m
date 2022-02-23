function [prior] = buildMoGPrior(sig,mu,w,theseVels)

% Reconstruct estimated prior from a set of variables

%%% might not want to hard code the upper range of the velocity support
%%% here
% xLin = 0:0.01:24;
if isempty(theseVels)
    xLin = [0.5 1 2 4 8 12];
    x = getLogXform(xLin,0.3);
else
    x = theseVels;
end

gauF = @(x,mu,sig) (1./(sig*sqrt(2*pi))).*exp(-0.5*((x-mu)./sig).^2);

numRuns = size(w,1);
numComp = size(w,2);

prior = zeros(numRuns,numel(x));
for jj = 1:numRuns
    for ii = 1:numComp
        prior(jj,:) = prior(jj,:) + w(jj,ii)*gauF(x,mu(jj,ii),sig(jj,ii));
    end
end

diffVel = diff(x);
auc = sum(repmat(diffVel,[numRuns 1]).*(prior(:,2:end)+prior(:,1:end-1))/2,2);

prior = prior./repmat(auc,[1 numel(x)]);

end