function [muF,sigF] = fitSingleGauss(support,prior)

% Fit a single Gaussian to an arbitary shaped prior

% Make sure support and prior are both row vectors
if ~isrow(support)
    support = support';
end
if ~isrow(prior)
    prior = prior';
end

lossFun = @(p) sum( (cumsum(prior) - normcdf(support,p(1),p(2))).^2 );

ub   = [ inf (support(end)-support(1))/2];
lb   = [-inf eps];
init = [0 1];

opts = optimset('display','off','tolx',1e-13,...
                'maxfunevals',1e4,'largescale', 'off');
            
pHat = fmincon(lossFun,init,[],[],[],[],lb,ub,[],opts);

muF  = pHat(1);
sigF = pHat(2);

end