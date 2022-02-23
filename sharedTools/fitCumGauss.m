function [pFxn,PSE] = fitCumGauss(x,p,varargin)

% Fit weibull function to psychometric data

if numel(varargin) == 0
    wts = ones(1,numel(p));
else
    wts = varargin{1};
end

lb = [0 0];
ub = [100 100];

opts = optimset('Algorithm', 'interior-point', 'Display', 'off', ...
    'MaxFunEvals', 5000, 'MaxIter', 500, 'GradObj', 'off');

% [mu sig]
parVec0 = [median(x) median(x)-min(x)];

lossFxn = @(prs) sum( wts.*(normcdf(x,prs(1),prs(2)) - p).^2 );

pars  = fmincon(lossFxn,parVec0,[],[],[],[],lb,ub,[],opts);

pFxn = normcdf(x,pars(1),pars(2));

PSE = pars(1);

end