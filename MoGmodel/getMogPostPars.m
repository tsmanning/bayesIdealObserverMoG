function [alphas,muTildes,wTildes,xHat] = getMogPostPars(w,gam,mus,sig,x)

% For a given mixture of Gaussians prior (w,gam,mus) and likelihood (sig,x), 
% returns parameters for posterior and BLS estimate

% Output mat dim 1: component n of MoG 
% Output mat dim 2: stimulus m

% Make sure prior pars are column vecs, likelihood pars are rows
if ~iscolumn(w)
    w = w';
end
if ~iscolumn(gam)
    gam = gam';
end
if ~iscolumn(mus)
    mus = mus';
end

if ~isrow(sig)
    sig = sig';
end
if ~isrow(x)
    x = x';
end

% Standard Normal
phi        = @(x) (1/sqrt(2*pi))*exp(-0.5*x.^2);

% Mixture component shrinkage factors
alphas     = (gam.^2) ./ (gam.^2 + sig.^2);

% Mixture component means
muTildes   = ( (sig.^2) ./ (sig.^2 + gam.^2) ).*mus;

% Adjusted mixture component weights (note this is x-dependent!)
wTildes    = (w./sqrt(gam.^2 + sig.^2)).* ...
             phi((x-mus) ./ sqrt(gam.^2 + sig.^2));
wTildes    = wTildes./sum(wTildes,1);

% Get BLS estimate
xHat       = sum(wTildes.*(alphas.*x + muTildes));

end