function [v] = getLinXform(vTilde,v0)

% Transform velocity back to linear domain from log domain
%
% Usage: [vTilde] = getLinXform(v,v0)

v = v0*(exp(vTilde) - 1);

end