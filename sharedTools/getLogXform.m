function [vTilde] = getLogXform(v,v0)

% Get log transformed speed
%
% Usage: [vTilde] = getLogXform(v,v0)

vTilde = log(1 + v/v0);

end