function [JSD] = getJSDiv(p,q)

% Calculate the Jensen-Shannon divergence between two probability distributions
%
% Usage: [JSD] = getJSDiv(p,q)
%
% Where p and q are either horizontal vectors of the same size 
% OR 
% matrices representing multiple different distribution comparisons, where 
% each row is a different comparison

%% Get distributions ready

% Make sure distributions have the same number of bins, otherwise return
% error

if numel(p) ~= numel(q)
    error('Number of elements in the two input distributions must be equal');
end

% Make sure we're working with horizontal vectors if both inputs are vectos
if iscolumn(p)
    p = p';
end
if iscolumn(q)
    q = q';
end

% Make sure distributions are normalized
p = p./sum(p,2);
q = q./sum(q,2);


%% Calculate JSD

M = 0.5 * (p + q);

% Ignore cases where p and q are both exactly zero
DPM = sum(p .* log(p./M),2,'omitnan');
DQM = sum(q .* log(q./M),2,'omitnan');

JSD = 0.5 * (DPM + DQM);

end