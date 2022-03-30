%determines how skewed the distribution of edge weights is
%
%Background: This function is utilized because sometimes the vast majority of non-zero edge weights in
% may be very close in value to zero (with a few outliers close to 1) i.e.
%highly skewed.  A more evenly distributed edge weights improves cluster
%detection.
%
%Practically rarely invoked on real netwokrs, but sometimes triggered in
%LFR benchmarks.
%
%This should be considered in conjunction with the othmer situation in which some edges are reweights
%which occurs when self-connections will set to 1 (i.e. main
%diagonal of adjacency matrix is all 1) which may be much higher than any non-self
%edges (different code, but mention here as it's related, and would yeild partition where each node is it's own cluster).  
function output= skewns(x)
output=(sum((x-mean(x)).^3)./length(x)) ./ (var(x,1).^1.5);
