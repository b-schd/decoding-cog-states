function [Networks]=getARNets(data, lWin)
% This script calls the getICov.R script located in EC_glasso to retrieve
% regularized inverse covariance and partial correlation matrices across windows of 
% channels recordings in "data". 

% INPUTS:
% data: Nxl data matrix, N channels, l-time samples
% lWin: time window to create network over (in samples)

% OUTPUT:
% Networks struct with the following fields: 
%   net: (NxNxT) set of reg. partial correlation matrices for each time window
%   config_net: (N(N-1)/2)xT vectorized set of pcm matrices

[nchans, T] = size(data); 
nWins = floor(T./(lWin));

As = zeros(nchans, nchans, nWins);

for w = 1:nWins 
    dat_win = data(:, (w-1)*lWin+[1:lWin]);
    [~, A, ~]= arfit(dat_win', 1,1, 'zero');
    A(logical(eye(nchans)))= 0;   
    
    As(:,:, w) = abs(A);

end

Networks.net = As; 

