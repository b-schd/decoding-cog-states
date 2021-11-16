% Create networks- Creates GLASSO networks and comuptes metrics on them
% Dependencies: 
%   * BCT Toolbox: https://sites.google.com/site/bctnet/
%   * "Rscript" command accessible from the commandline
%   * Library requirements in R: bootnet, R.matlab

% load parameters
params = proj_config();
savepth = fullfile(params.ddir, params.subj, 'networks');

% get subject data and event info
[dataStruct, Nsubj] = preprocessCCDT(params.ddir, params.subjChLoc, ...
    params.subj, params);
Ntrl = dataStruct.Ntrl;
fs = dataStruct.fs;
bDatwin = dataStruct.bDatwin;
datwin = dataStruct.datwin;

% Create networks from each time window of length Lwin, and compute control
% metrics
% [baselineNetworks] = getNets(params.glassoPath, bDatwin, params.Lwin*fs, params.gamma, params.beta); 
% [baselineMetrics] = getMets(baselineNetworks);'
% base_params = params;
% 
% save(fullfile(savepth, sprintf('baseline_%s_%s.mat', params.subj, params.sess)),...
%    'Networks', 'Metrics', 'base_params', '-V7.3');


% Compute networks for trial an windows specified in params
tic
clear Networks Metrics 
iev_params = params; 
for ii = 1:Ntrl
    ii
     Networks(ii) = getNets(params.glassoPath, squeeze(datwin(:,ii,:)), ...
        params.Lwin*fs, params.gamma, params.beta);
    Metrics(ii) = getMets(Networks(ii));
    
    save(fullfile(savepth, sprintf('iev-%d_%s%s_bp-%1.0f-%1.0f.mat', params.iev,...
        params.subj, params.sess, params.bpfilt)), ...
        'Networks', 'Metrics', 'iev_params', '-V7.3')
end
toc


