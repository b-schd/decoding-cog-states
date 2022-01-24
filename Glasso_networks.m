% Create networks- Creates GLASSO networks and comuptes metrics on them
% Dependencies: 
%   * BCT Toolbox: https://sites.google.com/site/bctnet/
%   * "Rscript" command accessible from the commandline
%   * Library requirements in R: bootnet, R.matlab

% load parameters
params = proj_config();
addpath(genpath('Functions'))
savepth = fullfile(params.ddir, params.subj, 'networks');

filtList =[[3 12]; [14 30]; [35 55]; [70 150]];

for i = 1:length(filtList)
    
    params.bpfilt = filtList(i,:); 

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
     Networks(ii) = getGLASSONets(params.glassoPath, squeeze(datwin(:,ii,:)), ...
        params.Lwin*fs, params.gamma, params.beta);
    Metrics(ii) = getMetsUND(Networks(ii));
    
    save(fullfile(savepth, sprintf('%s_iev-%d_bp-%1.0f-%1.0f.mat', params.subj,...
        params.iev, params.bpfilt)), ...
        'Networks', 'Metrics', 'iev_params', '-V7.3')
end
toc


end
