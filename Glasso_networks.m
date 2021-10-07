% Create networks- Creates GLASSO networks and comuptes metrics on them
% Dependencies: 
%   * BCT Toolbox: https://sites.google.com/site/bctnet/
%   * "Rscript" command accessible from the commandline
%   * Library requirements in R: bootnet, R.matlab


% load parameters
params = proj_config();

% get subject data and event info

% TODO: write a new script that extracts data into the correct trial format,
% and that also gets 10 sec of the baseline recording data out
[dataStruct, Nsubj] = preprocessCCDT(params.ddir, params.subjChLoc, ...
    params.subj, params);
Ntrl = dataStruct.Ntrl;
bDatwin = dataStruct.bDatwin;
datwin = dataStruct.datwin;

% Create networks from each time window of length Lwin, and compute control
% metrics
[baselineNetworks] = getNets(params.glassoPath, bDatwin, params.Lwin, params.gamma, params.beta); 
[baselineMetrics] = getMets(baselineNetworks);
% TODO: need to save the data for the baseline recording period

for ii = 1:Ntrl
    [trialNetworks] = getNets(params.glassoPath, squeeze(datwin(:,ii,:)), ...
        params.Lwin, params.gamma, params.beta);
    [trialMetrics] = getMets(trialNetworks);
    % TODO: need to save the data for the trial
end

savepth = fullfile(params.ddir, 'control', p.subj, sprintf('cntrl_data_%s_%s.mat', p.subj, p.sess)); 

save(savepth, 'params', 'Networks', 'metric_matrices')
