% Create networks- Creates GLASSO networks and comuptes metrics on them
% Dependencies: 
%   * BCT Toolbox: https://sites.google.com/site/bctnet/
%   * "Rscript" command accessible from the commandline
%   * Library requirements in R: bootnet, R.matlab

% load parameters
params = proj_config();
savepth = fullfile(params.ddir, params.subj, 'control', sprintf('cntrl_data_%s_%s.mat', params.subj, params.sess)); 

% get subject data and event info
[dataStruct, Nsubj] = preprocessCCDT(params.ddir, params.subjChLoc, ...
    params.subj, params);
Ntrl = dataStruct.Ntrl;
fs = dataStruct.fs;
bDatwin = dataStruct.bDatwin;
datwin = dataStruct.datwin;

% Create networks from each time window of length Lwin, and compute control
% metrics
[baselineNetworks] = getNets(params.glassoPath, bDatwin, params.Lwin*fs, params.gamma, params.beta); 
[baselineMetrics] = getMets(baselineNetworks);

save(savepth, 'params', 'baselineNetworks', 'baselineMetrics', '-V7.3')


clear trialNetworks trialMetrics 
for ii = 1:Ntrl
     trialNetworks(ii) = getNets(params.glassoPath, squeeze(datwin(:,ii,:)), ...
        params.Lwin*fs, params.gamma, params.beta);
    trialMetrics(ii) = getMets(trialNetworks(ii));
    
    save(savepth, 'trialNetworks', 'trialMetrics', '-V7.3', '-append')
end



