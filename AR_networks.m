%% AR_Networks
% Builds networks by simultaneously solving a least squares optimization to
% generate a coefficient matrix

% Dependencies: 
%   * BCT Toolbox: https://sites.google.com/site/bctnet/
%   * arfit

% load parameters
params = proj_config();
addpath(genpath('Functions'))

for psubj = {'HUP069', 'HUP142','HUP157', 'HUP133'} 

    params.subj = psubj{1};
    savepth = fullfile(params.ddir, params.subj, 'networks');
    filtList =[[3 12]; [14 30]; [35 55]; [70 150]];

    for iev = [1,2]
        params.iev = iev;
    
        for i = 1:length(filtList)
            params.bpfilt = filtList(i,:); 

% get subject data and event info
[dataStruct, Nsubj] = preprocessCCDT(params.ddir, params.subjChLoc, ...
    params.subj, params);
Ntrl = dataStruct.Ntrl;
fs = dataStruct.fs;
bDatwin = dataStruct.bDatwin;
datwin = dataStruct.datwin;

% Compute networks for trial an windows specified in params
clear Networks Metrics 
iev_params = params; 
for ii = 1:Ntrl
    ii
    Networks(ii) = getARNets(squeeze(datwin(:,ii,:)),params.Lwin*fs);
    Metrics(ii) = getMetsDIR(Networks(ii), {'strength', 'SSIndex',...
        'sourceIndex', 'sinkIndex'});
end

save(fullfile(savepth, sprintf('%s_net-ar_iev-%d_bp-%1.0f-%1.0f.mat', params.subj,...
    params.iev, params.bpfilt)), ...
    'Networks', 'Metrics', 'iev_params', '-V7.3')


        end
    end
end
