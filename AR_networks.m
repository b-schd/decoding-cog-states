%% AR_Networks
% Builds networks by simultaneously solving a least squares optimization to
% generate a coefficient matrix

% Dependencies: 
%   * BCT Toolbox: https://sites.google.com/site/bctnet/
%   * arfit

% load parameters
params = proj_config();
addpath(genpath('Functions'))

sloc = load(params. subjChLoc) ; sloc = sloc.patient_loc(1).session;

errors= {}; 
errctr = 1;

for psubj = {'HUP128' 'HUP145' 'HUP184'}

psubj
ind = strcmp(pt_fold(i_pt).name, {sloc.subjID});
roi = sloc(ind).roi(sloc(ind).roi > 0);
[~,~, roi] = unique(roi);
n_roi = length(unique(roi)); 

try
    params.subj = psubj{1};
    savepth = fullfile(params.ndir, params.subj, 'AR_networks');
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
clear Networks Metrics Metrics_roi
iev_params = params; 
Metrics_roi = struct();
for ii = 1:Ntrl
    ii
    net = getARNets(squeeze(datwin(:,ii,:)),params.Lwin*fs);
    Networks(ii).net = net;
    Metrics(ii) = getMetsDIR(Networks(ii).net, {'strength', 'SSIndex',...
        'sourceIndex', 'sinkIndex'});
    
    % Get network averaged by roi
    if size(roi,1) ~= size(Networks(ii).net,1), continue, end
    roi_net = zeros(n_roi, n_roi, size(Networks(ii).net,3)); 

    for iii = unique(roi)'
        for jjj = unique(roi)'
            if iii == jjj, continue, end
                roi_net(iii, jjj,:) = mean(net(roi==iii, roi==jjj,:),[1,2]); 
        end
    end 
    
    Networks(ii).net_roi = roi_net; 
    Metrics_roi(ii) = getMetsDIR(Networks(ii).net_roi, {'strength', 'SSIndex',...
        'sourceIndex', 'sinkIndex'});
    
end

save(fullfile(savepth, sprintf('%s_net-ar_iev-%d_bp-%1.0f-%1.0f.mat', params.subj,...
    params.iev, params.bpfilt)), ...
    'Networks', 'Metrics', 'Metrics_roi', iev_params', '-V7.3')


        end
    end
    
catch ME
    disp(ME.message)
    clear MEs
    MEs.message = ME.message;
    MEs.stack = ME.stack;
    MEs.ptID = psubj;
    MEs.iev = iev; 
    MEs.savepth = savepth; 
    errors{errctr} = MEs;
    errctr = errctr+1; 
    continue; 
end
    
    
end
