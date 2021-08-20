% Create networks- Creates GLASSO networks and comuptes metrics on them
% Dependencies: 
%   * BCT Toolbox: https://sites.google.com/site/bctnet/
%   * "Rscript" command accessible from the commandline
%   * Library requirements in R: bootnet, R.matlab


% load parameters
params = proj_config();

% get subject data and event info
[dat,Nsamp,fs,chn,chnm] = loadCCDTdata(params);
[ccdt,fs] = parseCCDTevents(params);

% Example: Split data into 5 second time windows
Lwin = fs*5; 
glassoPath = './'; 

% Create networks from each 5 second time window, and compute control
% metrics
[Networks] = getNets(glassoPath, dat', Lwin, .25, 1); 
[metric_matrices] = getMets(Networks);

savepth = fullfile(params.ddir, 'control', p.subj, sprintf('cntrl_data_%s_%s.mat', p.subj, p.sess)); 

save(savepth, 'params', 'Networks', 'metric_matrices')
