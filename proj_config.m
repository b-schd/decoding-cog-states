
function p = proj_config()

p.ddir = '/Users/bscheid/Documents/LittLab/PROJECTS/p15_Cog_Performance/CCDT/eeg/'; % data directory
p.subj = 'HUP069'; %'HUP142'; %; % subject
p.sess = '04Oct17'; %'30Jun17';% ; % session date
p.reref = 'noreref';
p.stime = []; % session time

p.rln = 1; % remove line noise? (0 or 1)
p.rrf = 2; % re-reference? (0 = none, 1 = common average re-reference, 2 = bipolar re-reference, 3 = laplacian re-reference)
p.outl = 1; % remove outlier channels? (0 or 1)
p.outlMetric = 'powGamma'; % outlier metric ('rms' = root-median-squared, 'powGamma' = high-gamma power)
p.outlThresh = 5; % outlier threshold (standard deviations)


end