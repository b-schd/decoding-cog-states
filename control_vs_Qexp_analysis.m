% parameters
cleanTrials = true; %if true, remove all trials with negative reaction times
freqBands = {'Theta/Alpha', 'Beta', 'Low Gamma', 'High Gamma'};
plrty = 1; %polarity can be 1 or 2 for channel location file

% %% LOAD DATA %% %
params = proj_config();
load(params.LTdata);

% Load DataStruct
[dataStruct, Nsubj] = preprocessCCDT(params.ddir, params.subjChLoc, ...
    params.subj, params); 

% Load Network and Metrics data for precue period
preCue = load(fullfile(params.ndir, params.subj, 'networks',...
    sprintf('iev_1_%s_.mat', params.subj)));
%{
preGo = load(fullfile(params.ndir, params.subj, 'networks',...
    sprintf('iev_2_%s_.mat', params.subj)));
%}

% Load SEEG location data 
load(params.subjChLoc);
db = CCDTdatabase;
subjInd = cellfun(@(x) strcmp(x,params.subj),db(:,1)); % find specified subj
if ~isempty(params.sess)
    subjInd2 = cellfun(@(x) strcmp(x,params.sess),db(:,2)); % find specified session
    subjInd = subjInd & subjInd2;
end
subjInd = find(subjInd);
if isempty(subjInd), error('subj/sess not found in database'); end
locInfo = patient_loc(plrty).session(subjInd);

% set up average controllability matrices
if cleanTrials
    trialIdxs = find(dataStruct.RT >= 0);
else
    trialIdxs = 1:dataStruct.Ntrl;
end
Ntrl = length(trialIdxs);
Nch = size(preCue.Metrics(1).aveCtrl,1);
precueControl = ones(Ntrl,Nch,size(preCue.Metrics(1).aveCtrl,2));
%pregoControl = ones(Ntrl,Nch);
for i = 1:Ntrl
    precueControl(i,:,:) = preCue.Metrics(trialIdxs(i)).aveCtrl;
    %pregoControl(i,:) = preGo.Metrics(trialIdxs(i)).aveCtrl(:,end);
end

% find channel indices for gray vs. white matter
grayIdxs = find(locInfo.type == 1);
whiteIdxs = find(locInfo.type == -1);
GWIdxs = {grayIdxs, whiteIdxs};
GWLabels = {'Gray', 'White'};

%find channel indices for 'significant' nodes
sigGWQexpIdxs = cell(length(freqBands),2);
if size(LTqexp{subjInd}, 1) ~= Nch, error(['Number of channels in ' ...
        'the Qexp and controllability datasets are not equal']); end
for i=1:length(freqBands)
    sigQexpIdxs = find(LTqexp{subjInd}(:,i,4) < 0.05);
    sigGWQexpIdxs{i,1} = intersect(sigQexpIdxs, GWIdxs{1});
    sigGWQexpIdxs{i,2} = intersect(sigQexpIdxs, GWIdxs{2});
end

boxplotBySubj(precueControl, {'All', 'Average Controllability'}, ...
    GWLabels, GWIdxs, Ntrl, params);
for i=1:length(freqBands)
    boxplotBySubj(precueControl, ...
        {string([freqBands{i}, ' Significant']), 'Average Controllability'}, ...
        GWLabels, sigGWQexpIdxs(i,:), Ntrl, params); 
end

% Plot controllability according to gray vs. white matter
function boxplotBySubj(data, titles, groupLbls, groupIdxs, Ntrl, params)
for i = 1:length(groupLbls)
    figure;
    thisData = data(:,groupIdxs{i},end);
    vecData = reshape(thisData,[],1);
    vecChLabels = repelem(groupIdxs{i},Ntrl);
    boxplot(vecData,vecChLabels)
    titleString = string([titles{2} ' for ' titles{1} ' ' groupLbls{i} ...
        ' Matter Channels in Precue Period' newline params.subj]);
    if ~isempty(params.sess)
        titleString = titleString + ' ' + params.sess;
    end
    title(titleString)
    xlabel('Channel Number')
    ylabel(titles{2})
end
end

