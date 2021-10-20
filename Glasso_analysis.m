% Glasso Analysis

% TODO: 
% - utilize baseline data? As a regressor (after metric calculation,
%       check variance of metric value per node)
% - Find electrodes (or lobes) of interest? Regularized regression
% - See if there is a trend b/w electrode distance and edgeweight
% - How to not throw away so much data? Do we _need_ all time windows? 
% - Other ways of splitting data over time - higher order fluctuations?
% - remember: other patients exist... so don't lose hope :) 
%      - maybe stratify by improvement, or real value of RT

% %% LOAD DATA %% %

params = proj_config();

% Load Network and Metrics data for baseline, pre cue, and pre go
Baseline = load(fullfile(params.ddir, params.subj, 'networks',...
    sprintf('baseline_%s_%s.mat', params.subj, params.sess)));

preCue = load(fullfile(params.ddir, params.subj, 'networks',...
    sprintf('iev_1_%s_%s.mat', params.subj, params.sess)));

preGo = load(fullfile(params.ddir, params.subj, 'networks',...
    sprintf('iev_2_%s_%s.mat', params.subj, params.sess)));

% Load DataStruct
[dataStruct, Nsubj] = preprocessCCDT(params.ddir, params.subjChLoc, ...
    params.subj, params); 

% Load Loc data for patient. 
load(params.subjChLoc)
i_pt = strcmp(params.subj, {patient_loc(1).session.subjID});
ptChLoc = patient_loc(1).session(i_pt);

ichan = (ptChLoc.type~=0);
coords = ptChLoc.coords(ichan,:); 
cnames = ptChLoc.names(ichan,:); 

%% Plot Baseline 

% %% SETTINGS %% %
metric = 'aveCtrl'; %aveCtrl, modalCtrl, strength, etc.

nodalMetricMean = mean(Baseline.Metrics.(metric),2); % Nodal mean

figure(1); clf; % "Unrolled" upper triangle of each network over time
imagesc(Baseline.Networks.config_pcm)
title('Basline configuration matrix'); xlabel('time windows'); ylabel('edges')

for i = 1:20
    figure(2); clf; 
    imagesc(Baseline.Networks.pcm(:,:,i))
    title('Basline pcm'); xlabel('time windows'); ylabel('edges')
    pause(0.1)
end

figure(3); clf;
imagesc(Baseline.Metrics.(metric))
title(sprintf('Baseline %s', metric)); xlabel('time windows'); ylabel('channels')

figure(4); clf; 
imagesc(nodalMetricMean)
title(sprintf('Baseline mean nodal %s', metric));  ylabel('channels')

figure(5); clf; hold on; grid on; % Plot metric on 3D axis with spatial coordinates
scatter3(coords(:,1), coords(:,2), coords(:,3), rescale(nodalMetricMean, 50, 350), nodalMetricMean, 'filled')
textscatter3(coords(:,1), coords(:,2), coords(:,3), cnames,'TextDensityPercentage', 80)
cb = colorbar; cb.Label.String = metric;


%% Plot Trials

% %% SETTINGS %% %
metric = 'aveCtrl';   % Select metric to plot
evnt = preGo;        % Select Event type (eg. preCue, preGo)

nTrials = length(evnt.Metrics); 

% Metric average across nodes for each time window and trial
avgMetricWind = cell2mat(arrayfun(@(x) mean(evnt.Metrics(x).(metric)), [1:nTrials], 'Uni', 0)');
avgMetricNode = cell2mat(arrayfun(@(x) mean(evnt.Metrics(x).(metric),2), [1:nTrials], 'Uni', 0));

% for i_trial= 1:nTrials
%     figure(10); clf; 
%     imagesc(trialMetrics(i_trial))
%     title(sprintf('Trial %d %s', i_trial, metric)); xlabel('time windows')
%     pause(0.1)
% end

figure(11); 
plot(avgMetricWind, '.-', 'linewidth', .5, 'markersize', 20)
title(sprintf('%s Trend in pre-windows', metric));
xlabel('trial #'); ylabel(metric)
legend(sprintfc('window %d',[1:size(avgMetricWind, 2)]))

figure(12); clf; hold on;
imagesc(evnt.Networks(1).pcm(:,:,1))

figure(13); clf; hold on;
i_RT = dataStruct.RT > 0;
x = mean(avgMetricWind(i_RT,:),2); y = dataStruct.RT(i_RT);
[r, pval] = corr(x,y);
pf = polyfit(x,y,1);
scatter(mean(avgMetricWind,2), dataStruct.RT)
plot(x,polyval(pf,x),'-')
xlabel(metric); ylabel('recation time'); ylim([0,1000])
title(sprintf('%s vs. positive RT, r = %0.02f (p = %0.2f)', metric, r, pval))

% Correlate Each channel metric trajectory with RT
figure(14);
metricCorr = corr(avgMetricNode(:, i_RT)', y);
[srt, i_srt] = sort(metricCorr);
plot(srt)
ylabel(sprintf('r', metric))
xticks([1:length(cnames)]); xlabel('Channels/Electrodes')
xticklabels(cnames(i_srt)); xtickangle(90); axis tight;
title(sprintf('Correlation coefficient of nodal %s and RT across trials', metric))

%3D plot of nodes that have highest correlation with RT
figure(15); clf; hold on; grid on; % Plot metric on 3D axis with spatial coordinates
scatter3(coords(:,1), coords(:,2), coords(:,3), rescale(metricCorr, 50, 350), metricCorr, 'filled')
textscatter3(coords(:,1), coords(:,2), coords(:,3), cnames,'TextDensityPercentage', 80)
cb = colorbar; cb.Label.String = metric;
title(sprintf('Nodal correlation of %s with RT', metric))


