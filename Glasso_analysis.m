% Glasso Analysis

% TODO: 
% - utilize baseline data? As a regressor (before or after metric calculation?)
% - Figure out patient mapping for patient_loc --> add an ID column?
% - Find electrodes (or lobes) of interest? 
% - See if there is a trend b/w electrode distance and edgeweight
% - How to not throw away so much data? Do we _need_ all time windows? 
% - Other ways of splitting data over time
% - remember: other patients exist... so don't lose hope :) 
%      - maybe stratify by improvement, or real value of RT


params = proj_config();
loadpth = fullfile(params.ddir, params.subj, 'control',...
    sprintf('cntrl_data_%s_%s.mat', params.subj, params.sess));

load(loadpth)
[dataStruct, Nsubj] = preprocessCCDT(params.ddir, params.subjChLoc, ...
    params.subj, params); 

%% Plot baseline 

figure(1); clf; 
imagesc(baselineNetworks.config_pcm)
title('Basline configuration matrix'); xlabel('time windows')

for i = 1:20
    figure(10); clf; 
    imagesc(baselineNetworks.pcm(:,:,i))
    title('Basline pcm'); xlabel('time windows')
    pause
end

metric = 'aveCtrl'; %aveCtrl, modalCtrl, strength, etc.

figure(2); clf;
imagesc(baselineMetrics.(metric))
title(sprintf('Baseline %s', metric)); xlabel('time windows'); ylabel('channels')

figure(3); clf; 
plot(mean(baselineMetrics.(metric),2))
title(sprintf('Baseline mean nodal %s', metric));  ylabel('channels')


%% Plot trials

metric = 'aveCtrl';
nTrials = length(trialMetrics); 

for i_trial= 1:nTrials
    figure(10); clf; 
    imagesc(trialMetrics(i_trial))
    title(sprintf('Trial %d %s', i_trial, metric)); xlabel('time windows')
    pause
end

avgMetric = cell2mat(arrayfun(@(x) mean(trialMetrics(x).(metric)), [1:nTrials], 'Uni', 0)');

figure(11); 
plot(avgMetric, '.-', 'linewidth', .5, 'markersize', 20)
title(sprintf('%s Trend in 5 pre-windows', metric));
xlabel('trial #'); ylabel(metric)
legend(sprintfc('%d',[1:5]))

figure(12); clf; hold on;
imagesc(trialNetworks(1).pcm(:,:,1))

figure(13); clf; hold on;
i_RT = dataStruct.RT > 0;
x = mean(avgMetric(i_RT,:),2); y = dataStruct.RT(i_RT);
[r, pval] = corr(x,y);
pf = polyfit(x,y,1);
scatter(mean(avgMetric,2), dataStruct.RT)
plot(x,polyval(pf,x),'-')
xlabel(metric); ylabel('recation time'); ylim([0,1000])
title(sprintf('%s vs. positive RT, r = %0.02f (p = %0.2f)', metric, r, pval))


