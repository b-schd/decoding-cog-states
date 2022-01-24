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

% - PCA on feature matrix, do certain windows cluster? <-- Brittany 
% - look at windows leading up to preCue or preGo <-- Brittany (deltas, N
% channels) (which channels or groups of channels change, reaction time)
% - post-go windows? prego/postgo transistion <-- 
% - linear regression stuff.... <- Grace 
% - Look at variance
% - 


% %% LOAD DATA %% %

params = proj_config();
net_type = 'ar';

% Load Network and Metrics data for baseline, pre cue, and pre go
try
Baseline = load(fullfile(params.ddir, params.subj, 'networks',...
    sprintf('net-%s_baseline_%s_%s.mat', params.subj, net_type, params.sess)));
catch
end

try
preCue = load(fullfile(params.ddir, params.subj, 'networks',...
    sprintf('net-%s_iev_1_%s_%s.mat', params.subj, net_type, params.sess)));
catch; end

try
preGo = load(fullfile(params.ddir, params.subj, 'networks',...
    sprintf('net-%s_iev_2_%s_%s.mat',params.subj, net_type, params.sess)));
catch; end

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


%% %% SETTINGS %% %%

metric = 'aveCtrl';   % Select metric to plot
evnt = preCue;        % Select Event type (eg. preCue, preGo)

nTrials = length(evnt.Metrics); 
i_RT = dataStruct.RT > 0;

%% load data to table for LME

clear preCueData preGoData
for i_tr = 1:nTrials; preCueData(:,:,i_tr) = preCue.Metrics(i_tr).(metric); end
[m,n,p]= size(preCueData); [x,y,z] = ndgrid(1:m,1:n,1:p);
preCue_tbl = array2table([x(:),y(:),z(:),dataStruct.RT(z(:)), dataStruct.DT(z(:)), preCueData(:)],...
    'VariableNames', {'channel', 'window', 'trial', 'RT', 'DT', metric});
preCue_tbl.event = repmat({'preCue'}, height(preCue_tbl), 1);

for i_tr = 1:nTrials; preGoData(:,:,i_tr) = preGo.Metrics(i_tr).(metric); end
[m,n,p]= size(preGoData); [x,y,z] = ndgrid(1:m,1:n,1:p);
preGo_tbl = array2table([x(:),y(:),z(:),dataStruct.RT(z(:)), dataStruct.DT(z(:)), preGoData(:)],...
    'VariableNames', {'channel', 'window', 'trial', 'RT', 'DT', metric});
preGo_tbl.event = repmat({'preGo'}, height(preGo_tbl), 1);

dataTbl = vertcat(preCue_tbl, preGo_tbl);

%writetable(dataTbl, fullfile(params.ddir, params.subj, 'networks', 'lmetable.csv')) 


%% Plot Baseline 

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


% Metric average across nodes for each time window and trial
avgMetricWind = cell2mat(arrayfun(@(x) mean(evnt.Metrics(x).(metric)), [1:nTrials], 'Uni', 0)');
avgMetricNode = cell2mat(arrayfun(@(x) mean(evnt.Metrics(x).(metric),2), [1:nTrials], 'Uni', 0));


figure(9); clf;
z= cell2mat({evnt.Metrics.(metric)});
imagesc(z)

% for i_trial= 1:nTrials
%     figure(10); clf; 
%     imagesc(evnt.Metrics(i_trial).(metric))
%     title(sprintf('Trial %d %s', i_trial, metric)); xlabel('time windows')
%     pause(0.1)
% end

figure(11); 
plot(avgMetricWind, '.-', 'linewidth', .5, 'markersize', 20)
title(sprintf('%s Trend in pre-windows', metric));
xlabel('trial #'); ylabel(metric)
legend(sprintfc('window %d',[1:size(avgMetricWind, 2)]))
yyaxis right 
plot(dataStruct.RT, '.-', 'linewidth', .5, 'markersize', 20)

figure(12); clf; hold on;
imagesc(evnt.Networks(1).pcm(:,:,1))

figure(13); clf; hold on;
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


%%
featMat1 = [avgMetricNode(:,i_RT); dataStruct.RT(i_RT)']';
featMat2 = [avgMetricWind(i_RT,:), dataStruct.RT(i_RT)];

% iteratively remove electrode with highest # of nans
aw = corrcoef(avgMetricNode(:,i_RT)').*~eye(size(avgMetricNode,1)); 
thresh = 0.6;

i_elecs = [1:size(avgMetricNode,1)];
mx = Inf;
while mx > 1
    sum(aw > thresh)
    [mx,i_m] = max(sum(aw > thresh));
    i_elecs(i_m) = [];
    aw = corrcoef(avgMetricNode(i_elecs,i_RT)').*~eye(length(i_elecs)); 
    imagesc(aw); colorbar
    pause(.01)
end
    


featMat3 = featMat1(:,[i_elecs, end]);
B = lasso(featMat3(:,1:end-1), featMat3(:,end));
figure(21); imagesc(B)

featMat4 = featMat3(:,[B(:,70)~=0; true]);
figure(22); crosscorr(mean(featMat4(:,1:end-1),2), featMat4(:,end))
lag = 15;
figure(23); scatter(mean(featMat4(1+lag:end,1:end-1),2), featMat4(1:end-lag,end)); title(sprintf('evnt %d, %s', params.iev, metric))

[rr, pp] = corr(mean(featMat4(1+lag:end,1:end-1),2), featMat4(1:end-lag,end))

featMat5 = [featMat4(1+lag:end,1:end-1), featMat4(1:end-lag,end)];

%%

% Look into linear change

wins = mean(preCueData,3);
wins = wins - mean(wins,2);
w_std = std(preCueData,[],3);

figure(40);clf;
imagesc(wins); ylabel('channels'), xlabel('window #'); title('averaged across trials and referenced using mean across windows')

[rr, pp] = corr(wins', [1:5]');
fprintf('%d channels at p < 0.05 \n', sum(pp<0.05));
find(pp<0.05)

figure(42); histogram(rr, 10); title('histogram of correlation b/w metric value and window number')

figure(41); clf;
errorbar(wins(pp<0.05,:)', w_std(pp<0.05,:)', 'linewidth', 2); 
xlabel('Window #'); ylabel('avgCtrl'); title('Channels with significant linear correlation')

avgslopes = [];
for i_chan = 1:length(wins)
    ps = polyfit([1:5], wins(i_chan,:), 1);
    avgslopes(i_chan) = ps(1);
end

slopes = [];
for i_chan = 1:length(wins)
    for i_trial = 1:nTrials
        ps = polyfit([1:5], preCueData(i_chan,:,i_trial), 1);
        slopes(i_chan, i_trial) = ps(1); 
    end
end

figure(38); clf; 
imagesc(slopes); xlabel('trial #'), ylabel('channel'); title('slopes')
colorbar

[rr, pp] = corr(slopes', dataStruct.RT);
fprintf('%d channel''s slopes significantly correlate with reaction time \n', sum(pp<0.05))
find(pp<0.05)

figure(39); clf;
histogram(rr, 15); title('histogram of correlation coeff. between slopes and RT'); 
xlabel('r'); ylabel('#channels')


%3D plot of nodes that have highest correlation with RT
figure(15); clf; hold on; grid on; % Plot metric on 3D axis with spatial coordinates
scatter3(coords(:,1), coords(:,2), coords(:,3), rescale(rr, 50, 350), rr, 'filled')
textscatter3(coords(:,1), coords(:,2), coords(:,3), cnames,'TextDensityPercentage', 80)
cb = colorbar; cb.Label.String = metric;
title(sprintf('Nodal correlation of slope with RT', metric))


win2 = diff(mean(preCueData,3),2)-min(diff(mean(preCueData,3),2),[],2)
imagesc(win2)


%% PCA

figure(3); 

[coeff, score, ~, ~, expl] = pca(avgMetricNode');

figure(1);
scatter(score(:,1), score(:,2), 50, rescale(dataStruct.RT,0,255), 'filled')

figure(2); plot(cumsum(expl))


