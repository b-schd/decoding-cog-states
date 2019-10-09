%% Define load and save directories
loadDir = '/Users/graceng/Documents/Med_School/MS1/Research/BCI/CCDT mat data';
saveDir = '/Users/graceng/Documents/Med_School/MS1/Research/BCI/CCDT mat data/Plots_030819/indiv';

%% Define analysis parameters
periods = {'precueone' 'precue' 'DTone' 'DTtwo' 'DTthree'};
minRT = 0; %the minimum reaction time in order for a trial to qualify as 
            % "no-error"
            
%typeTrial = "E"; %if typeTrial is "E", then we analyze the error trials
typeTrial = "noE"; %if typeTrial is "noE", the we analyze the non-error trials
            
DTThresh = 50; %trials with delay times within this threshold are considered in this analysis

DT = 1500; %the delay time of the trials to be considered in this analysis
typeDT = "min"; %if typeDT is 'min', then we consider trials with DTs greater than or equal to DT-DTThresh

%DT = 500;
%typeDT = "max"; %if typeDT is 'max', then we consider trials with DTs less
                    %than or equal to DT+DTThresh

freqInd = 2; %the index for the frequency band data 
                %(1 for theta, 2 for high gamma, 3 for low gamma)
freqType = 'HG';

%freqInd = 1;
%freqType = 'theta';

smoothFac = 5; %smoothing factor for plotting
sortType = "Sorted";
%sortType = "Chronological";

scaleMethod = "scale";
%scaleMethod = "norm";
subplotDim = {4, 5};
%save = "combined";
save = "indiv";

%% load data
load([loadDir '/allSubjBehavAndNetworkData'])
load([loadDir '/allSubjStats_n19'],['subjOrder'])

%% create and save sorted netStr plots
fig1 = figure;
fig2 = figure;
fig3 = figure;
for i=1:length(subjOrder)
    vRT = allSubjBehaviorStruct{i}.vRT;
    vDT = allSubjBehaviorStruct{i}.vDT;
    if typeDT == "min"
        indNoE = vRT>minRT&vDT>(DT-DTThresh);
    elseif typeDT == "max"
        indNoE = vRT>minRT&vDT<(DT+DTThresh);
    else
        error('typeDT must either be "min" or "max".');
    end
    
    if typeTrial == "E"
        ind = ~indNoE;
    elseif typeTrial == "noE"
        ind = indNoE;
    else
        error('typeTrial must either be "E" or "noE".');
    end
    netStr = zeros(sum(ind), length(periods));
    for j=1:length(periods)
        netStr(:,j) = mean(allSubjNetMetrics{i}.(['str_' periods{j}])((ind),:,freqInd),2);
    end
    if sortType == "Sorted"
        label = "Sorted by RT";
        [~, indSortRT] = sort(vRT(ind));
        plotNetStr = netStr(indSortRT,:);
    elseif sortType == "Chronological"
        label = "Chronological";
        plotNetStr = netStr;
    else
        error('sortType must be a cell array containing the strings "Sorted" or "Chronological"');
    end
    figure(fig1);
    hax = subplot(subplotDim{1},subplotDim{2}, i);
    imagesc(fgsmooth(plotNetStr,smoothFac))
    colorbar
    xticks([1 2 3 4 5])
    xticklabels(periods)
    xtickangle(45)
    xlabel('Periods')
    ylabel(strjoin([sprintf('%d',DT) 'ms trials, ' label])) 
    set(gcf,'color','w');
    hold on
    if save == "indiv"
        hfig = figure;
        hax_new = copyobj(hax, hfig);
        set(hax_new, 'Position', get(0, 'DefaultAxesPosition'));
        title(strjoin([sprintf('%s ',subjOrder{i}) freqType ' netStr, ' label]));
        filename = strjoin([saveDir '\' sprintf('%s',subjOrder{i}) '_' freqType 'netStr' sortType '_' sprintf('%d',DT) 'ms'],'');
        print(hfig, filename,'-dpng')
    else
        title(sprintf('%s ',subjOrder{i}))
    end

    %scale within period
    if scaleMethod == "scale"
        colmin = min(plotNetStr);
        colmax = max(plotNetStr);
        scalePlotNetStr = rescale(plotNetStr,'InputMin',colmin,'InputMax',colmax);
    elseif scaleMethod == "norm"
        scalePlotNetStr = normalize(plotNetStr,1);
    else
        error('scaleMethod must either be "scale" or "norm".')
    end
    figure(fig2);
    hax = subplot(subplotDim{1},subplotDim{2}, i);
    imagesc(fgsmooth(scalePlotNetStr,smoothFac))
    colorbar
    xticks([1 2 3 4 5])
    xticklabels(periods)
    xtickangle(45)
    xlabel('Periods')
    ylabel(strjoin([sprintf('%d',DT) 'ms Trials, ' label]))
    set(gcf,'color','w');
    hold on
    if save == "indiv"
        hfig = figure;
        hax_new = copyobj(hax, hfig);
        set(hax_new, 'Position', get(0, 'DefaultAxesPosition'));
        filename = strjoin([saveDir '\' sprintf('%s',subjOrder{i}) '_' freqType 'NetStr_' scaleMethod 'Period_' sortType '_' sprintf('%d',DT) 'ms'], '');
        title(hax_new, strjoin([sprintf('%s ',subjOrder{i}) freqType ' netStr (Scaled/Period, ' sortType, ')']))
        print(hfig, filename,'-dpng')
    else
        title(sprintf('%s ',subjOrder{i}))
    end

    %scale within trial
    if scaleMethod == "scale"
        rowmin = min(plotNetStr,[],2);
        rowmax = max(plotNetStr,[],2);
        scalePlotNetStr = rescale(plotNetStr,'InputMin',rowmin,'InputMax',rowmax);
    elseif scaleMethod == "norm"
        scalePlotNetStr = normalize(plotNetStr,2);
    else
        error('scaleMethod must either be "scale" or "norm".')
    end
    figure(fig3);
    hax = subplot(subplotDim{1},subplotDim{2}, i);
    imagesc(fgsmooth(scalePlotNetStr,smoothFac))
    colorbar
    xticks([1 2 3 4 5])
    xticklabels(periods)
    xtickangle(45)
    xlabel('Periods')
    ylabel(strjoin([sprintf('%d',DT) 'ms trials, ' label]))
    set(gcf,'color','w');
    hold on
    title(sprintf('%s ',subjOrder{i}))
    if save == "indiv"
        hfig = figure;
        hax_new = copyobj(hax, hfig);
        set(hax_new, 'Position', get(0, 'DefaultAxesPosition'));
        title(hax_new, strjoin([sprintf('%s ',subjOrder{i}) freqType ' netStr (Scaled/Trial, ' sortType ')']))
        filename = strjoin([saveDir '\' sprintf('%s',subjOrder{i}) '_' freqType 'NetStr_' scaleMethod 'Trial_' sprintf('%d',DT) 'ms'], '');
        print(hfig, filename,'-dpng')
    end
    
    clear netStrNoE vRT vDT sortRTNoE indSortRTNoE sortNetStrNoE scaleSortNetStrNoE
end

if save == "combined"
    figure(fig1);
    filename = strjoin([saveDir '/' freqType '_NetStr_' sprintf('%d',DT) 'ms_' typeTrial '_' sortType],'');
    sgtitle(strjoin([freqType ' netStr, ' label]));
    set(gcf, 'Position', get(0, 'Screensize'));
    print(filename,'-dpng')
    figure(fig2);
    filename = strjoin([saveDir '/' freqType '_NetStr_' sprintf('%d',DT) 'ms_scalePeriod_' typeTrial '_' sortType],'');
    sgtitle(strjoin([freqType ' netStr (Scaled Within Period), ' label]));
    set(gcf, 'Position', get(0, 'Screensize'));
    print(filename,'-dpng')
    figure(fig3);
    filename = strjoin([saveDir '/' freqType '_NetStr_' sprintf('%d',DT) 'ms_scaleTrial_' typeTrial '_' sortType],'');
    sgtitle(strjoin([freqType ' netStr (Scaled Within Trial), ' label]));
    set(gcf, 'Position', get(0, 'Screensize'));
    print(filename,'-dpng')
end