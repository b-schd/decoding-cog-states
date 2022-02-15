function NaN = dispPreCueStats(stats)
%function NaN = dispPreCueStats(stats)
    %   Display stats result for PreCue period
%load('/Users/graceng/Documents/Med_School/Research/BCI/control_data/SVM_Stats/SVM_HUP069_13-Feb-2022.mat')
y = zeros(16, 3);
y(1,:) = stats.preCue.HighGamma.strength.auc;
y(2,:) = stats.preCue.HighGamma.SSIndex.auc;
y(3,:) = stats.preCue.HighGamma.sourceIndex.auc;
y(4,:) = stats.preCue.HighGamma.sinkIndex.auc;
y(5,:) = stats.preCue.LowGamma.strength.auc;
y(6,:) = stats.preCue.LowGamma.SSIndex.auc;
y(7,:) = stats.preCue.LowGamma.sourceIndex.auc;
y(8,:) = stats.preCue.LowGamma.sinkIndex.auc;
y(9,:) = stats.preCue.Beta.strength.auc;
y(10,:) = stats.preCue.Beta.SSIndex.auc;
y(11,:) = stats.preCue.Beta.sourceIndex.auc;
y(12,:) = stats.preCue.Beta.sinkIndex.auc;
y(13,:) = stats.preCue.ThetaAlpha.strength.auc;
y(14,:) = stats.preCue.ThetaAlpha.SSIndex.auc;
y(15,:) = stats.preCue.ThetaAlpha.sourceIndex.auc;
y(16,:) = stats.preCue.ThetaAlpha.sinkIndex.auc;

%plot
xplot = 1:16;
xconf = [xplot xplot(end:-1:1)];
yhigh = y(:,2)';
ylow = y(:,3)';
yconf = [yhigh ylow(end:-1:1)];
yplot = y(:,1)';

figure
p = fill(xconf, yconf, 'red');
p.FaceColor = [1 0.8 0.8];
p.EdgeColor = 'none';

hold on
plot(xplot, yplot, 'ro')
yline(0.5)
hold off
end