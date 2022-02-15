function NaN = dispPreGoStats(stats)
%function NaN = dispPreGoStats(stats)
    %   Display stats result for PreGo period
%load('/Users/graceng/Documents/Med_School/Research/BCI/control_data/SVM_Stats/SVM_HUP069_13-Feb-2022.mat')
y = zeros(16, 3);
y(1,:) = stats.preGo.HighGamma.strength.auc;
y(2,:) = stats.preGo.HighGamma.SSIndex.auc;
y(3,:) = stats.preGo.HighGamma.sourceIndex.auc;
y(4,:) = stats.preGo.HighGamma.sinkIndex.auc;
y(5,:) = stats.preGo.LowGamma.strength.auc;
y(6,:) = stats.preGo.LowGamma.SSIndex.auc;
y(7,:) = stats.preGo.LowGamma.sourceIndex.auc;
y(8,:) = stats.preGo.LowGamma.sinkIndex.auc;
y(9,:) = stats.preGo.Beta.strength.auc;
y(10,:) = stats.preGo.Beta.SSIndex.auc;
y(11,:) =stats.preGo.Beta.sourceIndex.auc;
y(12,:) = stats.preGo.Beta.sinkIndex.auc;
y(13,:) = stats.preGo.ThetaAlpha.strength.auc;
y(14,:) = stats.preGo.ThetaAlpha.SSIndex.auc;
y(15,:) = stats.preGo.ThetaAlpha.sourceIndex.auc;
y(16,:) = stats.preGo.ThetaAlpha.sinkIndex.auc;

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