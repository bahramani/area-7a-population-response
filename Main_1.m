% Amirreza Bahramani
% Advanced Neuroscience
% March 2023
% HW 2
% Single Unit Analysis & Population Response Structure

close all;
clear;
clc;

%% UPLOAD DATA
load('UnitsData.mat')

%% PARAMETERS
pre = 1200; % ms
post = 2000; % ms
winSizeMovMean = 100;
saveFig = true;
% figure
% set(gcf, 'WindowState', 'maximized')
% figName = 'Distribution of Different Biases';
% title(figName)
% if saveFig
%     saveas(gcf,['pics' figName '.png']);
% end

%% STER 1: PSTH
%% PSTH for All Units in every Condition
numUnit = length(Unit);
numCnd = length(Unit(1).Cnd);
tVec = -pre:1:post;
spikeMat = cell(1,numUnit);
dataPSTH = struct('cnd1', {}, 'cnd2', {}, 'cnd3', {}, 'cnd4', {}, ...
                  'cnd5', {}, 'cnd6', {}, 'all', {});

for i = 1 : numUnit

    tempSpikeTimes = cellfun(@(x) x*1000,Unit(i).Trls,'un',0);
    spikeMat{i} = zeros(numel(tempSpikeTimes),pre+1+post);
    
    for j = 1:numel(tempSpikeTimes)
        for k = 1:numel(tempSpikeTimes{j})
            spikeMat{i}(j,floor(tempSpikeTimes{j}(k)+pre+1)) = 1;
        end
    end

    dataPSTH(i).all  = movmean((sum(spikeMat{i})/numel(spikeMat{i}(:,1)))*1000, winSizeMovMean)';
    dataPSTH(i).cnd1 = movmean((sum(spikeMat{i}(Unit(i).Cnd(1).TrialIdx,:))/length(Unit(i).Cnd(1).TrialIdx))*1000, 2*winSizeMovMean)';
    dataPSTH(i).cnd2 = movmean((sum(spikeMat{i}(Unit(i).Cnd(2).TrialIdx,:))/length(Unit(i).Cnd(2).TrialIdx))*1000, 2*winSizeMovMean)';
    dataPSTH(i).cnd3 = movmean((sum(spikeMat{i}(Unit(i).Cnd(3).TrialIdx,:))/length(Unit(i).Cnd(3).TrialIdx))*1000, 2*winSizeMovMean)';
    dataPSTH(i).cnd4 = movmean((sum(spikeMat{i}(Unit(i).Cnd(4).TrialIdx,:))/length(Unit(i).Cnd(4).TrialIdx))*1000, 2*winSizeMovMean)';
    dataPSTH(i).cnd5 = movmean((sum(spikeMat{i}(Unit(i).Cnd(5).TrialIdx,:))/length(Unit(i).Cnd(5).TrialIdx))*1000, 2*winSizeMovMean)';
    dataPSTH(i).cnd6 = movmean((sum(spikeMat{i}(Unit(i).Cnd(6).TrialIdx,:))/length(Unit(i).Cnd(6).TrialIdx))*1000, 2*winSizeMovMean)';

end

%% PSTH for each unit for all Conditions
unitNum = 425; % choose between 1 to 481

figure('Name','Single Unit Response', 'Position', [400 300 1000 600])
% Raster Plot
subplot(3,1,1)
hold on
for i = 1:192
    plot(tVec(spikeMat{unitNum}(i,:)~=0),spikeMat{unitNum}(i,spikeMat{unitNum}(i,:)~=0)*i,'k.', ...
        'MarkerSize',2)
end
axis([-pre+10 post+10 0.5 192+0.5])
xline(0,'-','Cue Onset', 'Color','#77AC30', 'LineWidth',1)
xline(900,'-','Target Onset','Color',"#A2142F",'LineWidth',1)
xlabel('Time [ms]')
ylabel('Trial')
title('Raster Plot')

% PSTH Plot
subplot(3,1,[2 3])
plot(tVec, dataPSTH(unitNum).all,'k','LineWidth',1.5)
axis([-pre+10 post+10 min(dataPSTH(unitNum).all)-0.5 max(dataPSTH(unitNum).all)+0.5])
xline(0,'-','Cue Onset', 'Color','#77AC30', 'LineWidth',1)
xline(900,'-','Target Onset','Color',"#A2142F",'LineWidth',1)
xlabel('Time [ms]')
ylabel('Firing Rate [Hz]')
title('Peri-Stimulus Time Histogram')
grid on
sgtitle(['Single Unit Response of Neuron #', num2str(unitNum), ' for All Conditions'])

figName = ['Single Unit Response of Neuron #', num2str(unitNum), ' for All Conditions'];
if saveFig
    saveas(gcf,['pics/04_' figName '.png']);
end

%% PSTH for each Condition for a Neuron
% generate 9 random integers between 1 and 481
unitNums = sort(randi([1, 481], 1, 9));

figure
set(gcf, 'WindowState', 'maximized')

for i = 1:length(unitNums)
    k = unitNums(i);

    maxValue = max(structfun(@max, dataPSTH(k)));
    minValue = min(structfun(@min, dataPSTH(k)));

    subplot(3,3,i)
    hold on

    plot(tVec, dataPSTH(k).cnd1, 'Color', "#D95319", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).cnd2, 'Color', "#EDB120", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).cnd3, 'Color', "#7E2F8E", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).cnd4, 'Color', "#77AC30", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).cnd5, 'Color', "#4DBEEE", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).cnd6, 'Color', "#A2142F", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).all  ,'k', 'LineWidth', 1.5)

    axis([-pre+10 post+10 minValue-0.5 maxValue+0.5])
    xline(0,'-','Cue Onset', 'Color','#77AC30', 'LineWidth',2)
    xline(900,'-','Target Onset','Color',"#A2142F",'LineWidth',2)
    legend({'3, -1', '3, +1', '6, -1', '6, +1', '9, -1', '9, +1', 'All'}, 'Location','northwest')
    xlabel('Time [ms]')
    ylabel('Firing Rate [Hz]')
    title(['PSTH of Neuron #', num2str(k)])
    grid on

end
figName = 'PSTH for each Condition for a Neuron';
if saveFig
    saveas(gcf,['pics/05_' figName '.png']);
end

%% Mean PSTH of All
figure
set(gcf, 'WindowState', 'maximized')
tmp = mean([dataPSTH(:).all], 2);
plot(tVec, tmp,'k', 'LineWidth', 2)
axis([-pre+10 post+10 min(tmp)-0.5 max(tmp)+0.5])
xline(0,'-','Cue Onset', 'Color','#77AC30', 'LineWidth',2)
xline(900,'-','Target Onset','Color',"#A2142F",'LineWidth',2)
xlabel('Time [ms]')
ylabel('Firing Rate [Hz]')
title('Mean PSTH of All Neurons in All Cinditions')
grid minor

figName = 'Mean PSTH of All Neurons in All Cinditions';
if saveFig
    saveas(gcf,['pics/06_' figName '.png']);
end

%% Mean PSTH for each Condition
figure
set(gcf, 'WindowState', 'maximized')

subplot(3,2,1)
plot(tVec, mean([dataPSTH(:).cnd1], 2), 'Color', "#D95319", 'LineWidth', 2)
axis([-pre+10 post+10 min(mean([dataPSTH(:).cnd1], 2))-0.5 max(mean([dataPSTH(:).cnd1], 2))+0.5])
xline(0,'-','Cue Onset', 'Color','#77AC30', 'LineWidth',2)
xline(900,'-','Target Onset','Color',"#A2142F",'LineWidth',2)
xlabel('Time [ms]')
ylabel('Firing Rate [Hz]')
title('PSTH of All Neurons for [3, -1]')
grid minor

subplot(3,2,2)
plot(tVec, mean([dataPSTH(:).cnd2], 2), 'Color', "#EDB120", 'LineWidth', 2)
axis([-pre+10 post+10 min(mean([dataPSTH(:).cnd2], 2))-0.5 max(mean([dataPSTH(:).cnd2], 2))+0.5])
xline(0,'-','Cue Onset', 'Color','#77AC30', 'LineWidth',2)
xline(900,'-','Target Onset','Color',"#A2142F",'LineWidth',2)
xlabel('Time [ms]')
ylabel('Firing Rate [Hz]')
title('PSTH of All Neurons for [3, +1]')
grid minor

subplot(3,2,3)
plot(tVec, mean([dataPSTH(:).cnd3], 2), 'Color', "#7E2F8E", 'LineWidth', 2)
axis([-pre+10 post+10 min(mean([dataPSTH(:).cnd3], 2))-0.5 max(mean([dataPSTH(:).cnd3], 2))+0.5])
xline(0,'-','Cue Onset', 'Color','#77AC30', 'LineWidth',2)
xline(900,'-','Target Onset','Color',"#A2142F",'LineWidth',2)
xlabel('Time [ms]')
ylabel('Firing Rate [Hz]')
title('PSTH of All Neurons for [6, -1]')
grid minor

subplot(3,2,4)
plot(tVec, mean([dataPSTH(:).cnd4], 2), 'Color', "#77AC30", 'LineWidth', 2)
axis([-pre+10 post+10 min(mean([dataPSTH(:).cnd4], 2))-0.5 max(mean([dataPSTH(:).cnd4], 2))+0.5])
xline(0,'-','Cue Onset', 'Color','#77AC30', 'LineWidth',2)
xline(900,'-','Target Onset','Color',"#A2142F",'LineWidth',2)
xlabel('Time [ms]')
ylabel('Firing Rate [Hz]')
title('PSTH of All Neurons for [6, +1]')
grid minor

subplot(3,2,5)
plot(tVec, mean([dataPSTH(:).cnd5], 2), 'Color', "#4DBEEE", 'LineWidth', 2)
axis([-pre+10 post+10 min(mean([dataPSTH(:).cnd5], 2))-0.5 max(mean([dataPSTH(:).cnd5], 2))+0.5])
xline(0,'-','Cue Onset', 'Color','#77AC30', 'LineWidth',2)
xline(900,'-','Target Onset','Color',"#A2142F",'LineWidth',2)
xlabel('Time [ms]')
ylabel('Firing Rate [Hz]')
title('PSTH of All Neurons for [9, -1]')
grid minor

subplot(3,2,6)
plot(tVec, mean([dataPSTH(:).cnd6], 2), 'Color', "#A2142F", 'LineWidth', 2)
axis([-pre+10 post+10 min(mean([dataPSTH(:).cnd6], 2))-0.5 max(mean([dataPSTH(:).cnd6], 2))+0.5])
xline(0,'-','Cue Onset', 'Color','#77AC30', 'LineWidth',2)
xline(900,'-','Target Onset','Color',"#A2142F",'LineWidth',2)
xlabel('Time [ms]')
ylabel('Firing Rate [Hz]')
title('PSTH of All Neurons for [9, +1]')
grid minor

figName = 'Mean PSTH for each Condition';
if saveFig
    saveas(gcf,['pics/07_' figName '.png']);
end

%% PSTH for each location and reward
figure
set(gcf, 'WindowState', 'maximized')

subplot(3,2,1)
tmp = mean([dataPSTH(:).all], 2);
plot(tVec, tmp, 'k', 'LineWidth', 2)
axis([-pre+10 post+10 min(tmp)-0.5 max(tmp)+0.5])
xline(0,'-','Cue Onset', 'Color','#77AC30', 'LineWidth',2)
xline(900,'-','Target Onset','Color',"#A2142F",'LineWidth',2)
xlabel('Time [ms]')
ylabel('Firing Rate [Hz]')
title('Mean PSTH of All Neurons in All Cinditions')
grid minor

subplot(3,2,2)
tmp = mean([mean([dataPSTH(:).cnd1], 2) , mean([dataPSTH(:).cnd2], 2)], 2);
plot(tVec, tmp, 'Color', "#D95319", 'LineWidth', 2)
axis([-pre+10 post+10 min(tmp)-0.5 max(tmp)+0.5])
xline(0,'-','Cue Onset', 'Color','#77AC30', 'LineWidth',2)
xline(900,'-','Target Onset','Color',"#A2142F",'LineWidth',2)
xlabel('Time [ms]')
ylabel('Firing Rate [Hz]')
title('PSTH of Cue REV=3')
grid minor

subplot(3,2,4)
tmp = mean([mean([dataPSTH(:).cnd3], 2) , mean([dataPSTH(:).cnd4], 2)], 2);
plot(tVec, tmp, 'Color', "#EDB120", 'LineWidth', 2)
axis([-pre+10 post+10 min(tmp)-0.5 max(tmp)+0.5])
xline(0,'-','Cue Onset', 'Color','#77AC30', 'LineWidth',2)
xline(900,'-','Target Onset','Color',"#A2142F",'LineWidth',2)
xlabel('Time [ms]')
ylabel('Firing Rate [Hz]')
title('PSTH of REV=6')
grid minor

subplot(3,2,6)
tmp = mean([mean([dataPSTH(:).cnd5], 2) , mean([dataPSTH(:).cnd6], 2)], 2);
plot(tVec, tmp, 'Color', "#7E2F8E", 'LineWidth', 2)
axis([-pre+10 post+10 min(tmp)-0.5 max(tmp)+0.5])
xline(0,'-','Cue Onset', 'Color','#77AC30', 'LineWidth',2)
xline(900,'-','Target Onset','Color',"#A2142F",'LineWidth',2)
xlabel('Time [ms]')
ylabel('Firing Rate [Hz]')
title('PSTH of REV=9')
grid minor

subplot(3,2,3)
tmp = mean([mean([dataPSTH(:).cnd1], 2) , mean([dataPSTH(:).cnd3], 2), mean([dataPSTH(:).cnd5], 2)], 2);
plot(tVec, tmp, 'Color', "#77AC30", 'LineWidth', 2)
axis([-pre+10 post+10 min(tmp)-0.5 max(tmp)+0.5])
xline(0,'-','Cue Onset', 'Color','#77AC30', 'LineWidth',2)
xline(900,'-','Target Onset','Color',"#A2142F",'LineWidth',2)
xlabel('Time [ms]')
ylabel('Firing Rate [Hz]')
title('PSTH of Cue Location=-1')
grid minor

subplot(3,2,5)
tmp = mean([mean([dataPSTH(:).cnd2], 2) , mean([dataPSTH(:).cnd4], 2), mean([dataPSTH(:).cnd6], 2)], 2);
plot(tVec, tmp, 'Color', "#4DBEEE", 'LineWidth', 2)
axis([-pre+10 post+10 min(tmp)-0.5 max(tmp)+0.5])
xline(0,'-','Cue Onset', 'Color','#77AC30', 'LineWidth',2)
xline(900,'-','Target Onset','Color',"#A2142F",'LineWidth',2)
xlabel('Time [ms]')
ylabel('Firing Rate [Hz]')
title('PSTH of Cue Location=+1')
grid minor

figName = 'PSTH for each location and reward';
if saveFig
    saveas(gcf,['pics/08_' figName '.png']);
end

%% STEP 2
timeBin = 50; %ms
windowLength = 9;

pVals = zeros(3, numUnit);
X = cell(1, numUnit);
y = cell(1, numUnit);
for i = 1:length(y)
    y{i} = zeros(192,3);
end

for i = 1:numUnit
    for j = 1:192
        X{i}(j,:) = PSTH(Unit(i).Trls(j), windowLength, floor((pre+post+1)/timeBin), linspace(-pre/1000, post/1000, floor((pre+post+1)/timeBin)));
    end
    y{i}(Unit(i).Cnd(1).TrialIdx,:) = [3, -1, 1].*ones(size(y{i}(Unit(i).Cnd(1).TrialIdx,:)));
    y{i}(Unit(i).Cnd(2).TrialIdx,:) = [3, 1, 2].* ones(size(y{i}(Unit(i).Cnd(2).TrialIdx,:)));
    y{i}(Unit(i).Cnd(3).TrialIdx,:) = [6, -1, 3].*ones(size(y{i}(Unit(i).Cnd(3).TrialIdx,:)));
    y{i}(Unit(i).Cnd(4).TrialIdx,:) = [6, 1, 4].* ones(size(y{i}(Unit(i).Cnd(4).TrialIdx,:)));
    y{i}(Unit(i).Cnd(5).TrialIdx,:) = [9, -1, 5].*ones(size(y{i}(Unit(i).Cnd(5).TrialIdx,:)));
    y{i}(Unit(i).Cnd(6).TrialIdx,:) = [9, 1, 6].* ones(size(y{i}(Unit(i).Cnd(6).TrialIdx,:)));

    mdl = fitglm(X{i}, y{i}(:,1));
    pVals(1, i) = coefTest(mdl);

    mdl = fitglm(X{i}, y{i}(:,2));
    pVals(2, i) = coefTest(mdl);

    mdl = fitglm(X{i}, y{i}(:,3));
    pVals(3, i) = coefTest(mdl);

end


%%
bestREV = find(pVals(1,:) < 0.01);
bestLoc = find(pVals(2,:) < 0.01);
bestCnd = find(pVals(3,:) < 0.01);

figure
set(gcf, 'WindowState', 'maximized')
hold on
histogram(pVals(1,:), 50)
histogram(pVals(2,:), 50)
histogram(pVals(3,:), 50)
xlabel('p-value')
ylabel('Neuron Count')
title('Performance of GLM')
legend('Reward Expected Value', 'Cue Location', 'Both(Condition)')
figName = 'Performance of GLM';
if saveFig
    saveas(gcf,['pics/09_' figName '.png']);
end

%% BEST REV
unitNums = sort(bestREV(randperm(length(bestREV), 9)));

figure
set(gcf, 'WindowState', 'maximized')

for i = 1:length(unitNums)
    k = unitNums(i);

    maxValue = max(structfun(@max, dataPSTH(k)));
    minValue = min(structfun(@min, dataPSTH(k)));

    subplot(3,3,i)
    hold on

    plot(tVec, dataPSTH(k).cnd1, 'Color', "#D95319", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).cnd2, 'Color', "#EDB120", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).cnd3, 'Color', "#7E2F8E", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).cnd4, 'Color', "#77AC30", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).cnd5, 'Color', "#4DBEEE", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).cnd6, 'Color', "#A2142F", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).all  ,'k', 'LineWidth', 1.5)

    axis([-pre+10 post+10 minValue-0.5 maxValue+0.5])
    xline(0,'-','Cue Onset', 'Color','#77AC30', 'LineWidth',2)
    xline(900,'-','Target Onset','Color',"#A2142F",'LineWidth',2)
    legend({'3, -1', '3, +1', '6, -1', '6, +1', '9, -1', '9, +1', 'All'}, 'Location','northwest')
    xlabel('Time [ms]')
    ylabel('Firing Rate [Hz]')
    title(['PSTH of Neuron #', num2str(k)])
    grid on

end
sgtitle('PSTH for Best REV-encoding Units')
figName = 'PSTH for Best REV-encoding Units';
if saveFig
    saveas(gcf,['pics/10_' figName '.png']);
end

%% BEST CND
unitNums = sort(bestCnd(randperm(length(bestCnd), 9)));

figure
set(gcf, 'WindowState', 'maximized')

for i = 1:length(unitNums)
    k = unitNums(i);

    maxValue = max(structfun(@max, dataPSTH(k)));
    minValue = min(structfun(@min, dataPSTH(k)));

    subplot(3,3,i)
    hold on

    plot(tVec, dataPSTH(k).cnd1, 'Color', "#D95319", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).cnd2, 'Color', "#EDB120", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).cnd3, 'Color', "#7E2F8E", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).cnd4, 'Color', "#77AC30", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).cnd5, 'Color', "#4DBEEE", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).cnd6, 'Color', "#A2142F", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).all  ,'k', 'LineWidth', 1.5)

    axis([-pre+10 post+10 minValue-0.5 maxValue+0.5])
    xline(0,'-','Cue Onset', 'Color','#77AC30', 'LineWidth',2)
    xline(900,'-','Target Onset','Color',"#A2142F",'LineWidth',2)
    legend({'3, -1', '3, +1', '6, -1', '6, +1', '9, -1', '9, +1', 'All'}, 'Location','northwest')
    xlabel('Time [ms]')
    ylabel('Firing Rate [Hz]')
    title(['PSTH of Neuron #', num2str(k)])
    grid on

end
sgtitle('PSTH for Best Cnd-encoding Units')
figName = 'PSTH for Best Cnd-encoding Units';
if saveFig
    saveas(gcf,['pics/11_' figName '.png']);
end

%% BEST Loc
unitNums = sort(bestLoc(randperm(length(bestLoc), 9)));

figure
set(gcf, 'WindowState', 'maximized')

for i = 1:length(unitNums)
    k = unitNums(i);

    maxValue = max(structfun(@max, dataPSTH(k)));
    minValue = min(structfun(@min, dataPSTH(k)));

    subplot(3,3,i)
    hold on

    plot(tVec, dataPSTH(k).cnd1, 'Color', "#D95319", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).cnd2, 'Color', "#EDB120", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).cnd3, 'Color', "#7E2F8E", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).cnd4, 'Color', "#77AC30", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).cnd5, 'Color', "#4DBEEE", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).cnd6, 'Color', "#A2142F", 'LineWidth', 0.5)
    plot(tVec, dataPSTH(k).all  ,'k', 'LineWidth', 1.5)

    axis([-pre+10 post+10 minValue-0.5 maxValue+0.5])
    xline(0,'-','Cue Onset', 'Color','#77AC30', 'LineWidth',2)
    xline(900,'-','Target Onset','Color',"#A2142F",'LineWidth',2)
    legend({'3, -1', '3, +1', '6, -1', '6, +1', '9, -1', '9, +1', 'All'}, 'Location','northwest')
    xlabel('Time [ms]')
    ylabel('Firing Rate [Hz]')
    title(['PSTH of Neuron #', num2str(k)])
    grid on

end
sgtitle('PSTH for Best Location-encoding Units')
figName = 'PSTH for Best Location-encoding Units';
if saveFig
    saveas(gcf,['pics/12_' figName '.png']);
end

%% PART 3

PCs = cell(1,7);
[~, PCs{1}, ~] = pca([dataPSTH(:).cnd1]);
[~, PCs{2}, ~] = pca([dataPSTH(:).cnd2]);
[~, PCs{3}, ~] = pca([dataPSTH(:).cnd3]);
[~, PCs{4}, ~] = pca([dataPSTH(:).cnd4]);
[~, PCs{5}, ~] = pca([dataPSTH(:).cnd5]);
[~, PCs{6}, ~] = pca([dataPSTH(:).cnd6]);
[~, PCs{7}, ~] = pca([dataPSTH(:).all]);

colorMap = ["#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F", 'k'];

figure
set(gcf, 'WindowState', 'maximized')
for i = 1:6
    hold on
    plot(PCs{i}(:,1),PCs{i}(:,2), 'Color', colorMap(i), 'LineWidth', 1);
    scatter(PCs{i}(1200,1),PCs{i}(1200,2), 100, 'green', 'filled')
    scatter(PCs{i}(2100,1),PCs{i}(2100,2), 100, 'red', 'filled')
end
plot(PCs{7}(:,1),PCs{7}(:,2), 'Color', 'k', 'LineWidth', 2);
scatter(PCs{7}(1200,1),PCs{7}(1200,2), 100, 'green', 'filled')
scatter(PCs{7}(2100,1),PCs{7}(2100,2), 100, 'red', 'filled')

xlabel('Principal Component 1')
ylabel('Principal Component 2')
title('PCA Analysis of All Units in 2D')
grid minor
figName = 'PCA Analysis of All Units in 2D';
if saveFig
    saveas(gcf,['pics/13_' figName '.png']);
end

%% 3D PCA
figure
set(gcf, 'WindowState', 'maximized')
for i = 1:6
    hold on
    plot3(PCs{i}(:,1), PCs{i}(:,2), PCs{i}(:,3), 'Color', colorMap(i), 'LineWidth', 1);
    scatter3(PCs{i}(1200,1),PCs{i}(1200,2),PCs{i}(1200,3), 100, 'green', 'filled')
    scatter3(PCs{i}(2100,1),PCs{i}(2100,2),PCs{i}(2100,3), 100, 'red', 'filled')
    scatter3(PCs{i}(1,1),PCs{i}(1,2),PCs{i}(1,3), 100, 'black', 'filled')
    scatter3(PCs{i}(end,1),PCs{i}(end,2),PCs{i}(end,3), 100, 'yellow', 'filled')
end
plot3(PCs{7}(:,1), PCs{7}(:,2), PCs{7}(:,3), 'Color', 'k', 'LineWidth', 1);
scatter3(PCs{7}(1200,1),PCs{7}(1200,2),PCs{7}(1200,3), 100, 'green', 'filled')
scatter3(PCs{7}(2100,1),PCs{7}(2100,2),PCs{7}(2100,3), 100, 'red', 'filled')

xlabel('Principal Component 1')
ylabel('Principal Component 2')
zlabel('Principal Component 3')
title('PCA Analysis of All Units in 3D')
grid minor
figName = 'PCA Analysis of All Units in 3D';
if saveFig
    saveas(gcf,['pics/14_' figName '.png']);
end

%% PART 4 Shuffling

addpath("CFR-master\")

matPSTH = zeros(numUnit, numCnd, 3201);
for i = 1:numUnit

    matPSTH(i, 1, :) = dataPSTH(i).cnd1;
    matPSTH(i, 2, :) = dataPSTH(i).cnd2;
    matPSTH(i, 3, :) = dataPSTH(i).cnd3;
    matPSTH(i, 4, :) = dataPSTH(i).cnd4;
    matPSTH(i, 5, :) = dataPSTH(i).cnd5;
    matPSTH(i, 6, :) = dataPSTH(i).cnd6;

end

%%
surrData = CFR(matPSTH,6,10);

%%
% Shuffle neurons
perm_neurons = randperm(numUnit);
dataPSTH_sh_neuron = struct();
for i = 1:numUnit
    dataPSTH_sh_neuron(i).cnd1 = dataPSTH(perm_neurons(i)).cnd1;
    dataPSTH_sh_neuron(i).cnd2 = dataPSTH(perm_neurons(i)).cnd2;
    dataPSTH_sh_neuron(i).cnd3 = dataPSTH(perm_neurons(i)).cnd3;
    dataPSTH_sh_neuron(i).cnd4 = dataPSTH(perm_neurons(i)).cnd4;
    dataPSTH_sh_neuron(i).cnd5 = dataPSTH(perm_neurons(i)).cnd5;
    dataPSTH_sh_neuron(i).cnd6 = dataPSTH(perm_neurons(i)).cnd6;
    dataPSTH_sh_neuron(i).all = dataPSTH(perm_neurons(i)).all;
end

% Shuffle conditions
dataPSTH_sh_cnd = struct();
for i = 1:numUnit
    perm_conditions = randperm(numCnd);
    dataPSTH_sh_cnd(i).cnd1 = dataPSTH(i).(sprintf('cnd%d', perm_conditions(1)));
    dataPSTH_sh_cnd(i).cnd2 = dataPSTH(i).(sprintf('cnd%d', perm_conditions(2)));
    dataPSTH_sh_cnd(i).cnd3 = dataPSTH(i).(sprintf('cnd%d', perm_conditions(3)));
    dataPSTH_sh_cnd(i).cnd4 = dataPSTH(i).(sprintf('cnd%d', perm_conditions(4)));
    dataPSTH_sh_cnd(i).cnd5 = dataPSTH(i).(sprintf('cnd%d', perm_conditions(5)));
    dataPSTH_sh_cnd(i).cnd6 = dataPSTH(i).(sprintf('cnd%d', perm_conditions(6)));
    dataPSTH_sh_cnd(i).all = dataPSTH(i).all;
end

%%
PCs = cell(1,7);
[~, PCs{1}, ~] = pca([dataPSTH_sh_neuron(:).cnd1]);
[~, PCs{2}, ~] = pca([dataPSTH_sh_neuron(:).cnd2]);
[~, PCs{3}, ~] = pca([dataPSTH_sh_neuron(:).cnd3]);
[~, PCs{4}, ~] = pca([dataPSTH_sh_neuron(:).cnd4]);
[~, PCs{5}, ~] = pca([dataPSTH_sh_neuron(:).cnd5]);
[~, PCs{6}, ~] = pca([dataPSTH_sh_neuron(:).cnd6]);
[~, PCs{7}, ~] = pca([dataPSTH_sh_neuron(:).all]);

colorMap = ["#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F", 'k'];

figure
set(gcf, 'WindowState', 'maximized')
for i = 1:6
    hold on
    plot(PCs{i}(:,1),PCs{i}(:,2), 'Color', colorMap(i), 'LineWidth', 1);
    scatter(PCs{i}(1200,1),PCs{i}(1200,2), 100, 'green', 'filled')
    scatter(PCs{i}(2100,1),PCs{i}(2100,2), 100, 'red', 'filled')
end
plot(PCs{7}(:,1),PCs{7}(:,2), 'Color', 'k', 'LineWidth', 2);
scatter(PCs{7}(1200,1),PCs{7}(1200,2), 100, 'green', 'filled')
scatter(PCs{7}(2100,1),PCs{7}(2100,2), 100, 'red', 'filled')

xlabel('Principal Component 1')
ylabel('Principal Component 2')
title('PCA Analysis of All Units with Neurons Shuffled')
grid minor
figName = 'PCA Analysis of All Units with Neurons Shuffled';
if saveFig
    saveas(gcf,['pics/15_' figName '.png']);
end

%%
PCs = cell(1,7);
[~, PCs{1}, ~] = pca([dataPSTH_sh_cnd(:).cnd1]);
[~, PCs{2}, ~] = pca([dataPSTH_sh_cnd(:).cnd2]);
[~, PCs{3}, ~] = pca([dataPSTH_sh_cnd(:).cnd3]);
[~, PCs{4}, ~] = pca([dataPSTH_sh_cnd(:).cnd4]);
[~, PCs{5}, ~] = pca([dataPSTH_sh_cnd(:).cnd5]);
[~, PCs{6}, ~] = pca([dataPSTH_sh_cnd(:).cnd6]);
[~, PCs{7}, ~] = pca([dataPSTH_sh_cnd(:).all]);

colorMap = ["#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F", 'k'];

figure
set(gcf, 'WindowState', 'maximized')
for i = 1:6
    hold on
    plot(PCs{i}(:,1),PCs{i}(:,2), 'Color', colorMap(i), 'LineWidth', 1);
    scatter(PCs{i}(1200,1),PCs{i}(1200,2), 100, 'green', 'filled')
    scatter(PCs{i}(2100,1),PCs{i}(2100,2), 100, 'red', 'filled')
end
plot(PCs{7}(:,1),PCs{7}(:,2), 'Color', 'k', 'LineWidth', 2);
scatter(PCs{7}(1200,1),PCs{7}(1200,2), 100, 'green', 'filled')
scatter(PCs{7}(2100,1),PCs{7}(2100,2), 100, 'red', 'filled')

xlabel('Principal Component 1')
ylabel('Principal Component 2')
title('PCA Analysis of All Units with Conditions Shuffled')
grid minor
figName = 'PCA Analysis of All Units with Conditions Shuffled';
if saveFig
    saveas(gcf,['pics/16_' figName '.png']);
end

















