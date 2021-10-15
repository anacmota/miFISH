clearvars -except Lamina Coordinates

% Read the probes coordinates
ProbesPosition = readtable(['datasets/' 'probes_position.csv']);
Index = ProbesPosition.chromstart;
% Read the HiC data
fileID = readtable(['datasets/' 'probes.HiC.csv']);

HiC_score = fileID.balanced;
HiC_probe2 = fileID.start_probe2;
HiC_probe1 = fileID.start_probe1;

% Create variable for contacts in between probes at distances less than 1µm
PDist_1_0 = nan(l, length(Coordinates(1,1,:)));

PairwiseLam = [];
PairwiseDist = [];
for kk = 1:length(HiC_probe2)
    for n = 1:length(Coordinates(1,1,:))
        % Lamina distance in between the probes
        PairwiseLam(kk, n) = abs(Lamina(n, Index == HiC_probe2(kk))-Lamina(n, Index == HiC_probe1(kk)));
        % Pairwise distances in between the selected probe 2 against probe
        % 1 for each cluster
        PairwiseDist(kk, n) = pdist2(Coordinates(Index == HiC_probe2(kk), :, n), Coordinates(Index == HiC_probe1(kk), :, n));

        % Pairwise distances equal to 0 means that probe 1 and 2 are the
        %same
        % Binary system to evaluate the number of contacts
        if PairwiseDist(kk, n) < 1 && PairwiseDist(kk, n)~=0
            PDist_1_0(kk, n) = 1;
        elseif PairwiseDist(kk, n) >= 1
            PDist_1_0(kk, n) = 0;
        end
    end
end
% Pairwise distances equal to 0 means that probe 1 and 2 are the same
PDist_1_0(PairwiseDist == 0) = nan;
PairwiseLam(PairwiseLam == 0) = nan;
PairwiseDist(PairwiseDist == 0) = nan;

% Calculate the median of all contacts
PairwiseDist_av = nanmedian(PairwiseDist, 2);

% Calculate the ratio of contacts below 1 µm
PDist_1_0_av = nansum(PDist_1_0, 2)/n;
PDist_1_0_av(all(isnan(PDist_1_0), 2)) = nan;

Variables = [PairwiseDist_av, PDist_1_0_av];
for t = 1:size(Variables, 2)
    % Run the loop for the variables defined
    % Plot HiC against the miFISH contacts in linear and log scale
    
figure

% Linear scale
subplot(1, 2, 1)
scatter(HiC_score, Variables(:, t), 30, nanmedian(PairwiseLam, 2) , 'filled')
xlabel('HiC contact frequency')
if t == 1
    ylabel('Pairwise distances (µm)')
else
    ylabel('FISH contact frequency')
end
NewHiC = HiC_score;
NewHiC(isnan(Variables(:, t))) = [];
xlim([0 max(NewHiC)+0.01])
title('Linear scale')
grid

% Logarithmic scale plot and regression line
subplot(1, 2, 2)
[RHOscc, RHOpcc] = regression_plot_log(HiC_score, Variables(:, t), nanmedian(PairwiseLam, 2));

text(5*10^(-4), min(Variables(:, t))+(max(Variables(:, t))-min(Variables(:, t)))*0.1, sprintf(['PCC = ' num2str(RHOpcc) '\nSCC = ' num2str(RHOscc)]), 'FontSize', 8)

xlabel('HiC contact frequency')
if t == 1
    ylabel('Pairwise distances (µm)')
else
    ylabel('FISH contact frequency')
end
title('Log scale')
end


function [RHOscc, RHOpcc] = regression_plot_log(x, y, lam)
% Calculate the Pearson and Spearman correlations in between HiC and miFISH
% Scatter plot and linear regression
% Only performed for log scale

if any(isnan(x))
    y(isnan(x)) = [];
    lam(isnan(x)) = [];
    x(isnan(x)) = [];
end

if any(isnan(y))
    x(isnan(y)) = [];
    lam(isnan(y)) = [];
    y(isnan(y)) = [];
end

X=log10(x); Y=log10(y); 

[RHOscc,~] = corr(X,Y,'Type','Spearman');
[RHOpcc,~] = corr(X,Y,'Type','Pearson');

%Use polyfit to compute a linear regression that predicts y from x:
p = polyfit(X,Y,1);
if RHOscc < 0
    yhat=10.^polyval(p,[max(X) min(X)]);  % evaluate at end points
    loglog([max(x) min(x)],yhat) % add fitted line
else
    yhat=10.^polyval(p,[min(X) max(X)]);  % evaluate at end points
    loglog([min(x) max(x)],yhat) % add fitted line
end
hold on
scatter(x, y, 30, lam, 'filled','MarkerFaceAlpha',0.6)
h = colorbar;
ylabel(h, 'lamina difference')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
grid
end