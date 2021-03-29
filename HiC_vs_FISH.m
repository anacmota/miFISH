clearvars -except Lamina Coordinates

ProbesPosition = readtable('probes_position.csv');
Index = ProbesPosition.chromstart;
fileID = readtable(['datasets/' 'probes.HiC.csv']);

HiC_score = fileID.balanced;
HiC_probe2 = fileID.start_probe2;
HiC_probe1 = fileID.start_probe1;

PDist_0_3 = nan(length(HiC_probe2), length(Coordinates(1,1,:)));
PDist_0_6 = nan(length(HiC_probe2), length(Coordinates(1,1,:)));
PDist_0_9 = nan(length(HiC_probe2), length(Coordinates(1,1,:)));
PDist_1_5 = nan(length(HiC_probe2), length(Coordinates(1,1,:)));
PDist_2 = nan(length(HiC_probe2), length(Coordinates(1,1,:)));
PDist_3 = nan(length(HiC_probe2), length(Coordinates(1,1,:)));
PairwiseLam = [];
PairwiseDist = [];
for kk = 1:length(HiC_probe2)
    for n = 1:length(Coordinates(1,1,:))
        PairwiseLam(kk, n) = abs(Lamina(n, Index == HiC_probe2(kk))-Lamina(n, Index == HiC_probe1(kk)));
        PairwiseDist(kk, n) = pdist2(Coordinates(Index == HiC_probe2(kk), :, n), Coordinates(Index == HiC_probe1(kk), :, n));
        
        if PairwiseDist(kk, n) < 0.3 && PairwiseDist(kk, n)~=0
            PDist_0_3(kk, n) = 1;
        elseif PairwiseDist(kk, n) >= 0.3
            PDist_0_3(kk, n) = 0;
        end
        
        if PairwiseDist(kk, n) < 0.6 && PairwiseDist(kk, n)~=0
            PDist_0_6(kk, n) = 1;
        elseif PairwiseDist(kk, n) >= 0.6
            PDist_0_6(kk, n) = 0;
        end
        
        if PairwiseDist(kk, n) < 0.9 && PairwiseDist(kk, n)~=0
            PDist_0_9(kk, n) = 1;
        elseif PairwiseDist(kk, n) >= 0.9
            PDist_0_9(kk, n) = 0;
        end
        
        if PairwiseDist(kk, n) < 1.5 && PairwiseDist(kk, n)~=0
            PDist_1_5(kk, n) = 1;
        elseif PairwiseDist(kk, n) >= 1.5
            PDist_1_5(kk, n) = 0;
        end
        
        if PairwiseDist(kk, n) < 2 && PairwiseDist(kk, n)~=0
            PDist_2(kk, n) = 1;
        elseif PairwiseDist(kk, n) >= 2
            PDist_2(kk, n) = 0;
        end
        
        if PairwiseDist(kk, n) < 3 && PairwiseDist(kk, n)~=0
            PDist_3(kk, n) = 1;
        elseif PairwiseDist(kk, n) >= 3
            PDist_3(kk, n) = 0;
        end
    end
end
PDist_0_3(PairwiseDist == 0) = nan;
PDist_0_6(PairwiseDist == 0) = nan;
PDist_0_9(PairwiseDist == 0) = nan;
PDist_1_5(PairwiseDist == 0) = nan;
PDist_2(PairwiseDist == 0) = nan;
PDist_3(PairwiseDist == 0) = nan;
PairwiseLam(PairwiseLam == 0) = nan;
PairwiseDist(PairwiseDist == 0) = nan;

PairwiseDist_av = nanmedian(PairwiseDist, 2);
PDist_0_3_av = nansum(PDist_0_3, 2)/n;
PDist_0_6_av = nansum(PDist_0_6, 2)/n;
PDist_0_9_av = nansum(PDist_0_9, 2)/n;
PDist_1_5_av = nansum(PDist_1_5, 2)/n;
PDist_2_av = nansum(PDist_2, 2)/n;
PDist_3_av = nansum(PDist_3, 2)/n;

PDist_0_3_av(all(isnan(PDist_0_3), 2)) = nan;
PDist_0_6_av(all(isnan(PDist_0_6), 2)) = nan;
PDist_0_9_av(all(isnan(PDist_0_9), 2)) = nan;
PDist_1_5_av(all(isnan(PDist_1_5), 2)) = nan;
PDist_2_av(all(isnan(PDist_2), 2)) = nan;
PDist_3_av(all(isnan(PDist_3), 2)) = nan;

%Variables = [nanmedian(PairwiseDist, 2), nansum(PDist_0_3, 2)/n, nansum(PDist_0_6, 2)/n, nansum(PDist_0_9, 2)/n, nansum(PDist_1_5, 2)/n, nansum(PDist_2, 2)/n, nansum(PDist_3, 2)/n];
Variables = [PairwiseDist_av, PDist_0_3_av, PDist_0_6_av, PDist_0_9_av, PDist_1_5_av, PDist_2_av, PDist_3_av];

for t = 1:size(Variables, 2)
figure
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


HiC_score(isnan(PairwiseDist_av)) = [];
PDist_0_9_av(isnan(PairwiseDist_av)) = [];
PairwiseLam(isnan(PairwiseDist_av), :) = [];
PairwiseDist_av(isnan(PairwiseDist_av)) = [];


Probes_coord = table(HiC_score, PairwiseDist_av, PDist_0_9_av, nanmedian(PairwiseLam, 2));
writetable(Probes_coord, 'HiC_PairwiseDist_LamDist_15_16.txt', 'WriteVariableNames',0)

function [RHOscc, RHOpcc] = regression_plot_log(x, y, lam)
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