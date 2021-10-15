function laminaDISTeigen(Coordinates, Lamina)

% Restore parameters used in the main script
% All the possible color combinations used in miFISH are stated together and
%placed according to their genomic order
OrderPairs = {'a488 a700'; 'a594  '; 'Cy5 a700'; 'a488 a594'; 'tmr  '; 'Cy5 a594'; 'tmr a488'; 'a700  '; ...
    'tmr Cy5'; 'a700 a594'; 'Cy5  '; 'tmr a700'; 'Cy5 a488'; 'ir800  '; 'a488  '; 'tmr a594'}; 

% Genomic distance from the consecutive probe
DistanceGenomic = [20 20 20 3 3 3 3 3 3 3 3 3 3 3 20];

% Read probes exact coordinates and eigen values from HiC datasets
ProbesPosition = readtable(['datasets/' 'probes_position.csv']);
Index = floor(ProbesPosition.chromstart/10^5);
fileID = fopen(['datasets/' 'eigen_values_HiC.txt'],'r');
EigenValues = fscanf(fileID,'%f')';
fclose(fileID);
EigenIndex = EigenValues(Index);

% Definition of compartments A and B
Compartments(:, EigenIndex(1:end) > 0.01) = {'A'};
Compartments(:, EigenIndex(1:end) < -0.005) = {'B'};
Compartments(:, EigenIndex(1:end) >= -0.005 & EigenIndex(1:end) <= 0.01) = {' '};


% pairwise distances of probes coordinates
PDISTall = [];
ProbesDistance = [];
for n = 1:length(Coordinates(1, 1, :))
    PDISTall(n, :) = pdist(Coordinates(:, :, n));
    ProbesDistance(:, :, n) = squareform(pdist(Coordinates(:, :, n)));
end 
% median distance for all possible combinations of probes
medAll = nanmedian(PDISTall, 1);


% Create vectors for cumulative genomic distance, Lamina distance and
% compartments associations (A-A) (A-B) (B-B)
DistanceGenomicCumulative = 0;
LamDif = [];
Comp = [];
for i = 1:length(OrderPairs)-1
    DistanceGenomicCumulative(i+1) = DistanceGenomicCumulative(end) + DistanceGenomic(i);
    LamDif = [LamDif Lamina(:, i) - Lamina(:, i+1:end)];
    Comp = [Comp, strcat(Compartments(i), Compartments(i+1:end))];
    interDistance(:, i) = ProbesDistance(i+1, i+0, :);
end
% pairwise distances of genomic distances
AllGENOMICdistances = pdist(DistanceGenomicCumulative');

% median of lamina position normalised
LamMedian = nanmedian(Lamina', 2);


%% All pairwise distances combination of probes with the corresponding eigenvector combination
hfig = figure;
set(hfig, 'Position', [0, 0, 1400, 800]);
subplot(2, 4, 1)
scatter(AllGENOMICdistances(strcmp(Comp, 'AA')), medAll(strcmp(Comp, 'AA')), 'filled')
hold on
scatter(AllGENOMICdistances(strcmp(Comp, 'BB')), medAll(strcmp(Comp, 'BB')), 'filled')
scatter(AllGENOMICdistances(strcmp(Comp, 'AB')|strcmp(Comp, 'BA')), medAll(strcmp(Comp, 'AB')|strcmp(Comp, 'BA')), 'filled')
xlabel('Genomic distance in Mbp')
ylabel('Observed distance in µm')
legend('AA', 'BB', 'AB')
grid on


%% All pairwise distances combination of probes (120 in total) represented according to the lamina distance difference
subplot(2, 4, 2)
count = nanmedian(LamDif, 1);
scatter(AllGENOMICdistances, medAll, 30, abs(count), 'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
xlabel('Genomic distance in Mbp')
ylabel('Observed distance in µm')
set(gca,'FontSize',12,'FontWeight','Bold', 'linewidth',1)
title('linear scale')
grid on

subplot(2, 4, 3)
scatter(AllGENOMICdistances, medAll, 30, abs(count), 'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
xlabel('Genomic distance in Mbp')
ylabel('Observed distance in µm')
title('logarithmic scale')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
h = colorbar;
ylabel(h, 'Lamina distance difference norm', 'FontSize', 12,'FontWeight','Bold')
set(h,'linew',2)
grid on


%% Representation of the median of lamina distance normalised for each probe and the color code is for the compartment classification
subplot(2, 4, 4)
scatter(DistanceGenomicCumulative(strcmp(Compartments, 'A')), nanmedian(Lamina(:, strcmp(Compartments, 'A'))), 'b', 'filled');
hold on
scatter(DistanceGenomicCumulative(strcmp(Compartments, 'B')), nanmedian(Lamina(:, strcmp(Compartments, 'B'))), 'r', 'filled');
scatter(DistanceGenomicCumulative(strcmp(Compartments, ' ')), nanmedian(Lamina(:, strcmp(Compartments, ' '))), 'g', 'filled');
xlabel('Genomic distance in Mbp')
ylabel('Lamina norm distance')
legend('A', 'B', 'NaN', 'Location','northwest')
grid

%% Boxplot representation for each eigenvector combination
subplot(2, 4, 5)
A = PDISTall(:, strcmp(Comp, 'AA'));
B = PDISTall(:, strcmp(Comp, 'BB'));
AB = PDISTall(:, strcmp(Comp, 'AB')|strcmp(Comp, 'BA'));
x = [ones(length(A(:)), 1); 2*ones(length(B(:)), 1); 3*ones(length(AB(:)), 1)];
scatter(x(:), [A(:); B(:); AB(:)],'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
h = boxplot([A(:); B(:); AB(:)], x, 'Labels', {'AA','BB', 'AB'});
ylabel('Observed distance in µm')
set(h,'linew',4)
grid on

%% Taking the 120 probes combinations, measure the median and make a boxplot for each eigenvector combination
subplot(2, 4, 6)
A = medAll(strcmp(Comp, 'AA'));
B = medAll(strcmp(Comp, 'BB'));
AB = medAll(strcmp(Comp, 'AB')|strcmp(Comp, 'BA'));
x = [ones(length(A(:)), 1); 2*ones(length(B(:)), 1); 3*ones(length(AB(:)), 1)];
color = [repmat([0 0 1], length(A(:)), 1); repmat([1 0 0], length(B(:)), 1); repmat([0.9294 0.6902 0.1294], length(AB(:)), 1)];
h = boxplot([A(:); B(:); AB(:)], x, 'Labels', {'AA','BB', 'AB'});
ylabel('Observed distance in µm')
set(h,{'linew'},{4})
hold on
scatter(x(:), [A(:); B(:); AB(:)], 30, color, 'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
grid on


%% Boxplot for all lamina distances measure for each probe
subplot(2, 4, 7:8)
h = boxplot(Lamina, 'positions', DistanceGenomicCumulative, 'labels', DistanceGenomicCumulative);
xlabel('Genomic distance in Mb')
ylabel('Lamina distance normalised')
grid on
set(gca,'FontSize',12,'FontWeight','Bold', 'linewidth',1)
set(h,'linew',4)


%% Representation of the Eigenvector distribution for all chr2, the probes are positioned at their corresponding genomic coordinate and the number indexed is the median lamina distance
figure
subplot(3, 1, 1)
plot(EigenValues)
hold on
scatter(Index, EigenIndex, 'filled')
text(Index, EigenIndex, num2str(round(LamMedian, 2)), 'FontWeight',  'bold')
ylabel('Eigen vector')
grid on

subplot(3, 1, 2)
plot(EigenValues)
hold on
scatter(Index, EigenIndex, 'filled')
axis([1550 2500 0.01 0.035])
text(Index(EigenIndex > 0.01), EigenIndex(EigenIndex > 0.01), num2str(round(LamMedian(EigenIndex > 0.01), 2)), 'FontWeight', 'bold')
ylabel('Eigen vector')
grid on

subplot(3, 1, 3)
plot(EigenValues)
hold on
scatter(Index, EigenIndex, 'filled')
axis([1550 2500 -0.02 -0.005])
text(Index(EigenIndex < -0.005), EigenIndex(EigenIndex < -0.005), num2str(round(LamMedian(EigenIndex < -0.005), 2)), 'FontWeight', 'bold')
xlabel('Genomic distance in 100kb')
ylabel('Eigen vector')
grid on    


%% all 3 and 20 distances grouped and show the behaviour of each 3Mbp and 20 Mbp
hfig = figure;
set(hfig, 'Position', [100, 100, 1000, 400]);
subplot(1, 3, 1)
x = DistanceGenomic.*ones(length(interDistance), 1);
x(x==3) = 1;
x(x==20) = 2;
scatter(x(:), interDistance(:),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
h = boxplot(interDistance, DistanceGenomic);
set(h,{'linew'},{4})
%ylim([0 10])
xlabel('Genomic distance in Mbp')
ylabel('Observed distance in µm')
hold off
set(gca,'FontSize',12,'FontWeight','Bold')
grid on
hold off

subplot(1, 3, 2)
P3 = (1:sum(sum(x==1)>0)).*ones(size(x, 1), 1);
scatter(P3(:), interDistance(x==1),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15)
hold on
h = boxplot(interDistance(x==1), P3(:));
set(h,{'linew'},{4})
xlabel('# probe')
ylabel('Observed distance in µm')
title('3 Mbp probes')
hold off
set(gca,'FontSize',12,'FontWeight','Bold')
grid on
hold off

subplot(1, 3, 3)
P20 = (1:sum(sum(x==2)>0)).*ones(size(x, 1), 1);
scatter(P20(:), interDistance(x==2),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15)
hold on
h = boxplot(interDistance(x==2), P20(:));
set(h,{'linew'},{4})
xlabel('# probe')
ylabel('Observed distance in µm')
title('20 Mbp probes')
hold off
set(gca,'FontSize',12,'FontWeight','Bold')
grid on
hold off
end