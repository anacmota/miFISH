clear all
close all
clc

clusters = 1;
ThresholdDetection = 0.55; % A channel combination is only considered when their distance is below 0.55µm

ID = readtable(['datasets/' 'wCentr.out.dilate5.iAM377.csv']);

OrderPairs = {'a488 a700'; 'a594  '; 'Cy5 a700'; 'a488 a594'; 'tmr  '; 'Cy5 a594'; 'tmr a488'; 'a700  '; ...
    'tmr Cy5'; 'a700 a594'; 'Cy5  '; 'tmr a700'; 'Cy5 a488'; 'ir800  '; 'a488  '; 'tmr a594'};

GenomicDistance = [0 20 40 60 63 66 69 72 75 78 81 84 87 90 93 113]; %[20 20 20 3 3 3 3 3 3 3 3 3 3 3 20];
selectedClusters = 0;

% the label 0 is indistinguishable clusters
Position0 = find(ID.Label == 0);
ID(Position0, :) = [];

% Reading the number of the last field from the long extension
LastIndex = ID.File(end);

for label = 1:4 % There are 4 labels to distinguish both alleles, even when in G2
    PositionNew = find(ID.Label == label);
    IDnew = ID(PositionNew, :);
    
    for Fields = 1:LastIndex % Select field of view
        
        PositionFile = find(IDnew.File == Fields); % Index of all dots belonging to this field
        nuclei = unique(IDnew.Nuclei(PositionFile));  % All nuclei in the current field
        
        for n = 1:length(nuclei) % Select Nuclei
            
            PositionNuclei = PositionFile(eq(nuclei(n), IDnew.Nuclei(PositionFile))); % Index of all dots belonging to this nuclei
            channel = IDnew.Channel(PositionNuclei); % Channel corresponding to nuclei index
            
            if length(unique(channel)) < 6
                continue % There is no point to continue if there aren't all channels
            end
            
            
            limit = 15;
            IndexPart1 = nan(length(OrderPairs), limit); % one of the channels of overlap plus singles
            IndexPart2 = nan(length(OrderPairs), limit); % the other channel of overlap
            overlapDist = nan(length(OrderPairs), limit); % overlap distance
            X = nan(1, length(OrderPairs)); % X axis
            Y = nan(1, length(OrderPairs)); % Y axis
            Z = nan(1, length(OrderPairs)); % Z axis
            
            for FlO = 1:length(OrderPairs) % Select channel pair
                channelPicked = strsplit(OrderPairs{FlO});
                
                if isempty(channelPicked{2}) % single dots
                    PositionChannels = PositionNuclei(strcmp(channelPicked{1}, channel));  % Index of all dots belonging to this channel
                    
                    %The brightest dots per channel might be the single probe
                    Index = find(mean(IDnew.Value(PositionChannels)) < IDnew.Value(PositionChannels)); % Find dots with higher intensity than average
                    IndexPart1(FlO, 1:length(Index)) = PositionChannels(Index);
                else
                    PositionChannels1 = PositionNuclei(strcmp(channelPicked{1}, channel)); % Index of all dots belonging to the channel 1
                    PositionChannels2 = PositionNuclei(strcmp(channelPicked{2}, channel)); % Index of all dots belonging to the channel 2
                    
                    if length(PositionChannels1) > limit % If there is more than 10 dots per channel in a cluster (rare)
                        [~, a] = maxk(IDnew.Value(PositionChannels1), limit); % Select the 10 brightest
                        PositionChannels1 = PositionChannels1(a);
                        
                    elseif length(PositionChannels2) > limit
                        [~, a] = maxk(IDnew.Value(PositionChannels2), limit);
                        PositionChannels2 = PositionChannels2(a);
                    end
                    
                    Dist = pdist2([IDnew.x(PositionChannels1)*0.130 IDnew.y(PositionChannels1)*0.130 IDnew.z(PositionChannels1)*0.200], ...
                        [IDnew.x(PositionChannels2)*0.130 IDnew.y(PositionChannels2)*0.130 IDnew.z(PositionChannels2)*0.200]);
                    
                    [A, B] = find(Dist < ThresholdDetection); % From all possible combinations, only below the threshold are taken
                    if length(A) > limit
                        [Dist, a] = mink(Dist(Dist < ThresholdDetection), limit);
                        A = A(a);
                        B = B(a);
                    end
                    IndexPart1(FlO, 1:length(A)) = PositionChannels1(A);
                    IndexPart2(FlO, 1:length(B)) = PositionChannels2(B);
                    overlapDist(FlO, 1:length(A)) = Dist(Dist < ThresholdDetection);
                end
            end
            
            Idx = correctCloudPnt(IndexPart1, IDnew); % Dots off the cloud by 5µm are removed
            if ~isempty(Idx)
                IndexPart1(Idx) = NaN;
                IndexPart2(Idx) = NaN;
                overlapDist(Idx) = NaN;
            end
            
            Position = [IDnew.x*0.130 IDnew.y*0.130 IDnew.z*0.20];
            IndexCluster = correctOverlap(overlapDist, IndexPart1, IndexPart2, IDnew.Value, Position, IDnew.Channel); % Correct the overlaps
             
            X(~isnan(IndexCluster)) = IDnew.x(IndexCluster(~isnan(IndexCluster)))*0.130; % Get the coordinates already converted to µm
            X = X - min(X) + 0.001; % Set a common reference for every cluster
            Y(~isnan(IndexCluster)) = IDnew.y(IndexCluster(~isnan(IndexCluster)))*0.130;
            Y = Y - min(Y) + 0.001;
            Z(~isnan(IndexCluster)) = IDnew.z(IndexCluster(~isnan(IndexCluster)))*0.200;
            Z = Z - min(Z) + 0.001;
            
            Cord = [X' Y' Z'];
            
            if any(sum(~isnan(Z))>15)
                selectedClusters = selectedClusters + 1;
                Coordinates(:, :, selectedClusters) = Cord; % Coordinates XYZ
                
                ToRemove = find(isnan(Cord(:, 1)));
                if isempty(ToRemove)
                    curves{selectedClusters} = fnplt(cscvn(Cord'));
                elseif ToRemove~=1 && ToRemove~=16
                    Cord(ToRemove, :) = (Cord(ToRemove-1, :)+Cord(ToRemove+1, :))./2;
                    curves{selectedClusters} = fnplt(cscvn(Cord'));
                else
                    selectedClusters = selectedClusters - 1;
                end
            end
            
            clusters = clusters + 1;
        end
    end
end

PeakGenomic = [];
PeakCurve = [];
for aa = 1:selectedClusters
    C1 = curves{aa};
    [k1(:, aa), t1(:, aa)] = curv_tors(C1(1, :),C1(2, :),C1(3, :));
    corr=1/k1(:, aa);
    PeakGenomic = [PeakGenomic find(corr~=0)];
    PeakCurve = [PeakCurve corr(corr~=0)];
end
% [~, idx] = sort(PeakGenomic);
% PeakGenomic = PeakGenomic(idx);
% PeakCurve = PeakCurve(idx);


Y = pdist(PeakGenomic', 'squaredeuclidean');
Z = linkage(Y,'average');
figure
dendrogram(Z, selectedClusters, 'Orientation','left')
xlabel('Height')
c = cophenet(Z,Y)


[idx,C] = kmeans([PeakGenomic' PeakCurve'], 4,'Distance','sqeuclidean',...
    'MaxIter',10000, 'Display','final');
figure
gscatter(PeakGenomic,PeakCurve,idx,'bgmr')
hold on
plot(C(:,1),C(:,2),'kx', 'MarkerSize', 8, 'linewidth',2)
set(gca, 'YScale', 'log')
legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster Centroid')
xlabel('Genomic distance (Mb)')
ylabel('1/curvature (µm)')
for aa = 1:selectedClusters
text(PeakGenomic(aa), PeakCurve(aa), ['  ' num2str(aa)], 'FontSize', 8)
end
set(gca, 'FontSize', 8)
xlim([0 130])


peakK = [];
figure
for aa = 1:selectedClusters
    [p, loc] = max(abs(k1(:, aa)));
    peakK(aa, 1:length(p)) = loc;
    
    if aa < 38
        subplot(6, 7, aa)
        plot(1/k1(:, aa),'k', 'linewidth',1)
        title(['    ' num2str(aa)], 'FontSize', 12)
%         xlabel('Genomic distance (Mb)')
%         ylabel('1/curvature (µm)')
    end
end
set(gca, 'FontSize', 12)
sgtitle('curvature')

%print -painters -dsvg output.svg


function [k, t] = curv_tors(x,y,z)

dx = gradient(x);
dy = gradient(y);
dz = gradient(z);

ddx = gradient(dx);
ddy = gradient(dy);
ddz = gradient(dz);

dddx = gradient(ddx);
dddy = gradient(ddy);
dddz = gradient(ddz);

alpha = (ddz.*dy - ddy.*dz).^2 + ...
     (ddx.*dz - ddz.*dx).^2 + ...
     (ddy.*dx - ddx.*dy).^2;

k = (alpha).^(1/2)./ ...
     (dx.^2+dy.^2+dz.^2).^(3/2);
 
 t = (dddx.*(ddz.*dy - ddy.*dz) + ...
      dddy.*(ddx.*dz - ddz.*dx) + ...
      dddz.*(ddy.*dx - ddx.*dy)) ... 
     ./alpha;

end

function Idx = correctCloudPnt(IndexPart1, IDnew)
IndexPart1new = IndexPart1(~isnan(IndexPart1));
X = IDnew.x(IndexPart1new)*0.13;
Y = IDnew.y(IndexPart1new)*0.13;
Z = IDnew.z(IndexPart1new)*0.2;
ptCloud = pointCloud([X Y Z]);
%plot3(X, Y, Z, 'ro')
[~,inlierIndices,outlierIndices] = pcdenoise(ptCloud);

for n = 1:length(outlierIndices)
    Dist = sqrt((X(outlierIndices(n)) - X(inlierIndices)) .^ 2 + (Y(outlierIndices(n)) - Y(inlierIndices)) .^ 2 + (Z(outlierIndices(n)) - Z(inlierIndices)) .^ 2);
    Remove(n) = any(Dist < 5);
end
outlierIndices(Remove) = [];
Idx = find(ismember(IndexPart1, IndexPart1new(outlierIndices)));
end