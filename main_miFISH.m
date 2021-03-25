clear all
close all
clc

ThresholdDetection = 0.55; % A channel combination is only considered when their distance is below 0.55µm

replicate = 1; % 1 or 2
ID = readtable(['datasets/' 'miFISH_chr2_rep' num2str(replicate) '.csv']); %337 or 377

OrderPairs = {'a488 a700'; 'a594  '; 'Cy5 a700'; 'a488 a594'; 'tmr  '; 'Cy5 a594'; 'tmr a488'; 'a700  '; ...
    'tmr Cy5'; 'a700 a594'; 'Cy5  '; 'tmr a700'; 'Cy5 a488'; 'ir800  '; 'a488  '; 'tmr a594'};

GenomicDistance = [0 20 40 60 63 66 69 72 75 78 81 84 87 90 93 113]; %[20 20 20 3 3 3 3 3 3 3 3 3 3 3 20];
DistGenomics = 0:3:115;
Contact = zeros(length(DistGenomics), 3);
ContactRatio = zeros(length(DistGenomics), 1);
selectedClusters = 0;

% the label 0 is indistinguishable clusters
Position0 = find(ID.Label == 0);
ID(Position0, :) = [];

% Reading the number of the last field from the long extension
LastIndex = ID.File(end);
clusters = 1;
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
            lam = nan(1, length(OrderPairs));
            
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
            
            Before(:, clusters) = sum(~isnan(IndexPart1), 2);
            
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
            
            lam(~isnan(IndexCluster)) = IDnew.lamin_dist_norm(IndexCluster(~isnan(IndexCluster)));
            Lamina(clusters, :) = lam;
            
            Coordinates(:, :, clusters) = [X' Y' Z']; % Coordinates XYZ
            
            NucleiCluster(clusters) = nuclei(n) + (Fields-1)*30; %giving a specific cellID but would be the same for different labels
           
            if any(sum(~isnan(Z))>14)
                selectedClusters = selectedClusters + 1;
                LaminaSelec(selectedClusters, :) = lam;
                Coordinates_15_16(:, :, selectedClusters) = [X' Y' Z'];
                %plot3D_splineCurve(Coordinates(:, :, clusters), OrderPairs, GenomicDistance);
            end
            
            clusters = clusters + 1;
        end
        
    end
end

sum(sum(~isnan(Coordinates(:, 1, :))))






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




function plot3D_splineCurve(Coordinates, OrderPairs, GenomicDistance, varargin)

ProbesPosition = readtable('probes_position.csv');
Index = floor(ProbesPosition.chromstart/10^5);
fileID = fopen('eigen_KR_2_100kb.txt','r');
EigenValues = fscanf(fileID,'%f')';
fclose(fileID);
EigenIndex = EigenValues(Index);

Compartments(:, EigenIndex > 0.01) = {'A'};
Compartments(:, EigenIndex < -0.005) = {'B'};
Compartments(:, EigenIndex >= -0.005 & EigenIndex <= 0.01) = {' '};


ToRemove = isnan(Coordinates(:, 1));
Coordinates(ToRemove, :) = [];
OrderPairs(ToRemove) = [];
Compartments(ToRemove) = [];
GenomicDistance(ToRemove) = [];

%% plot each dot in 3D coordinates
hfig = figure;
set(hfig, 'Position', [100, 100, 1200, 1500]);
subplot(1, 2, 1)

Cx = Coordinates(:, 1);
Cy = Coordinates(:, 2);
Cz = Coordinates(:, 3);
plot3(Cx, Cy, Cz, 'ro','LineWidth',2)
text(Cx, Cy, Cz, strcat({'  '}, OrderPairs))
axis equal
xlabel('X (µm)')
ylabel('Y (µm)')
zlabel('Z (µm)')
grid

%% make a spline curve from the points

subplot(1, 2, 2)
points = fnplt(cscvn([Cx, Cy, Cz]'));
lengthP = length(points);
colors_p = [linspace(0.1,0.9,lengthP)', linspace(0.1,0.9,lengthP)', linspace(0.9,0.1,lengthP)'];

for p = 1:2:lengthP-2
    plot3(points(1, p:p+2),points(2, p:p+2),points(3, p:p+2), 'color', colors_p(p, :), 'LineWidth', 1);
    hold on
end

plot3(Cx(strcmp(Compartments, 'A')), Cy(strcmp(Compartments, 'A')), Cz(strcmp(Compartments, 'A')), 'ro','MarkerSize', 5, 'LineWidth',2)
plot3(Cx(strcmp(Compartments, 'B')), Cy(strcmp(Compartments, 'B')), Cz(strcmp(Compartments, 'B')), 'bo','MarkerSize', 5,'LineWidth',2)
plot3(Cx(strcmp(Compartments, ' ')), Cy(strcmp(Compartments, ' ')), Cz(strcmp(Compartments, ' ')), 'go','MarkerSize', 5,'LineWidth',2)

text(Cx, Cy, Cz, strcat({'  '}, num2str(GenomicDistance')), 'FontSize',12, 'FontWeight','normal')

xlabel('X (µm)')
ylabel('Y (µm)')
zlabel('Z (µm)')
grid off
axis equal
set(gca,'FontSize',13,'FontWeight','normal', 'linewidth',1)
hold off
end