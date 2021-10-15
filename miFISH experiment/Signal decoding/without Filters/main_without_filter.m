clear all
close all
clc

% There are the replicate 1 or 2
replicate = 1; 
ID = readtable(['datasets/' 'miFISH_chr2_rep' num2str(replicate) '.csv']); %337 or 377

%All the possible color combinations used in miFISH are stated together and
%placed according to their genomic order
OrderPairs = {'a488 a700'; 'a594  '; 'Cy5 a700'; 'a488 a594'; 'tmr  '; 'Cy5 a594'; 'tmr a488'; 'a700  '; ...
    'tmr Cy5'; 'a700 a594'; 'Cy5  '; 'tmr a700'; 'Cy5 a488'; 'ir800  '; 'a488  '; 'tmr a594'};

% Genomic coordinates in Mb
GenomicDistance = [0 20 40 60 63 66 69 72 75 78 81 84 87 90 93 113];

% A channel combination is only considered when their distance is below 0.55µm
ThresholdDetection = 0.55; 

% the label 0 belongs to indistinguishable clusters. Removed from this
% analysis
Position0 = find(ID.Label == 0);
ID(Position0, :) = [];

clusters = 1;
selectedClusters = 0;

for label = 1:4 
    % There are 4 labels available to distinguish both alleles
    % Particularly relevant for G2/mitosis cells (not used in this
    % analysis)
    PositionNew = find(ID.Label == label);
    IDnew = ID(PositionNew, :);
    
    for Fields = 1:ID.File(end) 
        % Select field of view
        
        % Index of all dots belonging to this field
        PositionFile = find(IDnew.File == Fields); 
        % All nuclei in the current field
        nuclei = unique(IDnew.Nuclei(PositionFile));
        
        for n = 1:length(nuclei) 
            % Select Nuclei
            
            % Index of all dots belonging to this nuclei
            PositionNuclei = PositionFile(eq(nuclei(n), IDnew.Nuclei(PositionFile))); 
            % Channel corresponding to nuclei index
            channel = IDnew.Channel(PositionNuclei); 
            
            if length(unique(channel)) < 6
                continue 
                % The analysis does not proceed with incomplete clusters
            end
            
            limit = 10;
            % index of the dots present in the cluster
            IndexCluster = nan(length(OrderPairs), 1);
            
            % X axis coordinate
            X = nan(1, length(OrderPairs)); 
            % Y axis coordinate
            Y = nan(1, length(OrderPairs));
            % Z axis coordinate
            Z = nan(1, length(OrderPairs));
            % lamina distance
            lam = nan(1, length(OrderPairs));
            
            for FlO = 1:length(OrderPairs) 
                % Select channel pair
                
                %Splits the combinatory channels
                channelPicked = strsplit(OrderPairs{FlO});
                
                if isempty(channelPicked{2}) % single-color probes
                    PositionChannels = PositionNuclei(strcmp(channelPicked{1}, channel));  % Index of all dots belonging to this channel
                    
                    %The brightest dots per channel might be the single probe
                    Index = find(max(IDnew.Value(PositionChannels)));
                    IndexCluster(FlO) = PositionChannels(Index);
                
                else % dual-color probes
                    PositionChannels1 = PositionNuclei(strcmp(channelPicked{1}, channel)); % Index of all dots belonging to the channel 1
                    PositionChannels2 = PositionNuclei(strcmp(channelPicked{2}, channel)); % Index of all dots belonging to the channel 2
                    
                    if length(PositionChannels1) > limit 
                        % If there is more than 10 dots per channel in a cluster (rare)
                        
                        % Select the 10 brightest
                        [~, a] = maxk(IDnew.Value(PositionChannels1), limit); 
                        PositionChannels1 = PositionChannels1(a);
                        
                    elseif length(PositionChannels2) > limit
                        [~, a] = maxk(IDnew.Value(PositionChannels2), limit);
                        PositionChannels2 = PositionChannels2(a);
                    end
                    
                    % Compute pairwise distance of the 10 dots selected
                    Dist = pdist2([IDnew.x(PositionChannels1)*0.130 IDnew.y(PositionChannels1)*0.130 IDnew.z(PositionChannels1)*0.200], ...
                        [IDnew.x(PositionChannels2)*0.130 IDnew.y(PositionChannels2)*0.130 IDnew.z(PositionChannels2)*0.200]);
                    
                    % Find the closest dot
                    [A, B] = find(Dist == min(Dist,[],'all'));
                    
                    IndexPart1 = PositionChannels1(A);
                    IndexPart2 = PositionChannels2(B);
                    
                    if length(IndexPart1) > 1
                        % In case there are more than one pair at the same
                        % distance, it is chosen the brighest dot
                        [~, B] = max(IDnew.Value(IndexPart1));
                        IndexPart1 = IndexPart1(B);
                    end
                    
                    % From all possible combinations, only below the threshold are taken
                    if min(Dist,[],'all') > ThresholdDetection
                        IndexCluster(FlO) = NaN;
                    elseif min(Dist,[],'all') < ThresholdDetection
                        IndexCluster(FlO) = IndexPart1;
                    end
                end
                
            end
            
            
            Idx = correctCloudPnt(IndexPart1, IDnew); 
            % Dots off the cloud by 5µm are removed
            if ~isempty(Idx)
                IndexPart1(Idx) = NaN;
                IndexPart2(Idx) = NaN;
                overlapDist(Idx) = NaN;
            end
            
            Position = [IDnew.x*0.130 IDnew.y*0.130 IDnew.z*0.20];
            IndexCluster = correctOverlap(overlapDist, IndexPart1, IndexPart2, IDnew.Value, Position, IDnew.Channel); % Correct the overlaps
             
            X(~isnan(IndexCluster)) = IDnew.x(IndexCluster(~isnan(IndexCluster)))*0.130; % Get the coordinates already converted to µm
            % Set a common reference for every cluster
            X = X - min(X) + 0.001; 
            Y(~isnan(IndexCluster)) = IDnew.y(IndexCluster(~isnan(IndexCluster)))*0.130;
            Y = Y - min(Y) + 0.001;
            Z(~isnan(IndexCluster)) = IDnew.z(IndexCluster(~isnan(IndexCluster)))*0.200;
            Z = Z - min(Z) + 0.001;
            
            % Lamina distance normalised
            lam(~isnan(IndexCluster)) = IDnew.lamin_dist_norm(IndexCluster(~isnan(IndexCluster)));
            Lamina(clusters, :) = lam;
            
            % Coordinates XYZ
            Coordinates(:, :, clusters) = [X' Y' Z']; 
            
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
% Creates a cloud around all the dots and identifies the inliers and
% outliers to the is cloud of dots
% If the outliers are positioned at more than 5 µm away, the outlier dots are removed
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
% Plots all the dots assigned to the 16 probes
% Each probe is represented in different colors according to the
% compartment that the eigen value associates it to: red for A, blue for B
% and green for none
% A spline curve is draw to connect the 16 probes

ProbesPosition = readtable(['datasets/' 'probes_position.csv']);
Index = floor(ProbesPosition.chromstart/10^5);
fileID = fopen(['datasets/' 'eigen_values_HiC.txt'],'r');
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