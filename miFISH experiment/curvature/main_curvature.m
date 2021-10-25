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

for label = 1:4 % There are 4 labels to distinguish both alleles, even when in G2
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
            
            lam(~isnan(IndexCluster)) = IDnew.lamin_dist_norm(IndexCluster(~isnan(IndexCluster)));
            Lamina(clusters, :) = lam;
            
            if any(sum(~isnan(Z))>15)
                selectedClusters = selectedClusters + 1;
                Coordinates(:, :, selectedClusters) = Cord; % Coordinates XYZ
                
                % Interpolate the curve out of the 16 points
                curves{selectedClusters} = fnplt(cscvn(Cord'));
                LaminaInterp{selectedClusters} = interp1(GenomicDistance, lam, linspace(0,GenomicDistance(end),length(curves{selectedClusters})));
            end
            
            clusters = clusters + 1;
        end
    end
end



GenomeX = linspace(0,GenomicDistance(end),length(curves{1}));
for aa = 1:selectedClusters
    % Extract the curves 3D coordinates for each cluster
    C1 = curves{aa};
    % Calculate the curvature for each cluster
    k1(:, aa) = curv_tors(C1(1, :),C1(2, :),C1(3, :));
end


%% Represent the normalized curvature peaks
figure
peakK = [];
for aa = 1:selectedClusters
    [p, loc] = max(abs(k1(:, aa))); %1/k1(:, aa)
    peakK(aa, 1:length(p)) = loc;
    
    if aa < 38
        subplot(6, 7, aa)
        k = k1(:, aa)/max(k1(:, aa));
        k(k < 0.1) = 0;
        plot(GenomeX, k, 'Linewidth', 1)
        hold on
        title(['    ' num2str(aa)], 'FontSize', 12)
    end
end
set(gca, 'FontSize', 12)
sgtitle('Curvature normalised')

%% Scatter plot of all peaks > 0.1 from the 37 clusters
figure
for aa = 1:selectedClusters
    if aa < 38
        k = k1(:, aa)/max(k1(:, aa));
        k(k < 0.1) = 0;
        
        plot(GenomeX, k, 'k', 'Linewidth', 0.8)
        hold on
        
        P = find(k>0.1);
        Diff = [];
        for nn = 1:length(P)
        Diff(nn) = max(LaminaInterp{aa}(P(nn)-4:P(nn)+4))-min(LaminaInterp{aa}(P(nn)-4:P(nn)+4));
        end
        
        scatter(GenomeX(k>0.1), k(k>0.1), 60, Diff, 'filled')
    end
end
xlabel('Genomic distance (Mb)')
ylabel('Curvature peaks normalised')
xlim([60 93])
ylim([0 1.05])
h = colorbar;
ylabel(h,'Lamina difference within +/-3Mb')
set(gca, 'FontSize', 12)






function k = curv_tors(x,y,z)

dx = gradient(x);
dy = gradient(y);
dz = gradient(z);

ddx = gradient(dx);
ddy = gradient(dy);
ddz = gradient(dz);

alpha = (ddz.*dy - ddy.*dz).^2 + ...
     (ddx.*dz - ddz.*dx).^2 + ...
     (ddy.*dx - ddx.*dy).^2;

k = (alpha).^(1/2)./ ...
     (dx.^2+dy.^2+dz.^2).^(3/2);
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