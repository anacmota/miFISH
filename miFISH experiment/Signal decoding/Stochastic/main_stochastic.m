close all
clear all
clc

% There are the replicate 1 or 2
replicate = 1; 
ID = readtable(['datasets/' 'miFISH_chr2_rep' num2str(replicate) '.csv']); %337 or 377

%All the possible color combinations used in miFISH are stated together and
%placed according to their genomic order
OrderPairs = {'a488 a700'; 'a594  '; 'Cy5 a700'; 'a488 a594'; 'tmr  '; 'Cy5 a594'; 'tmr a488'; 'a700  '; ...
    'tmr Cy5'; 'a700 a594'; 'Cy5  '; 'tmr a700'; 'Cy5 a488'; 'ir800  '; 'a488  '; 'tmr a594'};

Channels = {'a488',  'tmr', 'a594', 'Cy5',  'a700', 'ir800'};

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
    PositionNew = find(ID.Label == label);
    IDnew = ID(PositionNew, :);
    
    for Fields = 1:ID.File(end) % Select field of view
        
        PositionFile = find(IDnew.File == Fields); % Index of all dots belonging to this field
        nuclei = unique(IDnew.Nuclei(PositionFile));  % All nuclei in the current field
        
        for n = 1:length(nuclei) % Select Nuclei
            
            PositionNuclei = PositionFile(eq(nuclei(n), IDnew.Nuclei(PositionFile))); % Index of all dots belonging to this nuclei
            channelAll = IDnew.Channel(PositionNuclei); % Channel corresponding to nuclei index
            
            if length(unique(channelAll)) < 6
                continue % There is no point to continue if there aren't all channels
            end
            
            
            Dots = []; % [c, x, y, z, i]
            c = 1; % color
            for FlO = 1:length(Channels) % Select channel pair
                channelPicked = Channels{FlO};
                
                PositionChannels = PositionNuclei(strcmp(channelPicked, channelAll));  % Index of all dots belonging to this channel
             
                x = IDnew.x(PositionChannels)*0.130;
                y = IDnew.y(PositionChannels)*0.130;
                z = IDnew.z(PositionChannels)*0.200;
                lam = IDnew.lamin_dist_norm(PositionChannels);
                i = IDnew.Value(PositionChannels);
                i = i-min(i(:));
                if max(i(:))>0
                    i = i./max(i(:));
                end
                
                Dots = [Dots; c*ones(length(x), 1), x, y, z, i, lam];
                
                c = c + 1;
            end
            [Coordinates(:, :, clusters), Lamina(clusters, :)] = probecombination(Dots);
            
            clusters = clusters + 1;
        end
    end
end

%% How many dots are in a set of 16 probes
figure
Coor(:, :) = Coordinates(:, 1, :);
[ii,~,kk] = unique(sum(~isnan(Coor))); % ii are the unique channels and kk allows the recreation of the variable like ii(kk)
freq=accumarray(kk,1); % 1 means that all cells have the same value and kk is counted

X = categorical(ii);
bar(X, freq) 
xlabel('# probes per cluster')
ylabel('# Clusters')
