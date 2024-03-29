clear all
close all
clc

% chromosome 1, 2 or 10
chromosome = 2;

% '' '_rep1' '_rep2' 
% chromosome 1 does not have replicate
replicate = '_rep1'; 

ID = readtable(['datasets/' 'iFISH_chr' num2str(chromosome) replicate '.csv']);

if chromosome == 2 % chr2
    Colours = {'a700'; 'ir800'; 'a488'; 'tmr'; 'Cy5'; 'a594'};
    RealName = {'AF700'; 'AF790'; 'AF488'; 'AT542'; 'AT647N'; 'AF594'};
elseif chromosome == 10 % chr10
    Colours = {'a594'; 'Cy5'; 'tmr'; 'a488'; 'ir800'; 'a700'};
    RealName = {'AF594'; 'AT647N'; 'AT542'; 'AF488'; 'AF790'; 'AF700'};
elseif chromosome == 1 % chr1
    Colours = {'a594'; 'Cy5'; 'tmr'; 'a488'; 'ir800'; 'a700'};
    RealName = {'AF594'; 'AT647N'; 'AT542'; 'AF488'; 'AF790'; 'AF700'};
end

clusters = 1;
% only 10 dots are accepted
X = nan(10, length(Colours), 1000);
Y = nan(10, length(Colours), 1000);
Z = nan(10, length(Colours), 1000);

% the label 0 are indistinguishable clusters
Position0 = find(ID.Label == 0);
ID(Position0, :) = [];

%Reading the number of the last field from the long extension
LastIndex = ID.File(length(ID.File));

for label = 1:2
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
            
            for FlO = 1:length(Colours)
                channelPicked = Colours{FlO};
                
                PositionChannels = PositionNuclei(strcmp(channelPicked, IDnew.Channel(PositionNuclei)));
                
                %The brightest 9 dots per channel might be the single probe
                if chromosome == 2 && strcmp(channelPicked, 'a700') == 1
                    [~, Index] = maxk(IDnew.Value(PositionChannels), 10); % Find the highest intensity per channel
                    IndexDots = PositionChannels(Index);
                else
                    [~, Index] = maxk(IDnew.Value(PositionChannels), 4); % Find the highest intensity per channel
                    IndexDots = PositionChannels(Index);
                end
                
                X(1:length(IndexDots), FlO, clusters) = IDnew.x(IndexDots)*0.130; % Get the coordinates already pixel converted to µm
                Y(1:length(IndexDots), FlO, clusters) = IDnew.y(IndexDots)*0.130;
                Z(1:length(IndexDots), FlO, clusters) = IDnew.z(IndexDots)*0.200;
            end
            
            clusters = clusters + 1;
        end
        
    end
end
disp('Number of clusters: ')
disp(clusters-1)

figure
for n = 1:length(Colours)
    for Flo = 1:length(Colours)
        Dist = [];
        for c = 1:clusters-1
            DistTemp = pdist2([X(:, n, c), Y(:, n, c), Z(:, n, c)], [X(:, Flo, c), Y(:, Flo, c), Z(:, Flo, c)]);
            DistTemp(DistTemp == 0) = nan;
            DistTemp(DistTemp < 0.5) = nan;
            Dist = [Dist DistTemp(:)'];
        end
        Distance(Flo, :) = Dist;
        Pattern(Flo, :) = Flo*ones(1, length(Distance(Flo, :)));
    end
subplot(2, 3, n)
h = boxplot(Distance(:), Pattern(:), 'Labels', Colours);
xlabel('Probes')
ylabel(['Pairwise distances against ' Colours{n}])
ylim([0 4])
set(gca,'FontSize',12,'FontWeight','Bold', 'linewidth',1)
set(h,{'linew'},{4})
grid
end