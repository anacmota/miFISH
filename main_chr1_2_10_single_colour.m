clear all
close all
clc

filename = 370;
ID = readtable(['datasets/' 'wCentr.out.dilate5.iAM' num2str(filename) '.csv']);

if filename == 370 || filename == 454 % chr2
    Colours = {'a700'; 'ir800'; 'a488'; 'tmr'; 'Cy5'; 'a594'};
elseif filename == 371 || filename == 455 % chr10
    Colours = {'a594'; 'Cy5'; 'tmr'; 'a488'; 'ir800'; 'a700'};
elseif filename == 456 % chr1
    Colours = {'a594'; 'Cy5'; 'tmr'; 'a488'; 'ir800'; 'a700'};
end

clusters = 1;
X = nan(10, length(Colours), 100); %%%%%%%% only 10 are accepted
Y = nan(10, length(Colours), 100);
Z = nan(10, length(Colours), 100);
Lam = nan(10, length(Colours), 100);

% the label 0 are indistinguishable clusters
Position0 = find(ID.Label == 0);
ID(Position0, :) = [];

LastIndex = ID.File(end);
for label = 1:2
    PositionNew = find(ID.Label == label);
    IDnew = ID(PositionNew, :);
    
    for Fields = 1:LastIndex
        
        PositionFile =  find(IDnew.File == Fields);
        nuclei = unique(IDnew.Nuclei(PositionFile));  % Select Nuclei
        
        ForNucleiSurface = 0;
        for n = 1:length(nuclei)
            
            PositionNuclei = PositionFile(eq(nuclei(n), IDnew.Nuclei(PositionFile)));
            channel = IDnew.Channel(PositionNuclei);
            
            if length(unique(channel)) < 6
                continue
            end
            
            for FlO = 1:length(Colours)
                channelPicked = Colours{FlO};
                
                PositionChannels = PositionNuclei(strcmp(channelPicked, IDnew.Channel(PositionNuclei)));
                
                %The brightest 9 dots per channel might be the single probe
                if filename == 370 && strcmp(channelPicked, 'a700') == 1
                    [~, Index] = maxk(IDnew.Value(PositionChannels), 10); % Find the highest intensity per channel
                    IndexDots = PositionChannels(Index);
                else
                    [~, Index] = maxk(IDnew.Value(PositionChannels), 4); % Find the highest intensity per channel
                    IndexDots = PositionChannels(Index);
                end
                
                X(1:length(IndexDots), FlO, clusters) = IDnew.x(IndexDots)*0.130; % Get the coordinates already pixel converted to µm
                Y(1:length(IndexDots), FlO, clusters) = IDnew.y(IndexDots)*0.130;
                Z(1:length(IndexDots), FlO, clusters) = IDnew.z(IndexDots)*0.200;
                Lam(1:length(IndexDots), FlO, clusters) = IDnew.lamin_dist_norm(IndexDots);
            end
            
            cell_ID = IDnew.cell_ID(IndexDots(~isnan(IndexDots)));
            clusters = clusters + 1;
        end
        
    end
end
disp(clusters-1)

figure
for n = 1:length(Colours)
    for Flo = 1:length(Colours)
        Dist = [];
        for c = 1:clusters-1
            DistTemp = pdist2([X(:, n, c), Y(:, n, c), Z(:, n, c)], [X(:, Flo, c), Y(:, Flo, c), Z(:, Flo, c)])./Radius(c);
            DistTemp(DistTemp == 0) = nan;
            Dist = [Dist DistTemp(:)'];
        end
        Distance(Flo, :) = Dist;
        Pattern(Flo, :) = Flo*ones(1, length(Distance(Flo, :)));
    end
subplot(2, 3, n)
h = boxplot(Distance(:), Pattern(:), 'Labels', Colours);
xlabel('Probes')
ylabel(['Pairwise distances against ' Colours{n}])
ylim([0 1])
set(gca,'FontSize',12,'FontWeight','Bold', 'linewidth',1)
set(h,{'linew'},{4})
grid
end

Pattern = [];
for Flo = 1:length(Colours)  
    for c = 1:clusters-1
        LamMedian(Flo, c) = nanmedian(Lam(:, Flo, c));
        DistMedian(Flo, c) = nanmedian(pdist([X(:, Flo, c), Y(:, Flo, c), Z(:, Flo, c)]))./Radius(c);
    end
    Pattern(Flo, :) = Flo*ones(1, length(LamMedian(Flo, :)));
end
figure
subplot(1, 3, 1)
h = boxplot(LamMedian(:), Pattern(:), 'Labels', Colours);
xlabel('Probes')
ylabel('Lamina distance')
ylim([0 1])
set(gca,'FontSize',12,'FontWeight','Bold', 'linewidth',1)
set(h,{'linew'},{4})
grid

subplot(1, 3, 2)
h = boxplot(DistMedian(:), Pattern(:), 'Labels', Colours);
xlabel('Probes')
ylabel('Intra-probe distance')
ylim([0 1])
set(gca,'FontSize',12,'FontWeight','Bold', 'linewidth',1)
set(h,{'linew'},{4})
grid

subplot(1, 3, 3)
[RHOscc, RHOpcc] = regression_plot(DistMedian(:), LamMedian(:));
text(0.8, 0.8, sprintf(['PCC = ' num2str(RHOpcc) '\nSCC = ' num2str(RHOscc)]), 'FontSize', 8)
xlabel('Intra-probe distance')
ylabel('Lamina distance')
set(gca,'FontSize',12,'FontWeight','Bold', 'linewidth',1)