clear all
close all
clc

filename = 'AT647N-AT542'; % AT647N-AT542 and AF594-AF488 the probes overlapping
ID = readtable(['datasets/' 'iFISH_single_colour_overlap_' filename '.csv']);
OrderPairs = {'a594'; 'tmr'; 'Cy5'; 'a488'};

% the label 0 are indistinguishable clusters
Position0 = find(ID.Label == 0);
ID(Position0, :) = [];

%Reading the number of the last field from the long extension
LastIndex = ID.File(end);
clusters = 1;
for label = 1:2
    %There are 4 labels to distinguish both alleles, the label 2 is analyzed
    %separately afterwards label 1
    PositionNew = find(ID.Label == label);
    IDnew = ID(PositionNew, :);
    
    for Fields = 1:LastIndex % Select field of view
        
        PositionFile = find(IDnew.File == Fields); % Index of all dots belonging to this field
        nuclei = unique(IDnew.Nuclei(PositionFile));  % All nuclei in the current field
        
        for n = 1:length(nuclei)
            
            PositionNuclei = PositionFile(eq(nuclei(n), IDnew.Nuclei(PositionFile)));
            channel = IDnew.Channel(PositionNuclei);
            
            %plotMAXchannel(channel, IDnew.Value(PositionNuclei), IDnew.x(PositionNuclei), IDnew.y(PositionNuclei), IDnew.z(PositionNuclei))
            
            X = nan(1, 4);
            Y = nan(1, 4);
            Z = nan(1, 4);
            
            for FlO = 1:length(OrderPairs)
                channelPicked = strsplit(OrderPairs{FlO});
                PositionChannels = PositionNuclei(strcmp(channelPicked{1}, IDnew.Channel(PositionNuclei)));
                
                DotsCount(clusters, FlO) = length(PositionChannels);
                
                if ~isempty(PositionChannels)
                    
                    [~, Index] = max(IDnew.Value(PositionChannels)); % Find the highest intensity per channel
                    
                    X(FlO) = IDnew.x(PositionChannels(Index))*0.130; % Get the coordinates already pixel converted to µm
                    Y(FlO) = IDnew.y(PositionChannels(Index))*0.130;
                    Z(FlO) = IDnew.z(PositionChannels(Index))*0.200;
                    
                    fwhm(clusters, FlO) = IDnew.FWHM(PositionChannels(Index));
                end
            end
            
            if strcmp(filename, 'AT647N-AT542')

                Dist(clusters) = pdist2([X(2), Y(2), Z(2)], [X(3), Y(3), Z(3)]);
                DistW1(clusters) = pdist2([X(1), Y(1), Z(1)], [X(3), Y(3), Z(3)]);
                DistW2(clusters) = pdist2([X(2), Y(2), Z(2)], [X(4), Y(4), Z(4)]);
                
            elseif strcmp(filename, 'AF594-AF488')
                
                Dist(clusters) = pdist2([X(1), Y(1), Z(1)], [X(4), Y(4), Z(4)]);
                DistW1(clusters) = pdist2([X(1), Y(1), Z(1)], [X(3), Y(3), Z(3)]);
                DistW2(clusters) = pdist2([X(2), Y(2), Z(2)], [X(4), Y(4), Z(4)]);
                
            end
            
            clusters = 1 + clusters;
        end
    end
end

Dist(isnan(Dist)) = [];
DistW1(isnan(DistW1)) = [];
DistW2(isnan(DistW2)) = [];
Overlap = Dist<0.5;
CorrectOverlapPercentageMax = sum(Overlap)/length(Dist)*100
IncorrectOverlapPercentageMax1 = sum(DistW1(DistW1<0.5))/length(DistW1)*100
IncorrectOverlapPercentageMax2 = sum(DistW2(DistW2<0.5))/length(DistW2)*100

figure(1)
subplot(2, 2, 1)
histogram(Dist, 30)
legend(num2str(length(Dist)))
title(['mean = ' num2str(round(mean(Dist), 2))])

if strcmp(filename, 'AT647N-AT542')
    xlabel('Cy5-tmr channels distance in same cluster (µm)')
elseif strcmp(filename, 'AF594-AF488')
    xlabel('a594-a488 channels distance in same cluster (µm)')
end

subplot(2, 2, 2)
histogram(Dist(Overlap), 20)
legend(num2str(length(Dist(Overlap))))
if strcmp(filename, 'AT647N-AT542')
    xlabel('Cy5-tmr channels distance below 0.5µm (µm)')
elseif strcmp(filename, 'AF594-AF488')
    xlabel('a594-a488 channels distance below 0.5µm (µm)')
end

subplot(2, 2, 3)
histogram(DistW1, 30)
legend(num2str(length(DistW1)))
xlabel('a594-Cy5 channels distance (µm)')
title(['mean = ' num2str(round(mean(DistW1), 2))])

subplot(2, 2, 4)
histogram(DistW2, 30)
legend(num2str(length(DistW2)))
xlabel('tmr-a488 channels distance(µm)')
title(['mean = ' num2str(round(mean(DistW2), 2))])




figure(2)
subplot(2, 2, 1)
histogram(DotsCount(:, 1))
xlabel('a594 channel dots per cluster')

subplot(2, 2, 2)
histogram(DotsCount(:, 2))
xlabel('tmr channel dots per cluster')

subplot(2, 2, 3)
histogram(DotsCount(:, 3))
xlabel('Cy5 channel dots per cluster')

subplot(2, 2, 4)
histogram(DotsCount(:, 4))
xlabel('a488 channel dots per cluster')




fwhm(fwhm < 0.05) = nan;

figure(3)
subplot(2, 2, 1)
histogram(fwhm(:, 1))
xlabel('a594 channel FWHM')
xlim([2 5])

subplot(2, 2, 2)
histogram(fwhm(:, 2))
xlabel('tmr channel FWHM')
xlim([2 5])

subplot(2, 2, 3)
histogram(fwhm(:, 3))
xlabel('Cy5 channel FWHM')
xlim([2 5])

subplot(2, 2, 4)
histogram(fwhm(:, 4))
xlabel('a488 channel FWHM')
xlim([2 5])