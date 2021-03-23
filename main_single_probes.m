clear all
close all
clc

filename = '378'; % 378 and 379
ID = readtable(['datasets/' 'iAM' filename '.csv']);
OrderPairs = {'a594'; 'tmr'; 'Cy5'; 'a488'};

% the label 0 are indistinguishable clusters
Position0 = find(ID.Label == 0);
ID(Position0, :) = [];

%Reading the number of the last field from the long extension
LastIndex = strsplit(ID.File{length(ID.File)},'/');
Extension = strjoin(LastIndex(1:length(LastIndex)-1), '/');
LastIndex = strsplit(LastIndex{length(LastIndex)},'.');

clusters = 1;
for label = 1:2
    %There are 4 labels to distinguish both alleles, the label 2 is analyzed
    %separately afterwards label 1
    PositionNew = find(ID.Label == label);
    IDnew = ID(PositionNew, :);
    
    for Fields = 1:str2double(LastIndex{1})
        
        if Fields < 10
            extension = [Extension, '/00', num2str(Fields), '.NM'];
        elseif Fields > 9 && Fields < 100
            extension = [Extension, '/0', num2str(Fields), '.NM'];
        else
            extension = [Extension, num2str(Fields), '.NM'];
        end
        
        PositionFile =  find(strcmp(extension, IDnew.File));
        nuclei = unique(IDnew.Nuclei(PositionFile));  % Select Nuclei
        
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
                    
                    X(FlO) = IDnew.x(PositionChannels(Index))*0.130; % Get the coordinates already pixel converted to �m
                    Y(FlO) = IDnew.y(PositionChannels(Index))*0.130;
                    Z(FlO) = IDnew.z(PositionChannels(Index))*0.200;
                    
                    fwhm(clusters, FlO) = IDnew.FWHM(PositionChannels(Index));
                end
            end
            
            if strcmp(filename, '378')

                Dist(clusters) = pdist2([X(2), Y(2), Z(2)], [X(3), Y(3), Z(3)]);
                DistW1(clusters) = pdist2([X(1), Y(1), Z(1)], [X(3), Y(3), Z(3)]);
                DistW2(clusters) = pdist2([X(2), Y(2), Z(2)], [X(4), Y(4), Z(4)]);
                
            elseif strcmp(filename, '379')
                
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

if strcmp(filename, '378')
    xlabel('Cy5-tmr channels distance in same cluster (�m)')
elseif strcmp(filename, '379')
    xlabel('a594-a488 channels distance in same cluster (�m)')
end

subplot(2, 2, 2)
histogram(Dist(Overlap), 20)
legend(num2str(length(Dist(Overlap))))
if strcmp(filename, '378')
    xlabel('Cy5-tmr channels distance below 0.5�m (�m)')
elseif strcmp(filename, '379')
    xlabel('a594-a488 channels distance below 0.5�m (�m)')
end

subplot(2, 2, 3)
histogram(DistW1, 30)
legend(num2str(length(DistW1)))
xlabel('a594-Cy5 channels distance (�m)')
title(['mean = ' num2str(round(mean(DistW1), 2))])

subplot(2, 2, 4)
histogram(DistW2, 30)
legend(num2str(length(DistW2)))
xlabel('tmr-a488 channels distance(�m)')
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
xaxis([2 5])

subplot(2, 2, 2)
histogram(fwhm(:, 2))
xlabel('tmr channel FWHM')
xaxis([2 5])

subplot(2, 2, 3)
histogram(fwhm(:, 3))
xlabel('Cy5 channel FWHM')
xaxis([2 5])

subplot(2, 2, 4)
histogram(fwhm(:, 4))
xlabel('a488 channel FWHM')
xaxis([2 5])