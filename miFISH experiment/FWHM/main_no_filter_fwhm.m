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
            FWHM1 = nan(limit, 6);
            
            n = 1;
            for FlO = 1:length(OrderPairs)
                channelPicked = strsplit(OrderPairs{FlO});
                
                if isempty(channelPicked{2})==1
                    PositionChannels = PositionNuclei(strcmp(channelPicked{1}, IDnew.Channel(PositionNuclei)));

                    FWHM1(1:length(PositionChannels), n) = IDnew.FWHM(PositionChannels);
                    n = n + 1;
                end
                
            end

            FWHM = [FWHM; FWHM1];
        end
        
    end
end

subplot(2, 3, 1)
histogram(FWHM(:, 1))
xlabel('a594 fwhm')
title(['median = ' num2str(round(nanmedian(FWHM(:, 1)), 2))])

subplot(2, 3, 2)
histogram(FWHM(:, 2))
xlabel('tmr fwhm')
title(['median = ' num2str(round(nanmedian(FWHM(:, 2)), 2))])

subplot(2, 3, 3)
histogram(FWHM(:, 3))
xlabel('a700 fwhm')
title(['median = ' num2str(round(nanmedian(FWHM(:, 3)), 2))])

subplot(2, 3, 4)
histogram(FWHM(:, 4))
xlabel('Cy5 fwhm')
title(['median = ' num2str(round(nanmedian(FWHM(:, 4)), 2))])

subplot(2, 3, 5)
histogram(FWHM(:, 5))
xlabel('ir800 fwhm')
title(['median = ' num2str(round(nanmedian(FWHM(:, 5)), 2))])

subplot(2, 3, 6)
histogram(FWHM(:, 6))
xlabel('a488 fwhm')
title(['median = ' num2str(round(nanmedian(FWHM(:, 6)), 2))])