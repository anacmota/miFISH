function IndexFinal = correctOverlap(overlapDist, IndexPart1, IndexPart2, Intensities, Position, Channel)
%% If the overlap is less than 0.25µm means that the pairing is highly trusted

TrustedDotsIndex = find(overlapDist < 0.25); % Find all pairs below 0.25
[~, idx] = sort(overlapDist(TrustedDotsIndex)); % Even within those below 0.25, the shortest would be the best
TrustedDotsIndex(idx) = TrustedDotsIndex;

IndexInCheck = IndexPart1(TrustedDotsIndex); % Dot index from part1

        
for check = 1:length(TrustedDotsIndex)
    if IndexPart1(TrustedDotsIndex(check)) == IndexInCheck(check)
        [IndexPart1, IndexPart2, overlapDist] = removeRepeat(IndexInCheck(check), IndexPart1, IndexPart2, TrustedDotsIndex(check), overlapDist);
        [IndexPart1, IndexPart2, overlapDist] = mergeChannel(IndexInCheck(check), Position, Channel, IndexPart1, IndexPart2, overlapDist);
        
    end
end

% Same is done for the part2
TrustedDotsIndex = find(overlapDist < 0.25); % Find all pairs below 0.25
[~, idx] = sort(overlapDist(TrustedDotsIndex)); % Even within those below 0.25, the shortest would be the best
TrustedDotsIndex(idx) = TrustedDotsIndex;
IndexInCheck = IndexPart2(TrustedDotsIndex);
for check = 1:length(TrustedDotsIndex)
    if IndexPart2(TrustedDotsIndex(check)) == IndexInCheck(check)
    [IndexPart2, IndexPart1, overlapDist] = removeRepeat(IndexInCheck(check), IndexPart2, IndexPart1, TrustedDotsIndex(check), overlapDist);
    [IndexPart2, IndexPart1, overlapDist] = mergeChannel(IndexInCheck(check), Position, Channel, IndexPart2, IndexPart1, overlapDist);
    end
end
overlapDist(isnan(IndexPart1)) = NaN;


%% Give the only chance of overlap for certain pairs the priority

RowInCheck = (sum(~isnan(IndexPart1), 2) == 1); % Rows where only one combination was detected

minInCheck = overlapDist(RowInCheck, :);
IndexInCheck = IndexPart1(RowInCheck, :);
minInCheck = minInCheck';
IndexInCheck = IndexInCheck';
minInCheck = minInCheck(~isnan(IndexInCheck));
IndexInCheck = IndexInCheck(~isnan(IndexInCheck)); % Index present in the referred row
[minInCheck,a] = sort(minInCheck(:, 1));
IndexInCheck = IndexInCheck(a, :);

for n = 1:size(IndexInCheck, 1)
    TrustedDotsIndex = find(overlapDist == minInCheck(n));
    
    if IndexPart1(TrustedDotsIndex) == IndexInCheck(n)
        [IndexPart1, IndexPart2, overlapDist] = removeRepeat(IndexInCheck(n), IndexPart1, IndexPart2, TrustedDotsIndex, overlapDist);
        [IndexPart1, IndexPart2, overlapDist] = mergeChannel(IndexInCheck(n), Position, Channel, IndexPart1, IndexPart2, overlapDist);
    end
end

% Same is done for the part2

RowInCheck = (sum(~isnan(IndexPart2), 2) == 1);

minInCheck = overlapDist(RowInCheck, :);
IndexInCheck = IndexPart2(RowInCheck, :);
minInCheck = minInCheck';
IndexInCheck = IndexInCheck';
minInCheck = minInCheck(~isnan(IndexInCheck));
IndexInCheck = IndexInCheck(~isnan(IndexInCheck)); % Index present in the referred row
[minInCheck,a] = sort(minInCheck(:, 1));
IndexInCheck = IndexInCheck(a, :);

for n = 1:size(IndexInCheck, 1)
    TrustedDotsIndex = find(overlapDist == minInCheck(n));
    if IndexPart2(TrustedDotsIndex) == IndexInCheck(n)
        [IndexPart2, IndexPart1, overlapDist] = removeRepeat(IndexInCheck(n), IndexPart2, IndexPart1, TrustedDotsIndex, overlapDist);
        [IndexPart2, IndexPart1, overlapDist] = mergeChannel(IndexInCheck(n), Position, Channel, IndexPart2, IndexPart1, overlapDist);
    end
end

overlapDist(isnan(IndexPart1)) = NaN; % Remove distance values from discarded combinations


%% Pick the lowest distance overlaps per probe unless index is repeated
Test{1} = IndexPart1;
Test{2} = IndexPart2;

RowInCheck = (sum(~isnan(overlapDist), 2) > 1);
BlockDist = overlapDist(RowInCheck, :);

for kk = 1:2
    minValues = nanmin(BlockDist, [], 2); % From remaining combinations, look for the minimum distance at each row
    % Check if the lowest distance value per row is not repeated
    minIndex = Test{kk}(ismember(overlapDist, minValues)); % minimum value per row is associated to the dot index
    
    [a, ~, c]=unique(minIndex);
    a_counts = accumarray(c,1); % count if there is repetition of a dot index
    
    while all(a_counts == 1) == 0 % in case is found repetition in at least one case
        RepeatedIndex = a(a_counts > 1); % Take all repeated indexes
        for n = 1:length(RepeatedIndex)
            RepeatInd = ismember(Test{kk}(RowInCheck, :), RepeatedIndex(n));
            minDist = min(BlockDist(RepeatInd)); % finds the minimum distance
            
            if kk == 1
                [IndexPart1, IndexPart2,overlapDist] = removeRepeat(RepeatedIndex(n), IndexPart1, IndexPart2, find(eq(minDist, overlapDist)), overlapDist);
                Test{1} = IndexPart1;
                Test{2} = IndexPart2;
            elseif kk == 2
                [IndexPart2, IndexPart1, overlapDist] = removeRepeat(RepeatedIndex(n), IndexPart2, IndexPart1, find(eq(minDist, overlapDist)), overlapDist);
                Test{1} = IndexPart1;
                Test{2} = IndexPart2;
            end
            
            RepeatInd(sum(ismember(BlockDist, minDist), 2)>0, :) = zeros(1, size(Test{kk}, 2));
            A = sort(BlockDist(sum(RepeatInd, 2)>0, :));
            
            i = 0;
            B = 0;
            while (B == 0 && i+1 < length(A))
                i = i + 1;
                B = any(sum(ismember(overlapDist, A(1+i))));
            end
            if B>0
                a(eq(RepeatedIndex(n), a)) = Test{kk}(ismember(overlapDist, A(1+i)));
            else
                a(eq(RepeatedIndex(n), a)) = [];
            end
        end
        [a, ~, c]=unique(a);
        a_counts = accumarray(c,1);
    end
end

overlapDist(isnan(IndexPart1)) = NaN;

% Now the rows are replaced with the minimum distance only

for n = 1:size(IndexPart1, 1)
    M = nanmin(overlapDist(n, :)); % The index value and the position in the row
    if ~isnan(M)
        filter = zeros(size(overlapDist));
        filter(n, :) = ones(1, size(overlapDist, 2));
        TrustedDotsIndex = find(eq(overlapDist, M) & filter);
        TrustedDotsIndex = TrustedDotsIndex(1);
        
        Index = IndexPart1(TrustedDotsIndex);
        if ~isnan(Index)
            [IndexPart1, IndexPart2,overlapDist] = removeRepeat(Index, IndexPart1, IndexPart2, TrustedDotsIndex, overlapDist);
        end

        TrustedDotsIndex = find(eq(overlapDist, M) & filter);
        Index = IndexPart2(TrustedDotsIndex);
        if ~isnan(Index)
            [IndexPart2, IndexPart1, overlapDist] = removeRepeat(Index, IndexPart2, IndexPart1, TrustedDotsIndex, overlapDist);
        end
    end
end





%% The single probes should be chosen by the maximum intensity
Doubles = IndexPart1(sum(~isnan(IndexPart1), 2)==1, :);
for n = 1:size(Doubles, 1)
    [IndexPart1, IndexPart2, overlapDist] = mergeChannel(Doubles(~isnan(Doubles(n, :))), Position, Channel, IndexPart1, IndexPart2, overlapDist);
    [IndexPart2, IndexPart1, overlapDist] = mergeChannel(Doubles(~isnan(Doubles(n, :))), Position, Channel, IndexPart2, IndexPart1, overlapDist);
end

Singles = IndexPart1(sum(~isnan(IndexPart1), 2)>1, :); % up until now only the rows with several indexes are attributed to single probes

for n = 1:size(Singles, 1)
    Index = Singles(n, :);
    Index = Index(~isnan(Index));
    [~, b] = max(Intensities(Index)); % the row index of the maximum intensity dot. there is no need to check for index overlap because the dots come from different channels
    IndexPart1(logical(sum(eq(IndexPart1, Index(b)), 2)), :) = [Index(b) nan(1, size(IndexPart1, 2)-1)];
end








All = [IndexPart1(:); IndexPart2(:)];
[~, ~, c]=unique(All);
a_counts = accumarray(c,1);
if any(a_counts > 1)
    disp('overlap error')
end


%% From matrix turn into vector of indexes

IndexFinal = nan(1, size(IndexPart1, 1)); % the final vector is nan
IndexPart1 = IndexPart1'; % transpose since matlab reads through columns
IndexFinal(sum(~isnan(IndexPart1))>0) = IndexPart1(~isnan(IndexPart1)); % find no nan values in each column and add into the vector
end






function [IndexPart1, IndexPart2, varargout] = removeRepeat(IndexInCheck, IndexPart1, IndexPart2, Right, varargin)
Remove = find(IndexPart1 == IndexInCheck); % Matrix position of all dots with the same index of part1
Remove(eq(Remove, Right)) = []; % Remove the position from the right dot

IndexPart2(Remove) = NaN; % Matrix position of the repeated dot is removed in part2
IndexPart1(Remove) = NaN;  % Matrix position of the repeated dot is removed in part1
if ~isempty(varargin)
    varargin{1}(Remove) = NaN;
end

% The channel in question can be present in both Part1 and 2
Remove = find(IndexPart2 == IndexInCheck); % Matrix position of all dots with the same index of part1
Remove(eq(Remove, Right)) = []; % Remove the position from the right dot

IndexPart2(Remove) = NaN;
IndexPart1(Remove) = NaN;
if ~isempty(varargin)
    varargin{1}(Remove) = NaN;
end

if ~isempty(varargin)
    varargin{1}(logical(sum(eq(IndexPart1, IndexInCheck), 2)), :) = [varargin{1}(Right) nan(1, size(IndexPart1, 2)-1)];
end
IndexPart2(logical(sum(eq(IndexPart1, IndexInCheck), 2)), :) = [IndexPart2(Right) nan(1, size(IndexPart1, 2)-1)];
IndexPart1(logical(sum(eq(IndexPart1, IndexInCheck), 2)), :) = [IndexInCheck nan(1, size(IndexPart1, 2)-1)]; % The row represents one probe, so all others present in there are removed

varargout = varargin;
end





function [IndexPartInCheck, IndexPartOther, overlapDist] = mergeChannel(IndexInCheck, Position, Channel, IndexPartInCheck, IndexPartOther, overlapDist)
if ~isnan(IndexInCheck) && ~isempty(IndexInCheck)
    SameChannel = strcmp(Channel, Channel(IndexInCheck));
    Index = find(SameChannel);
    Dist = pdist2(Position(IndexInCheck, :), Position(SameChannel, :));
    RemoveCloseOnes = Index(Dist < 0.25 & Dist ~= 0);
   
    if ~isempty(RemoveCloseOnes)
        IndexPartOther(ismember(IndexPartInCheck, RemoveCloseOnes)) = NaN;
        overlapDist(ismember(IndexPartInCheck, RemoveCloseOnes)) = NaN;
        IndexPartInCheck(ismember(IndexPartInCheck, RemoveCloseOnes)) = NaN;
    end
end
end