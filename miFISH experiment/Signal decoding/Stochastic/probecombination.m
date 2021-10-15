function [Coordinates, Lamina] = probecombination(Dots)
%
% Given dots from one chromosome and nuclei, figure out how they should
% be mapped to a set of multi-colour probes.
%
% -> Input
% 1/ A list of dots, given as a matrix where each row is
% [c, x, y, z, i] 
% c: channel -- identified by a number
% x, y, z : cartesian coordinates of the fitted and shift correcte dots
% Please make sure that you multiply z with the correct value, 200 
% or 300 nm is typically used in the lab
% 2/ A list of probes simply saying what channels they have, i.e.,
% something like
% probe #1 is expected to have one dot in channel 1
% probe #2 is expected to have one dot in channel 1 and one in channel 2
% We encode that as
% P = {[1], [1, 2], [2, 3], ...}
% 
% -> Output
% A mapping, M saying something like
% probe #1 consists of dot #17
% probe #2 consists of dot #2 and #89
% ....
% M = {[17], [2, 89], ... }
% (i.e. probe centric)
%
% Depending on the objective function, optimization method and dots
% probes can have more or less than the expected number of dots as
% the objective function might "merge" two or mor close by dots into one.

% Name, short name, plot style
channels = {...
    {'a488', 'b', 'ko'}, ...
    {'TMR', 'g', 'rx'}, ...
    {'a594', 'o', 'g+'}, ...
    {'Cy5', 'r', 'b*'}, ...
    {'a700', 'p', 'cs'}, ...
    {'ir800', 'v', 'm>'}};
  

% Associate a numerical id to each channel
for cc = 1:numel(channels)
    p.(channels{cc}{2}) = cc;
end

% Use the channels ids to define the probes
Probes = {[p.g, p.o], [p.b], [p.v], [p.r, p.b], ...
    [p.g, p.p], [p.r], [p.p, p.o], [p.g, p.r], ...
    [p.p], [p.g, p.b], [p.r, p.o], [p.g], ...
    [p.b, p.o], [p.r, p.p], [p.o], [p.b, p.p]};


% fprintf("> Probes:\n");
% for pp = 1:numel(Probes)
%     probe = Probes{pp};
%     fprintf('Probe #%d: ', pp);
%     for cc = 1:numel(probe)
%         fprintf('%s (%s) ', ...
%             channels{probe(cc)}{1}, ...
%             channels{probe(cc)}{2});
%     end
%     fprintf('\n');
% end

fprintf('>> Finding mapping ...\n');    
[Mapping, score] = find_probecombination(Dots, Probes);

[Coordinates, Lamina] = plotDots(Dots, channels, Mapping);

% fprintf(">> Found the following mapping with score = %f\n", score);
% fprintf("> Probes:\n");
% for pp = 1:numel(Probes)
%     
%     probe = Probes{pp};
%     probedots = Mapping{pp};
%     
%     fprintf('Probe #%d: ', pp);
%     for cc = 1:numel(probe)
%         fprintf('%s (%s) ', ...
%             channels{probe(cc)}{1}, ...
%             channels{probe(cc)}{2});        
%     end
%     
%     fprintf('\n');
%     fprintf('\t %d dots: ', numel(probedots));    
%     
%     for mm = 1:numel(probedots)
%         dotid = probedots(mm);
%         dotchannel = Dots(dotid,1);
%         fprintf('#%d/%s ', dotid, channels{dotchannel}{1});
%     end
%     fprintf('\n');
% end

end

function [Coordinates, Lamina] = plotDots(D, Channels, Mapping)

% figure
% hold on
for kk = 1:size(D,1)
    X = D(kk, 2:4);
    cId = D(kk, 1);
    style = Channels{cId}{3};
    %plot3(X(1), X(2), X(3), style,'DisplayName',Channels{cId}{1});
end

Coordinates = nan(numel(Mapping), 3);
Lamina = nan(1, numel(Mapping));
if numel(Mapping) > 0
    for mm = 1:numel(Mapping)
        points = D(Mapping{mm}, 2:4);
        for aa = 1:size(points, 1)
            Lamina(mm) = D(Mapping{mm}(1), 6);
            Coordinates(mm, :) = points(1, :);
            for bb = aa+1:size(points, 1)
%                 plot3([points(aa, 1), points(bb,1)], ...
%                 [points(aa, 2), points(bb,2)], ...
%                 [points(aa, 3), points(bb,3)], 'k');
                Coordinates(mm, :) = (points(1, :)+points(2, :))/2;
            end
        end
    end    
end


% view(3)
% axis equal
% grid on

end

function [M, score] = find_probecombination(D, P)
%% function [M, s] = find_probecombination(D, P)
% Find the the best mapping, M, reported by it's score
% given the dots in D and probes in P
% see probecombination for the format of the input and ouput data


% DM: dot centric mapping. DM(kk) says which
% probe that dot kk belongs to

nrestarts = 20; niter = 100;

best_score = inf; % best score, the lower the better
best_M = cell([numel(P),1]); 

for rr = 1:nrestarts
M0 = cell([numel(P),1]);  % Zero guess, no dots assigned
M = cell([numel(P),1]); 
score = inf;

for kk = 1:niter     
    score0 = objfun(D, P, M0);    
    if score0 < score
        %keyboard
        % save best
        score = score0; 
        M = M0; 
        % Proceed
    else
        % Revert one step
        M0 = M;
    end
    M0 = generateRandomProbe(D, P, M0);
end
if(score < best_score)
    best_score = score;
    best_M = M;
end
end

score = best_score;

% convert dot mapping to probe centric format
M = best_M;

end

function p = objfun(Dots, Probes, ProbeDotIDs)
% Objective function to MINIMIZE
%
% Dots: All dots
% Probe: Cell structure with one array per probe saying what channels
%        the dots are expected to come from
% DotMap: Says which probe each dots belongs to (to be evaluated here)

maxEl=2;%max number of elements of a cell element in your cell
DotMap = cell2mat(cellfun(@(x) [x nan(1,maxEl-numel(x))],ProbeDotIDs, 'uni', 0));

if isempty(ProbeDotIDs{1})
    p = inf;
    return
end


%% Errors

% Distance between dots in the same probe
e_probe_intra_dist = 0;
for pp = 1:numel(Probes) % Per probe
    dots_xyz = Dots(ProbeDotIDs{pp}, 2:4);
    e_probe_intra_dist = e_probe_intra_dist + sum(pdist(dots_xyz));
end

% Repeated probes
[~, ~, c]=unique(DotMap);
a_counts = accumarray(c,1);
e_repeat_dots = sum(a_counts > 1);


% 1 = the top brighest are used, 0 = the bottom dimmest are used
e_brightness = 0;
for pp = 1:numel(Probes) % Per probe
    bright_selected = Dots(ProbeDotIDs{pp}, 5);
    e_brightness = e_brightness + sum(1-bright_selected);
end

% penalise missing dots
e_ndot = sum(sum(isnan(DotMap), 2)>1);

% closeness of dots from same channel
e_merge_dot = 0;
idx = DotMap(~isnan(DotMap));
for kk = 1:6 % Per channel
    Dist = pdist(Dots(idx(Dots(idx, 1) == kk), 2:4));
    if any(Dist < 0.25)
        e_merge_dot = e_merge_dot + 1;
    end
end

% Finally: Define what is important by weighting together the terms
% brightness error is from 0 to 18
% repeat dots error is from 0 to 7
% overlap distance error is from 0 to 8 -> most important
p = e_brightness + 3*e_repeat_dots + 10*e_probe_intra_dist + 10*e_ndot + 10*e_merge_dot; 

assert(isfinite(p))
assert(~isnan(p))
    
end



function M0 = generateRandomProbe(D, P, M0)
for probe = 1:numel(P)
    for n = 1:10
        dot1idx = find(ismember(D(:, 1),P{probe}));
        Dtry = randi(length(dot1idx)); % pick a random dot from a valid group
        Dtry = dot1idx(Dtry);
        Dtry2 = [];
        
        if numel(P{probe}) > 1 % for overlaping probes
            dot2 = P{probe}(~ismember(P{probe}, D(Dtry)));
            dot2idx = find(D(:, 1)==dot2);
            [minVal, idx] = min(pdist2(D(Dtry, 2:4), D(dot2idx, 2:4))); % find the closest dot
            if minVal < 0.55
                Dtry2 = dot2idx(idx); % assign if are located below 550nm
                break
            else
                Dtry = [];
                Dtry2 = [];
            end
        else
            break
        end
    end
    M0{probe} = [Dtry Dtry2]; % One probe assigned
end
end