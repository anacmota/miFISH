function [Contact, ContactRatio] = plot3D_splineCurve(Coordinates, OrderPairs, GenomicDistance, varargin)

ProbesPosition = readtable('probes_position.csv');
Index = floor(ProbesPosition.chromstart/10^5);
fileID = fopen('eigen_KR_2_100kb.txt','r');
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