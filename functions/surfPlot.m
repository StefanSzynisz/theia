function surfPlot(X_coord,Y_coord,Data,plotName,zLabel, path, resolution,...
    left_pos,bott_pos,minValue, maxValue, colorBarSwitch, varargin)
% Usage: surfPlot(X_coord,Y_coord,Data,plotName,zLabel, path, resolution,...
%    left_pos,bott_pos,minValue, maxValue, colorBarSwitch, varargin)
% varagin are optional series for X- and Y-ticks, e.g.
% 0:0.2:1, 0:0.2:1

if (length(varargin)==2)
    Xticks = varargin{1};
    Yticks = varargin{2};
end

minX = min(X_coord); maxX = max(X_coord);
minY = min(Y_coord); maxY = max(Y_coord);

range_x = minX:resolution:maxX;
range_y = minY:resolution:maxY;

[A, B] = meshgrid(range_x,range_y);

C = griddata(X_coord,Y_coord,Data,A,B,'cubic');

width = 9;
height = 0.8 * width;
fig = get(groot,'CurrentFigure');  % Get the current figure handl
set(fig, 'Name',plotName, 'NumberTitle','on', 'Visible', 'on', ...
        'Units', 'centimeters', 'Position', [left_pos bott_pos width height]);

%fig.Renderer='Painters';
y = surf(A,B,C);
view(2); % is the top view of the surface
set(y, 'EdgeColor', 'none');
% y = contourf(A,B,C);
% Setting only x- and y- limits of the axis enables 2D top-view plot:
axis([minX maxX minY maxY]);
if minValue == maxValue   % if they are the same, then set the lowest to zero
   minValue = 0;
end
caxis([minValue maxValue]);  % Extent of values in the plot


if (length(varargin)==2)
    set(gca,'Xtick',Xticks)
    set(gca,'Ytick',Yticks)
end
colormap(jet);
if strcmp(colorBarSwitch,'colorbar_on')
    colorbar;
end

hTitle = title(plotName);

zlabel(zLabel);

set(gca,'FontSize',11,'FontName','Arial', 'fontweight','normal');
set( hTitle, ...
        'FontName'   , 'Helvetica',...
        'fontweight', 'bold',...
        'FontSize'   , 11);


fileName = strcat(plotName,'_Surf.png');
saveas(y, fullfile(path,fileName));
%saveas(y, fullfile(path,fileName),'epsc');
end
