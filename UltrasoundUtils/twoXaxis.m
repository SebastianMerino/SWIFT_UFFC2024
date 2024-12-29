function [ax1,ax2]=twoXaxis(x1,y1,x2,y2,er1,er2,color1,color2)
if nargin == 0
    % Demo
    x1 = 0:0.1:10;
    y1 = sin(x1);
    color1 = [0 0.4470 0.7410];
    
    x2 = 0:0.2:20;
    y2 = cos(x2);
    color2 = [0.8500 0.3250 0.0980];
elseif nargin==4 || nargin==6
    % Default colors
    color1 = [0 0.4470 0.7410];
    color2 = [0.8500 0.3250 0.0980];
end

% Create the first axes
ax1 = axes;
ax1.Position = [0.13 0.2 0.775 0.6];
plot(x1, y1, '-', 'LineWidth', 1.5, 'Color',color1);
if nargin>4
    hold on, errorbar(x1, y1, er1, '.', 'Color',color1, 'CapSize',2);
    hold off
end
xlabel(ax1, 'X-axis 1');
ylabel('Y-axis');
ax1.XAxis.Color = color1;
ax1.Box = 'off';
ax1.GridColor = [0.15,0.15,0.15];
grid on;

% Create the second axes (overlayed on the first)
ax2 = axes;
ax2.Position = [0.13 0.2 0.775 0.6];
plot(x2, y2, '--', 'LineWidth', 1.5, 'Color',color2); % Second plot
if nargin>4
    hold on, errorbar(x2, y2, er2, '.', 'Color',color2, 'CapSize',2);
    hold off
end
ax2.Color = 'none';         % Make the background transparent
ax2.XAxisLocation = 'top';  % Position the second x-axis on top
ax2.YAxis.Visible = 'off';  % Hide the second y-axis
ax2.XAxis.Color = color2;
ax2.Box = 'off';
xlabel(ax2, 'X-axis 2');


end