function h = fillcircle( x, y, r, label, varargin )

p = inputParser;
addParameter(p,'FaceColor','none');
addParameter(p,'EdgeColor',[0 0 0]);
addParameter(p,'LineWidth',1);
addParameter(p,'TextColor',[0 0 0]);
addParameter(p,'TextSize',10);
parse(p,varargin{:});
params = p.Results;

h = fill( x+r*cos(-pi:0.01:pi), y+r*sin(-pi:0.01:pi), '', ...
      'FaceColor', params.FaceColor, ...
      'EdgeColor', params.EdgeColor, ...
      'LineWidth', params.LineWidth   );

% Add label is provided
text( x, y, label, ...
      'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
      'Color', params.TextColor, 'FontSize', params.TextSize, 'FontWeight', 'bold' );
