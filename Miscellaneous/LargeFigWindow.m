% Opens a large figure window with height and width given by the specified
% fractions of the screen size (defaults are 0.5 for both). Returns a
% handle to the new window.
%
% Tom Okell, June 2022
%
% function h = LargeFigWindow(HeightFrac, WidthFrac)

function h = LargeFigWindow(HeightFrac, WidthFrac)
  
  % Deal with optional arguments
  if nargin < 2; WidthFrac  = 0.5; end
  if nargin < 1; HeightFrac = 0.5; end
  
  % Get the screen size
  scrsz = get(0,'ScreenSize');
  
  % When the screen size call fails, this returns [1 1 1 1] so adapt here
  if (scrsz(3) == 1) || (scrsz(4) == 1)
      figure;
  else
      
      % Set up a figure window in the bottom left corner of the screen
      % figure('Position', [left bottom width height]);
      h = figure('Position',[1 1 scrsz(3)*WidthFrac scrsz(4)*HeightFrac]);
  end