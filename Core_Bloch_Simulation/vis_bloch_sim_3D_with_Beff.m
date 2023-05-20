% Visualise a bloch simulation in 3D with the effective B field
%
% Tom Okell, May 2023
%
% h = vis_bloch_sim_3D_with_Beff(t, M, P, G, RF, OffRes, LHR, LegendOn, NewWin, PlotBeff)

function h = vis_bloch_sim_3D_with_Beff(t, M, P, G, RF, OffRes, LHR, LegendOn, NewWin, PlotBeff)

  if nargin < 7; LHR = true, end
  if nargin < 8; LegendOn = true, end
  if nargin < 9; NewWin = true, end
  if nargin < 10; PlotBeff = true, end
  
  % If using left hand rule for rotations, modify the RF phase to be
  % negative.
  if LHR; RF(2,:) = - RF(2,:); end
  
  % Create a new figure window if necessary
  if NewWin; h = figure; end;
  hold on
  
  % Set variables for arrows
  stemRatio = 0.8; LineWidth = 0.015; ArrowHeadWidth = 0.05;
  
  % Plot the initial magnetisation
  %quiver3(0,0,0,M(1,1),M(2,1),M(3,1),0,'Color',[0.8,0.5,0],'MaxHeadSize',5
  %00);
  arrow3D([0,0,0],M(:,1),  [0.8,0.5,0],stemRatio, LineWidth, ArrowHeadWidth); 
  
  % Plot the final magnetisation
  arrow3D([0,0,0],M(:,end),[0.3,0,0.8],stemRatio, LineWidth, ArrowHeadWidth);
  
  % Plot a line linking the magnetisation vector tips
  plot3(M(1,:)', M(2,:)', M(3,:)', 'c-', 'linewidth', 1); 
  
  % Calculate Beff and its components at the final time point and plot
  if PlotBeff
      g = GetGamma; Beff_scaling = 150;
      
      % Gradient component
      BG = [0,0,P(:,end)'*G(:,end)];
      arrow3D([0,0,0],BG*Beff_scaling,[0,0.8,0],stemRatio, LineWidth, ArrowHeadWidth);
      
      % RF component
      B1 = [ RF(1,end)*cos( RF(2,end) ); RF(1,end)*sin( RF(2,end) ); 0];
      arrow3D([0,0,0],B1*Beff_scaling,[0.8,0,0],stemRatio, LineWidth, ArrowHeadWidth);
      
      % Off resonance component
      BOR = [0,0,2*pi*OffRes(end)/g];
      arrow3D([0,0,0],BOR*Beff_scaling,[0,0,0.8],stemRatio, LineWidth, ArrowHeadWidth);
      
      % Total Beff
      Beff = Calc_Beff( P(:,end), G(:,end), RF(:,end), OffRes(:,end), g );
      arrow3D([0,0,0],Beff*Beff_scaling,[0.8,0.8,0],stemRatio, LineWidth, ArrowHeadWidth);
  end
  
  % Show grid and have equal length axes greater than one to help visualisation
  grid on
  axis square
  axis([-1 1 -1 1 -1 1]*1); 
  
  % Label the axes
  xlabel('M_{x''}');
  ylabel('M_{y''}');
  zlabel('M_{z''}');
    
  % View from standard orientation before using quiver to ensure arrows can
  % be seen clearly and set lighting
  view(45,15); %lighting phong; camlight head;
  
  % Legend
  if LegendOn
    legend('M_i','M_f','M_{path}','B_G','B_1','B_{OR}','B_{eff}')
  end
  
 