%% Visualise CASL and PCASL inversion processes
%
% Tom Okell, May 2023

%% Set simulation parameters
set_VEPCASL_defaults_gaussianhann;
v = 10/100; % Velocity cm/s -> m/s
T = 1000e-3; % Total simulation time in s
zmax = v*T/2; % Maximum z location to simulate
% NB. this ensures the spin crosses the labelling plane halfway through the simulation
dt = 10e-6; % Time increment for the simulation 
% (needs to be small enough to well characterise RF and gradient pulses)

%% Simulate inversion due to CASL
G_amp = -0.8; % mT/m. 
% NB. use a negative gradient amplitude here to aid visualisation, but a
% positive value works just as well
RF_amp = 0.0014; % mT
[Mc, tc, Pc, Gc, RFc] = test_CASL_seq(v, zmax, G_amp, RF_amp, ORFreq, dt, T1,T2,false,true);

%% Set PCASL sequence parameters
set_VEPCASL_defaults_gaussianhann;
v = 10/100; T = 1000e-3; zmax = v*T/2; dt = 10e-6;
cycle = 1; % PCASL tag

%% Simulate inversion due to PCASL
% As above, use negative Gz during RF pulses as above to show the inversion more
% clearly in this view, although the process is identical with postive Gz
% too.
NegGz = true; 
[Mp, tp, Pp, Gp, RFp] = test_VEPCASL_seq_meanGz(v, zmax, z_offset, meanGz, G_amp, RF_shape, ...
                                                RF_shape_params, RF_amp, RF_dur, ...
                                                RF_sep, Pa, Pb, cycle, Ps, dt, T1, ...
                                                T2, ORFreq, Unipolar, NegGz);

%% Plot x, y and z components separately
figure; subplot(1,2,1);
plot(tc*1000,Mc'); legend('x','y','z'); xlabel 'Time/ms'; ylabel 'Magnetisation/M_0'
xlim([min(tc) max(tc)]*1000); ylim([-1 1]); grid on; title 'CASL'

subplot(1,2,2); plot(tp*1000,Mp'); legend('x','y','z'); xlabel 'Time/ms'; ylabel 'Magnetisation/M_0'
xlim([min(tp) max(tp)]*1000); ylim([-1 1]); grid on; title 'PCASL'

%% Plot 3D versions
% Find times for plots
Idx(1) = find(tp > 450e-3,1); % With the spin before the labelling plane
Idx(2) = find(tp > 500e-3,1); % Spin at the centre of the labelling plane
Idx(3) = find(tp > 550e-3,1); % Spin beyond the labelling plane

LargeFigWindow(1,0.3);
% CASL
for ii = 1:3
  pIdx = (ii-1)*4+1;
  % NB. we use 7 rows here so each subplot is double the height of the
  % legend added at the bottom.
  subplot(7,2,[pIdx pIdx+2])
  ThisIdx = 1:5:Idx(ii);
  vis_bloch_sim_3D_with_Beff(tc(ThisIdx),Mc(:,ThisIdx),Pc(:,ThisIdx),Gc(:,ThisIdx),RFc(:,ThisIdx),zeros(size(ThisIdx)),true,false,false,true);  
  if ii == 1; title 'CASL'; end

end

% Repeat for PCASL
for ii = 1:3
  pIdx = 4*ii - 2;
  subplot(7,2,[pIdx pIdx+2])
  ThisIdx = 1:5:Idx(ii);
  vis_bloch_sim_3D_with_Beff(tp(ThisIdx),Mp(:,ThisIdx),Pp(:,ThisIdx),Gp(:,ThisIdx),RFp(:,ThisIdx),zeros(size(ThisIdx)),true,false,false,false);  
  if ii == 1; title 'PCASL'; end
end

% Create an additional dummy plot with the legend, but make the axis small so
% it's not visible and only the legend shows
hSubplot = subplot(7,2,13:14); 
vis_bloch_sim_3D_with_Beff(tp(ThisIdx),Mp(:,ThisIdx),Pp(:,ThisIdx),Gp(:,ThisIdx),RFp(:,ThisIdx),zeros(size(ThisIdx)),true,true,false,true);  
pos = get(gca,'position');
set(gca,'position',[(pos(1)+pos(3)/2) (pos(2)+pos(4)*3/4) 0.01 0.01])
hLegend = findobj(gcf, 'Type', 'Legend'); axis off;
hLegend.Location = 'best'; hLegend.Orientation = 'horizontal';

% Note that in the legend we have:
% - M_i: initial magnetisation
% - M_f: final magnetisation
% - M_path: the path of the tip of the magnetisation
% - B_G: Gradient component of the effective field (G.r)
% - B_1: RF component of the effective field
% - B_OR: Off-resonance component of the effective field
% - B_eff: Total effective field, about which the magnetisation precesses
% 
% These effective field components are only plotted for CASL, since they
% vary in a more complex way during the PCASL pulse sequence, making them
% harder to interpret in a few frames.