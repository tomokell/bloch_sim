% This script runs the test_VEPCASL_seq_meanGz_angle for a variety of
% vessel angles and spin velocities
%
% Tom Okell, May 2023

%% Use default VEPCASL settings
set_VEPCASL_defaults_gaussianhann;
dt = dt/3;

%% Define the velocity and position vector of the spin, vessel positions and the
% minimum and maximum z position and offset of the labelling plane in the
% z direction
v = [5 20 40 80] / 100; %  cm/s to m/s conversion
ang = 0:5:85;  xory = 'x'; % Angles to be tested in xz plane

%% Define spin position to be tested
Ps = Pa; % The spin will pass through the "vessel A" position within the labelling plane

%% For each vessel angle and spin velocity, calculate the final Mz
FinalMz = []; IE = [];
for cycle = 1:4
  disp( ['Cycle ' num2str(cycle) ' of 4'] );
  
  for jj = 1:size(v,2)
    disp(['Velocity step ' num2str(jj) ' of ' num2str(size(v,2))])
  
    parfor ii = 1:size(ang,2)
      disp(['Angle step ' num2str(ii) ' of ' num2str(size(ang,2))])
    
      % Run the simulation
      [M, t, P, G, RF, T] = test_VEPCASL_seq_meanGz(v(jj), zmax, z_offset, meanGz, G_amp, RF_shape, ...
                                           RF_shape_params, RF_amp, RF_dur, ...
                                           RF_sep, Pa, Pb, cycle, Ps, dt, T1, ...
                                           T2, ORFreq, Unipolar, NegGz, ang(ii), xory)


      % Record the final z magnetisation
      FinalMz(ii,jj, cycle) = M(3,end);

      % Calculate the inversion efficiency
      %T_final = 2*zmax / (cos(ang(ii)*pi/180) * v(jj)); T_inv = T_final/2;
      IE(ii,jj,cycle) = InvEff(FinalMz(ii,jj,cycle),T, T/2, T1);
      
    end    
  end
end

%% Calculate the contrast
contr_ns      = IE(:,:,1) - IE(:,:,2); % Non-selective contrast
contr_vs_tag  = IE(:,:,3) - IE(:,:,2); % Vessel-selective contrast when tagged
contr_vs_cntl = IE(:,:,4) - IE(:,:,2); % Vessel-selective contrast when controlled

%% Plot IE for each cycle against angle for v = 20 cm/s
Idx = find(v==20/100);
figure
plot(ang, squeeze(IE(:,Idx,:)),'linewidth',2);
legend({'Tag All','Control All','Vessel Selective Tag', 'Vessel Selective Control'}, ...
	'Location','Best');
xlabel 'Vessel angle/degrees'
ylabel 'Inversion efficiency'
axis tight; ylim([0 1]);

%% Plot the contrast against angle for each velocity
figure
plot(ang, contr_vs_tag,'LineWidth',2);
legend(num2str(100*v'),'Location','Best');
xlabel 'Vessel angle/degrees'
ylabel 'Contrast'
axis tight; ylim([0 1]);