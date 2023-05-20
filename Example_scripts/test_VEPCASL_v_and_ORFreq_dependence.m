% This script runs the test_VEPCASL_seq_meanGz function for a variety of
% velocities and off-resonance frequencies
%
% Tom Okell, May 2023

%% Set the VEPCASL defaults
set_VEPCASL_defaults_gaussianhann;

%% Define the off-resonance frequencies to be tested
ORFreq = -200:10:200;

%% Define velocities to be tested
v = [3 10 30 60 100]/100; % cm/s to m/s conversion

%% Initialise
FinalMz = zeros(2,2,size(v,2)); IE = zeros(2,2,size(v,2));

%% For each off-resonance frequency and velocity, calculate IE
for ORCount = 1:size(ORFreq,2)
  for cycle = 1:2

    % Loop through a number of velocity values, calculating the final Mz in
    % each case
    parfor VCount=1:size(v,2)
  
      disp(['Simulation ' num2str(VCount + (cycle-1)*size(v,2) + (ORCount-1)*2*size(v,2)) ...
		  ' of ' num2str(size(ORFreq,2)*2*size(v,2))])
      % Run the simulation
      [M, t, P, G, RF] = test_VEPCASL_seq_meanGz(v(VCount), zmax, z_offset, meanGz, G_amp, RF_shape, ...
                                                  RF_shape_params, RF_amp, RF_dur, ...
                                                  RF_sep, Pa, Pb, cycle, ...
                                                  Ps, dt, T1, T2, ORFreq(ORCount));

      % Record the final z magnetisation
      FinalMz(cycle,ORCount,VCount) = M(3,end);

    end
    
    % Calculate the inversion efficiency
    IE(cycle,ORCount,:) = reshape(InvEff(squeeze(FinalMz(cycle,ORCount,:)), ...
                                     2*zmax./v', zmax./v', T1), 1,1,size(v,2));
    
  end    
end

%% Calculate the contrast
contr = squeeze(IE(1,:,:) - IE(2,:,:));

%% Repeat for IE tag and control
figure
plot(ORFreq', squeeze(IE(1,:,:)), 'linewidth', 2);
hold on
plot(ORFreq', squeeze(IE(2,:,:)), '-.', 'linewidth', 2);
xlabel 'Off resonance frequency/Hz'
ylabel 'Inversion efficiency'
ylim([0 1])
legend(num2str(v'*100),'Location',[0.45 0.2 0.15 0.3])

%% Plot the tagging efficiency as a function of velocity
figure
plot(ORFreq', contr, 'linewidth', 2);
xlabel 'Off resonance frequency/Hz'
ylabel 'Contrast'
ylim([0 1])
legend(num2str(v'*100),'Location','South')

% NB. Note there is some asymmetry in these plots, especially at low
% velocities, which arises from the effective labelling plane position
% shifting in the z direction with off-resonance, so when we correct for T1
% decay back to when the spin crossed the nominal labelling plane location,
% this causes some asymmetry depending on whether it is a positive or
% negative off-resonance frequency.