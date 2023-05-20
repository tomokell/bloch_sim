% This script tests the perturbation of the magnetisation near the PCASL
% labelling plane in label and control modes.
%
% Tom Okell, May 2023

%% Set up VEPCASL defaults
set_VEPCASL_defaults_gaussianhann;

%% Set up simulation duration
T = 0.5; % s
t = 0:dt:T;

%% Set up the off resonance array
OffRes = ones(size(t)) * ORFreq;

%% Determine the z positions of static spins to simulate

% Calculate the position estimated to give the greatest saturation
z_sat = pi / (g * meanGz  * RF_sep); % in m

% Determine z positions to test
z = (-z_sat*2.2):(z_sat/30):(z_sat*2.2);

%% Run simulations for the tag and control all cycles
Mz_final = []; SimNo = 1;
for cycle = 1:2  % Loop through tag and control cycles

    for meanGzCount = 1:size(meanGz,2)  % Loop through mean Gz values

        % Generate the VEPCASL sequence
        [G, RF] = make_VEPCASL_seq_meanGz_Gamp(RF_shape, RF_shape_params, RF_amp, ...
    	    RF_dur, RF_sep, dt, T, meanGz(meanGzCount), G_amp, Pa, Pb, z_offset, cycle);

        parfor zCount = 1:size(z,2)  % Loop through z positions
            disp(['cycle ' ns(cycle) ', MeanGz ', ns(meanGzCount) ', zCount ' ns(zCount)])

            % Create a position array for this spin
            P = generate_position('static', [0 0 z(zCount)]', 0, t);

            % Run the Bloch simulation
            M = bloch_sim(P, G, RF, OffRes, t, T1, T2);

            % Record the final magnetisation
            Mz_final(zCount, meanGzCount, cycle) = M(3,end);
    	end
    end
end


%% Plot the results against each other
% Tag Condition
figure;
plot(z*100,Mz_final(:,:,1),'linewidth', 2);
xlabel 'z position/cm'
ylabel 'Final z magnetisation'
ylim([0 1.2]); xlim(100*[min(z) max(z)]); 
title 'Label condition'

% Control condition
figure;
plot(z*100,Mz_final(:,:,2),'linewidth', 2);
xlabel 'z position/cm'
ylabel 'Final z magnetisation'
ylim([0 1.2]); xlim(100*[min(z) max(z)])
title 'Control condition'