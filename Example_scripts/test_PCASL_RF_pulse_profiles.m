% This script tests the effective flip angle at different z-locations due
% to various PCASL RF pulses
%
% Tom Okell, May 2023

%% Set up VEPCASL defaults
set_VEPCASL_defaults_gaussianhann;
cycle = 1; % Tag all
dt = 1e-6;

%% Determine the z positions of static spins to simulate
z_max = 5/100;
z = -z_max:(z_max/200):z_max;

% Set up simulation duration
T = RF_sep;
t = 0:dt:T;

% Set up the off resonance array
OffRes = ones(size(t)) * ORFreq;

%% Loop through RF pulse types
RF_shapes = {'hard', 'gaussianhann', 'hanning', 'hamming', 'gaussian_dsv'};
RF_shape_param_sets = {[], [4.0000e-04], [], [], []};

DesiredFlipAngle = 20;

% Array to save out the RF pulses
RF_pulses_for_plotting = [];

% Final effective flip angle array
Final_Ang = []; 

% Loop through RF pulse shapes and z locations
for ii = 1:length(RF_shapes)
    RF_shape = RF_shapes{ii};
    RF_shape_params = RF_shape_param_sets{ii};

    % Set up the sequence
    [G, RF] = make_VEPCASL_seq_meanGz_Gamp(RF_shape, RF_shape_params, RF_amp, ...
        RF_dur, RF_sep, dt, T, meanGz, G_amp, ...
        Pa, Pb, z_offset, cycle);

    % Modulate the RF amplitude to achieve the desired flip angle
    FlipAng = FlipAngle(RF,dt);
    RF(1,:) = RF(1,:) * DesiredFlipAngle/FlipAng;

    % Simulate each z position
    parfor zCount = 1:size(z,2)  % Loop through z positions
        disp(['zCount ' ns(zCount)])

        % Create a position array for this spin
        P = generate_position('static', [0 0 z(zCount)]', 0, t);

        % Run the Bloch simulation
        M = bloch_sim(P, G, RF, OffRes, t, T1, T2);

        % Record the final magnetisation angle
        Mfinal = M(:,end);
        Mxy = sqrt(Mfinal(1)^2 + Mfinal(2)^2);
        Mz = Mfinal(3);
        Final_Ang(ii,zCount) = atan2(Mxy, Mz) * 180 / pi;
    end

    RF_pulses_for_plotting(ii,:) = RF(1,:);
end

%% Plot the RF pulses
figure;
plot(t,RF_pulses_for_plotting', 'linewidth', 2);
xlabel 'time/s'
ylabel 'RF pulse amplitude/mT'
xlim([-10*dt RF_dur+10*dt])
legend(RF_shapes,'Interpreter','none')
title 'RF pulse shapes'

%% Plot the results
figure;
plot(z*100,Final_Ang, 'linewidth', 2);
set(gca, 'XTick',-5:5);
xlabel 'z position/cm'
ylabel 'Flip angle/degrees'
ylim([0 1.1*max(Final_Ang(:))]); xlim(100*[min(z) max(z)]);
legend(RF_shapes,'Interpreter','none')
title 'RF profiles'